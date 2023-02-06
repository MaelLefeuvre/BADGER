use std::{str, env, io::{self, BufRead, Write}};

mod error;
use error::{HeaderError, GenotypeError, MainError};

const GENOTYPES_IDX          : usize = 9;
//const DEFAULT_BUFFER_CAPACITY: usize = 2 * 1024 * 10usize.pow(6); // 1GB MAX RSS
const DEFAULT_BUFFER_CAPACITY: usize = 8192 ; // 8Kb
const HEADER_PEEK_LEN        : usize = 6;

/// Read standard input until we found the last header line of the VCF (i.e.: the line starts with '#CHROM').  
/// Returns the field indices of the samples that are requested for duplication.
/// 
/// # Side effects
/// Since we want to duplicate sample columns, the last header line where samples are declared is modified
/// accordingly: duplicated samples are thus appended to the end of the header line. To distinguish these
/// entries from their original values, the sample name is modified according to the provided `tamper` 
/// closure or function.
fn tamper_header<W, R, F>(stdin: &mut R, stdout: &mut W, tamper: F, requested_samples: &[String]) -> Result<Vec<usize>, HeaderError>
where   R: BufRead,
        W: Write,
        F: FnMut(&str) -> String
{
    use HeaderError::*;
    let mut peeker: [u8; HEADER_PEEK_LEN] = [0,0,0,0,0,0];
    let mut buffer = Vec::with_capacity(DEFAULT_BUFFER_CAPACITY);
    loop {
        // ---- Check the first n bytes of each line, and check if it starts with
        //      '#CHROM'. Leave each line untampered and flush them to stdout
        //      until this is the case.
        stdin.by_ref().read_exact(&mut peeker).map_err(PeekHeader)?;
        buffer.extend(&peeker);
        stdin.read_until(b'\n', &mut buffer).map_err(ReadStdin)?;
        if  &peeker == b"#CHROM" { // We've arrived at destination!
            buffer.pop(); // Remove '\n'

            // ---- Parse the samples id of the header.
            let samples: Vec<String> = str::from_utf8(&buffer)
                .unwrap()
                .split('\t')
                .skip(GENOTYPES_IDX)
                .map(String::from)
                .collect();
            
            // ---- Match our requested_samples id with those in the header
            //      to obtain the field indices (shifted by -GENOTYPES_IDX)
            let indices = requested_samples.iter()
                .flat_map(|req_sample| samples
                    .iter()
                    .position(|vcf_sample| vcf_sample == req_sample)
                )
                .collect::<Vec<usize>>();
            
            // ---- Modify the requested's sample ID using our tamper function.
            //      This is so we can distinguish 'dopplegangers' from their original
            //      counterpart.
            samples.iter()
                .enumerate()
                .filter(|(i, _)| indices.contains(i))
                .map(|(_, sample)| sample.as_str())
                .map(tamper)
                .for_each(|sample| { buffer.push(b'\t'); buffer.extend(sample.as_bytes()) });
            
            // ---- Add a final '\n', flush the buffer a final time and leave.
            buffer.push(b'\n');
            stdout.write_all(&buffer).map_err(WriteStdout)?;
            return Ok(indices)
        }

        // ---- No luck... Flush the buffer to stdout and carry on to the next line.
//        if buffer.len() >= DEFAULT_BUFFER_CAPACITY {
            stdout.write_all(&buffer).map_err(WriteStdout)?;
            buffer.clear();
//        }

    }
}

/// Read standard input and duplicate the genotypes fields at the requested indices.
/// 
/// Genotype fields are duplicated in their order of appearance.
fn duplicate_genotypes<R, W>(stdin: &mut R, stdout: &mut W, indices: &[usize]) -> Result<(), GenotypeError>
where   R: BufRead,
        W: Write
{
    use GenotypeError::*;
    let mut buffer            = Vec::with_capacity(DEFAULT_BUFFER_CAPACITY);
    let mut sample_buffer     = [0,0,0,0];
    let mut duplicate_samples = Vec::with_capacity(indices.len());

    while stdin.fill_buf().map(|b| !b.is_empty()).map_err(ReadLine)? {
        // ---- Skip fields CHROM -> FORMAT, and fill buffer w/ it.
        for _ in 0..GENOTYPES_IDX { 
            stdin.read_until(b'\t', &mut buffer).map_err(ReadInfo)?;
        };
        
        // ---- Sequentially read genotype fields within that line.
        let mut sample_counter = 0;
        'genotypes:  loop {
            stdin.read_exact(&mut sample_buffer).map_err(ReadGenotype)?;

            // ---- If that column was requested for duplication, store the 
            //      alleles in a separate buffer
            if indices.contains(&sample_counter) {
                duplicate_samples.push(sample_buffer);
            }
            sample_counter+=1;

            buffer.extend(sample_buffer); // Flush alleles to stdout.

            // ---- If we're at the end of the line, it's time to append our
            //      duplicated columns to our buffer and flush its contents 
            //      to stdout.
            if sample_buffer.ends_with(&[b'\n']) {
                *buffer.last_mut().unwrap() = b'\t';
                for sample in duplicate_samples.iter() {
                    buffer.extend(&sample[0..3]);
                    buffer.push(b'\t');
                }
                *buffer.last_mut()
                    .expect("Unexpectedly empty genotype buffer") = b'\n';
                
                break 'genotypes
            }
        }
        duplicate_samples.clear();
//        if buffer.len() >= DEFAULT_BUFFER_CAPACITY {
            stdout.write_all(&buffer).map_err(WriteStdout)?;
            buffer.clear();
//        }
    }

    // ---- Flush one last time before exiting.
    stdout.write_all(&buffer).map_err(WriteStdout)?;
    buffer.clear();
    Ok(())
}

/// Function to manipulate the duplicate sample name, and thus distinguish them within the VCF.
fn tamper_sample_name(name: &str) -> String {
    name.to_ascii_uppercase().replace("PED", "ped")
}

/// Read a VCF file from standard input, and duplicate any sample column which was provided as 
/// positional arguments
/// 
/// # Usage:
/// bcftools view pedigree.vcf.gz | duplicate-vcf-samples ped1_g2-b1-i1 ped2_g2-b1-i1 ped3_g2-b1-i1 | bgzip > pedigree-twins.vcf.gz
/// 
/// Columns are duplicated @ the end of the vcf's, in order of appearance. sample names are modified using
/// the [`tamper_sample_name`] function to distinguish them from their original.
/// 
/// # @TODO:
/// 1. While this has been thoroughly tested through command line comparison w/ the output of bcftools merge,
///    There isn't proper unit testing  ensuring allele information is properly duplicated.
/// 2. Flushing buffers at each line is highly inefficient. Should instead flush when buffer 
///    size >= DEFAULT_BUFFER_CAPACITY
fn main() -> Result<(), MainError> {
    let requested_samples = env::args().collect::<Vec<String>>();
    let mut stdin         = io::stdin().lock();
    let mut stdout        = io::stdout().lock();

    tamper_header(&mut stdin, &mut stdout, tamper_sample_name, &requested_samples)
        .map(|indices| duplicate_genotypes(&mut stdin, &mut stdout, &indices) )??;

    Ok(())
}

#[cfg(test)]
mod tests {

    const HEADER: &str = r###"##contig=<ID=1,assembly=b37,length=249250621>
##contig=<ID=5,assembly=b37,length=180915260>
##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant the line represents">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ped1_g1-b1-i1	ped1_g1-b2-i1	ped1_g2-b1-i1
1	752721	.	A	G	100	PASS	NS=3;AA=.|||;VT=SNP;	GT	0|1	1|1	0|1
5	654645	.	C	T	100	PASS	NS=3;AA=.|||;VT=SNP;	GT	1|1	1|0	0|1"###;

    use std::io::Cursor;
    use super::*;

    fn mock_stream() -> (Cursor<Vec<u8>>, Vec<usize>) {
        let mut mock_stdin  = Cursor::new(HEADER);
        let mut mock_stdout = Cursor::new(vec![]);
        let request_samples = vec!["ped1_g2-b1-i1", "ped1_g1-b2-i1"]
            .into_iter()
            .map(String::from)
            .collect::<Vec<String>>();
        let output = tamper_header(&mut mock_stdin, &mut mock_stdout, tamper_sample_name, &request_samples,)
            .expect("Failed to flush header");
        (mock_stdout, output)
    }

    #[test]
    fn flush_header() {
        let (header, _) = mock_stream();
        let bytes       = header.into_inner(); 
        // Ensure there are four newlines (i.e. the number of header lines.)
        assert_eq!(bytes.iter().filter(|byte| **byte == b'\n').count(), 4);
    }

    #[test]
    fn fetch_samples_indices() {
        let (_, samples) = mock_stream();
        // Ensure there are two samples in the output
        assert_eq!(samples, vec![2, 1]);
    }

    #[test]
    fn check_duplicated_header() {
        let (header, _) = mock_stream();
        let header = header.into_inner();
        // Ensure there are two samples in the output
        let last_header_line = str::from_utf8(&header)
            .expect("Invalid UTF-8 in header")
            .split('\n')
            .nth_back(1)
            .expect("Empty header");

        let header_samples = &last_header_line
            .split('\t')
            .collect::<Vec<&str>>()[GENOTYPES_IDX..];

        let wanted_samples = vec![
            "ped1_g1-b1-i1", "ped1_g1-b2-i1", "ped1_g2-b1-i1", "ped1_G1-B2-I1", "ped1_G2-B1-I1"
        ];
        assert_eq!(header_samples, wanted_samples);
    }
}
