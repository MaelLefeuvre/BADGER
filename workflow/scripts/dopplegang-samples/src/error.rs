use thiserror::Error;

use super::HEADER_PEEK_LEN;

#[derive(Debug, Error)]
pub enum HeaderError {
    #[error("Failed to peek the first {} of the vcf header [{0}]", HEADER_PEEK_LEN)]
    PeekHeader(#[source] std::io::Error),

    #[error("Failed to read header line: [{0}]")]
    ReadStdin(#[source] std::io::Error),

    #[error("Failed to flush header contents to stdout")]
    WriteStdout(#[source] std::io::Error)
}

#[derive(Debug, Error)]
pub enum GenotypeError {
    #[error("Failed to read genotype line: [{0}]")]
    ReadInfo(#[source] std::io::Error),

    #[error("Failed to read genotype line: [{0}]")]
    ReadGenotype(#[source] std::io::Error),

    #[error("Failed to read the next line of the VCF")]
    ReadLine(#[source] std::io::Error),

    #[error("Failed to flush genotypes contents to stdout")]
    WriteStdout(#[source] std::io::Error)
}

#[derive(Debug, Error)]
pub enum MainError {
    #[error("Failed to parse VCF header: [{0}]")]
    TamperHeader(#[from] HeaderError),

    #[error("Failed to duplicate genotypes of sample(s): [{0}]")]
    DuplicateGenotypes(#[from] GenotypeError)
}