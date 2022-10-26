localrules: generate_bam_list


def get_samples_ids(wildcards):
    # Run through the initial samples files and extract pedigree ids 
    with checkpoints.get_samples.get().output[0].open() as f:
        samples = str.split(f.readline().replace('\n', ''), '\t')
        ids     = set([sample.split('_')[1] for sample in samples])
        return sorted(expand("{generation}_{ids}", ids=ids, generation="{generation}"), key=str.casefold)


def get_pileup_input_bams(wildcards):
    # Run through the initial samples files and extract pedigree ids 
    samples_ids = get_samples_ids(wildcards)

    # Return a list of input bam files for pileup
    return expand(
        "results/02-preprocess/06-mapdamage/{sample_id}/{sample_id}.srt.rmdup.rescaled.bam",
        sample_id = samples_ids,
    )


rule generate_bam_list:
    """
    Create an ordered list of bam files using the original samples file
    """
    input:
        bamlist = get_pileup_input_bams
    output:
        bamlist = "results/04-kinship/READ/{generation}/{generation}.bam.list"
    priority: 15
    shell: """
        ls {input.bamlist} > {output.bamlist}
    """

def get_samtools_optargs(wildcards):
    optargs = ""
    if config['variant-calling']['pileup']['disable-BAQ']:
        optargs += "-B "


rule samtools_pileup:
    input:
        bamlist   = rules.generate_bam_list.output.bamlist,
        targets   = os.path.splitext(config["kinship"]["targets"])[0] + ".ucscbed",
        reference = config["reference"],
        fai       = config["reference"] + ".fai",
    output:
        pileup    = "results/03-variant-calling/01-pileup/{generation}/{generation}.pileup"
    params:
        min_MQ    = config['variant-calling']['pileup']['min-MQ'],
        min_BQ    = config['variant-calling']['pileup']['min-BQ'],
        optargs   = get_samtools_optargs
    log:   "logs/03-variant-calling/{generation}-samtools_pileup.log"
    conda: "../envs/samtools-1.15.yml"
    shell: """
        samtools mpileup \
           -R {params.optargs}\
           -q {params.min_MQ} \
           -Q {params.min_BQ} \
           -l {input.targets} \
           -f {input.reference} \
           -b {input.bamlist} | LC_ALL=C sort -k1,1n -k2,2n > {output.pileup} 2> {log} \
        """

# ----------------------------------------------------------------------------------------------- #
#                               Random haploid calling with PileupCaller                          #
# ----------------------------------------------------------------------------------------------- #

def parse_pileup_caller_flags(wildcards):
    args=""
    if config['variant-calling']['pileupCaller']["skip-transitions"]:
        args+="--skipTransitions "

    seed = config['variant-calling']['pileupCaller']['seed']
    if seed is not None:
        args += f"--seed {seed} "

    match config['variant-calling']['pileupCaller']["mode"]:
        case "randomHaploid":
            args += "--randomHaploid "
        case "majorityCall":
            args += "--majorityCall "
        case _:
            raise RuntimeError("Incorrect pileupCaller mode selected.")
    return args


rule pileup_caller:
    input:
        pileup            = rules.samtools_pileup.output.pileup,
        bamlist           = rules.generate_bam_list.output.bamlist,
        targets           = os.path.splitext(config["kinship"]["targets"])[0] + ".snp",

    output:
        plink             = multiext("results/03-variant-calling/02-pileupCaller/{generation}/{generation}", ".bed", ".bim", ".fam"),
        sample_names_file = "results/03-variant-calling/02-pileupCaller/{generation}/{generation}-names.txt"
    params:
        basename          = "results/03-variant-calling/02-pileupCaller/{generation}/{generation}",
        optargs           = parse_pileup_caller_flags,
        min_depth         = config['variant-calling']["pileupCaller"]["min-depth"],
        seed              = config['variant-calling']['pileupCaller']["seed"],
        sample_names      = lambda wildcards: expand(get_samples_ids(wildcards), generation = wildcards.generation)
    log:   "logs/03-variant-calling/{generation}-pileup_caller.log"
    conda: "../envs/sequencetools-1.5.2.yml"
    shell: """
        echo {params.sample_names} | tr ' ' '\n' > {output.sample_names_file}
        pileupCaller {params.optargs} \
            --snpFile <( awk '($2<=22)' {input.targets}) \
            --sampleNameFile {output.sample_names_file} \
            --plinkOut {params.basename} \
            --minDepth {params.min_depth} \
            --samplePopName {wildcards.generation} < {input.pileup} > {log} 2>&1 
    """


rule plink_bfile_to_tped:
    input:
        bfile  = multiext("{directory}/{file}", ".bed", ".bim", ".fam"),
    output:
        tplink = multiext("{directory}/{file}", ".tped", ".tfam") 
    params:
        basename = "{directory}/{file}"
    log:   "logs/03-variant-calling/{directory}/{file}-plink_bfile_to_tped.log"
    conda: "../envs/plink-1.9.yml"
    shell: """
        plink --bfile {params.basename} --out {params.basename} --recode transpose tab --allow-no-sex --keep-allele-order > {log} 2>&1
    """


# ----------------------------------------------------------------------------------------------- #
#                                   Random haploid calling with ANGSD                             #
# ----------------------------------------------------------------------------------------------- #

rule ANGSD_random_haploid:
    """
    Perform random pseudo-haploid variant calling with ANGSD
    """
    input:
        pileup    = rules.samtools_pileup.output.pileup,
        bamlist   = rules.generate_bam_list.output.bamlist,
        reference = config['reference'],
        fai       = config["reference"] + ".fai",
    output:
        haplos    = "results/03-variant-calling/02-ANGSD/{generation}/{generation}.haplo.gz"
    params:
        out       = "results/03-variant-calling/02-ANGSD/{generation}/{generation}"
    threads: 4
    priority: 15
    conda: "../envs/angsd-0.939.yml"
    shell: """
        angsd -pileup {input.pileup}         \
        -fai {input.fai}                     \
        -nInd $(cat {input.bamlist} | wc -l) \
        -rf <(seq 1 22)                      \
        -dohaplocall 1                       \
        -doCounts 1                          \
        -nthreads {threads}                  \
        -out {params.out}
    """


rule ANGSD_haplo_to_plink:
    """
    Convert a .haplo.gz ANGSD file to a set of PLINK .tped / .tfam files 
    """
    input:
        haplos  = rules.ANGSD_random_haploid.output.haplos,
        bamlist = rules.generate_bam_list.output.bamlist 
    output:
        tped = "results/03-variant-calling/02-ANGSD/{generation}/{generation}.tped",
        tfam = "results/03-variant-calling/02-ANGSD/{generation}/{generation}.tfam",
    params:
        outputname = "results/03-variant-calling/02-ANGSD/{generation}/{generation}"
    conda: "../envs/angsd-0.939.yml"
    priority: 15
    shell: """
        haploToPlink {input.haplos} {params.outputname}
        sed -i 's/N/0/g' {output.tped}
        cat {input.bamlist} | xargs basename -a \
                            | grep -oP '^[^.]+(?=(\.[^.]+)*(\.bam$))' \
                            | awk 'BEGIN{{OFS="\t"}}{{print "{wildcards.generation}", $1, 0, 0, 0, 0}}' \
                            > {output.tfam}
    """