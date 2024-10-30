# ---- Set config variables
configfile: "./config/config.yml"

# ------------------------------------------------------------------------------------------------------------------- #
# ---- Generic / Utility rules 

rule samtools_faidx:
    """
    Generate a `.fai` fasta index using samtools faidx
    """
    input:
        fasta = ReferenceGenome.get_path()
    output:
        fai   = ReferenceGenome.get_path() + ".fai"
    resources:
        cores = lambda w, threads: threads
    log:       f"logs/00-preprocess-reference/samtools_faidx/{ReferenceGenome.get_path()}.log"
    benchmark: f"benchmarks/00-preprocess-reference/samtools_faidx/{ReferenceGenome.get_path()}.tsv"
    conda:     "../envs/samtools-1.15.yml"
    threads:   1
    shell: """
        samtools faidx {input.fasta} --fai-idx {output.fai}
    """


# ------------------------------------------------------------------------------------------------------------------- #
# ---- 00. Download reference genome from Ensembl

#module netrules:
#    snakefile: "00-netrules.smk"
#    config: config
#
#use rule download_reference_genome from netrules


# ------------------------------------------------------------------------------------------------------------------- #
# ---- 01. Decompress the reference genome.
#          @TODO: This could be a generic rule.

rule decompress_reference_genome:
    """
    Decrompress a reference file from .fa.gz -> .fa 
    """
    input:
        refgen = ReferenceGenome.get_path() + ".gz"
    output:
        refgen = ReferenceGenome.get_path()
    resources:
        cores = lambda w, threads: threads
    log:       f"logs/00-preprocess-reference/decompress_reference_genome/{ReferenceGenome.get_name()}.log"
    benchmark: f"benchmarks/00-preprocess-reference/decompress_reference_genome/{ReferenceGenome.get_name()}.tsv"
    conda:     "../envs/coreutils-9.1.yml"
    threads:   1
    shell: """
        gunzip -v {input.refgen} > {log} 2>&1
    """


# ------------------------------------------------------------------------------------------------------------------- #
# ---- 02. Index the reference genome.

rule index_reference_genome:
    """
    Generate the Burrows-Wheeler transform on the reference genome.
    """
    input:
        refgen = ReferenceGenome.get_path()
    output:
        amb = ReferenceGenome.get_path() + ".amb",
        ann = ReferenceGenome.get_path() + ".ann",
        bwt = ReferenceGenome.get_path() + ".bwt",
        pac = ReferenceGenome.get_path() + ".pac",
        sa  = ReferenceGenome.get_path() + ".sa"
    resources:
        cores = lambda w, threads: threads
    log:       f"logs/00-preprocess-reference/index_reference_genome-{ReferenceGenome.get_name()}.log"
    benchmark: f"benchmarks/00-preprocess-reference/index_reference_genome-{ReferenceGenome.get_name()}.tsv"
    conda:     "../envs/bwa-0.7.17.yml"
    threads:   1
    shell: """
        bwa index {input.refgen} > {log} 2>&1
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 03. Split the reference genome on a per-chromosome basis

rule split_reference_genome:
    """
    Split the reference genome according to chromosome.
    """
    input:
        refgen   = ReferenceGenome.get_path()
    output:
        splitted  = expand(f"{os.path.dirname(ReferenceGenome.get_path())}/splitted/{{chr}}.fasta", chr=range(1, 23))
    resources:
        cores = lambda w, threads: threads
    log:       f"logs/00-preprocess-reference/split_reference_genome-{ReferenceGenome.get_name()}.log"
    benchmark: f"benchmarks/00-preprocess-reference/split_reference_genome-{ReferenceGenome.get_name()}.tsv"
    conda:     "../envs/coreutils-9.1.yml"
    threads:   1
    shell: """
        curr_wd=`pwd`
        cd $(dirname {output.splitted[0]}) 2> {log}
        csplit -s -z $curr_wd/{input.refgen} '/>/' '{{*}}'
        for i in xx* ; do                                  \
            n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ;         \
            mv "$i" "$n.fasta" ;                           \
        done 2>> $curr_wd/{log}
    """