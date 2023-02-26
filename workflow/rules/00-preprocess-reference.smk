# ---- Set config variables
configfile: "./config/config.yml"

# ------------------------------------------------------------------------------------------------------------------- #
# ---- Generic / Utility rules 

rule samtools_faidx:
    """
    Generate a `.fai` fasta index using samtools faidx
    """
    input:
        fasta = "{directory}/{fasta}.fa"
    output:
        fai   = "{directory}/{fasta}.fa.fai"
    resources:
        cores = lambda w, threads: threads
    log:       "logs/generics/{directory}/samtools_faidx-{fasta}.log"
    benchmark: "benchmarks/generics/{directory}/samtools_faidx-{fasta}.tsv"
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
        # refgen = rules.download_reference_genome.output.refgen
        refgen = "data/refgen/GRCh37/{reference}.fa.gz"
    output:
        refgen = protected("data/refgen/GRCh37/{reference}.fa")
    resources:
        cores = lambda w, threads: threads
    log:       "logs/00-preprocess-reference/decompress_reference_genome/{reference}.log"
    benchmark: "benchmarks/00-preprocess-reference/decompress_reference_genome/{reference}.tsv"
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
        refgen = config["reference"]
    output:
        amb = config["reference"] + ".amb",
        ann = config["reference"] + ".ann",
        bwt = config["reference"] + ".bwt",
        pac = config["reference"] + ".pac",
        sa  = config["reference"] + ".sa"
    resources:
        cores = lambda w, threads: threads
    log:       "logs/00-preprocess-reference/index_reference_genome.log"
    benchmark: "benchmarks/00-preprocess-reference/index_reference_genome.tsv"
    conda:     "../envs/bwa-0.7.17.yml"
    threads:   1
    shell: """
        bwa index {input.refgen} >2 {log}
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 03. Split the reference genome on a per-chromosome basis

rule split_reference_genome:
    """
    Split the reference genome according to chromosome.
    """
    input:
        refgen   = config["reference"]
    output:
        splitted = expand("data/refgen/GRCh37/splitted/{chr}.fasta", chr=range(1,23))
    resources:
        cores = lambda w, threads: threads
    log:       "logs/00-preprocess-reference/split_reference_genome.log"
    benchmark: "benchmarks/00-preprocess-reference/split_reference_genome.tsv"
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