from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider(retry=config['FTP']['retries']) # Anonymous 

# ---- Set config variables
configfile: "./config/config.yml"

rule download_reference_genome:
    """
    Download a reference genome from a predefined ftp URL
    human_g1k_v37 is malformed... see: https://github.com/hammerlab/biokepi/issues/117
    """
    input:
        refgen = FTP.remote("ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/{reference}.fa.gz")
    output:
        refgen = "data/refgen/GRCh37/{reference}.fa.gz"
    shell: """
        mv {input.refgen} {output.refgen}
    """


rule decompress_reference_genome:
    """
    Decrompress a reference file from .fa.gz -> .fa 
    """
    input:
        refgen = rules.download_reference_genome.output.refgen
    output:
        refgen = protected("data/refgen/GRCh37/{reference}.fa")
    log: "logs/refgen/decompress_reference_genome/{reference}.log"
    shell: """
        gunzip -v {input.refgen} > {log} 2>&1
    """

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
	    sa  = config["reference"] + ".sa",
    conda: "../envs/bwa-0.7.17.yml"
    log: "logs/refgen/index_reference_genome/index_reference_genome.log"
    shell: """
	    bwa index {input.refgen} >2 {log}
	"""

rule samtools_faidx:
    """
    Generate a `.fai` fasta index using samtools faidx
    """
    input:
        fasta = "{directory}/{fasta}.fa"
    output:
        fai   = "{directory}/{fasta}.fa.fai"
    conda: "../envs/samtools-1.15.yml"
    shell: """
        samtools faidx {input.fasta} --fai-idx {output.fai}
    """

rule split_reference_genome:
    """
    Split the reference genome according to chromosome.
    """
    input:
        refgen   = config["reference"]
    output:
        splitted = expand("data/refgen/GRCh37/splitted/{chr}.fasta", chr=range(1,23))
    shell: """
        curr_wd=`pwd`
        cd $(dirname {output.splitted[0]})
        csplit -s -z $curr_wd/{input.refgen} '/>/' '{{*}}'
        for i in xx* ; do                                  \
            n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ;         \
            mv "$i" "$n.fasta" ;                           \
        done
    """