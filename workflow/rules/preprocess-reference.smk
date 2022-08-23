from os.path import basename, dirname 
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider(retry=config['FTP']['retries']) # Anonymous 

# ---- Set config variables
configfile: "./config/config.yml"

reference = "Homo_sapiens.GRCh37.dna.primary_assembly"

rule download_reference_genome:
    """
    Download a reference genome from a predefined ftp URL
    human_g1k_v37 is malformed... see: https://github.com/hammerlab/biokepi/issues/117
    """
    input:
        remote = FTP.remote(expand("{url}/{reference}.gz", 
            url       = config["FTP"]["refgen"], 
            reference = lambda wildcards: basename(config["refgen"])
        ))
    output:
        reference = config["refgen"]
    params:
        url      = config["FTP"]["refgen"],
    shell: """
        mv {remove} {output.reference}
    """


rule decompress_reference_genome:
    """
    Decrompress a reference file from .fa.gz -> .fa 
    """
    input:
        reference = rules.download_reference_genome.output.reference
    output:
        reference = config["refgen"]
    shell: """
        gunzip -v {input.reference}
    """

rule index_reference_genome:
    """
    Generate the Burrows-Wheeler transform on the reference genome.
    """
    input:
        reference = rules.decompress_reference_genome.output.reference
    output:
        bwt = multiext(config["refgen"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    conda: "../envs/bwa-0.7.17.yml"
    shell: """
	    bwa index {input.reference}
	"""

rule split_refgen:
    """
    Split the reference genome according to chromosome.
    """
    input:
        reference   = rules.decompress_reference_genome.output.reference
    output:
        splitted = expand(dirname(config["refgen"]) + "/splitted/{chr}.fasta", chr=range(1,23))
    shell: """
        curr_wd=`pwd`
        cd $(dirname {output.splitted[0]})
        csplit -s -z $curr_wd/{input.reference} '/>/' '{{*}}'
        for i in xx* ; do                                  \
            n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ;         \
            mv "$i" "$n.fasta" ;                           \
        done
    """
