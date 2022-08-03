from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider(retry=config['FTP']['retries']) # Anonymous 

# ---- Set config variables
configfile: "./config/config.yml"

ref_url  = config["refgen"]["url"]
reference="Homo_sapiens.GRCh37.dna.primary_assembly"

rule download_refgen:
    """
    Download a reference genome from a predefined ftp URL
    human_g1k_v37 is malformed... see: https://github.com/hammerlab/biokepi/issues/117
    """
    input:
        refgen = FTP.remote(expand("ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/{reference}.fa.gz", reference=reference))
    output:
        refgen = "data/refgen/{reference}.fa.gz"
    shell:
        "mv {input.refgen} {output.refgen}"


rule decompress_refgen:
    """
    Decrompress a reference file from .fa.gz -> .fa 
    """
    input:
        refgen = rules.download_refgen.output.refgen
    output:
        refgen = "data/refgen/{reference}.fa"
    shell:
        "gunzip -v {input.refgen}"

rule index_refgen:
    """
    Generate the Burrows-Wheeler transform on the reference genome.
    """
    input:
        refgen = rules.decompress_refgen.output.refgen
    output:
        amb = "data/refgen/{reference}.fa.amb",
	    ann = "data/refgen/{reference}.fa.ann",
	    bwt = "data/refgen/{reference}.fa.bwt",
	    pac = "data/refgen/{reference}.fa.pac",
	    sa  = "data/refgen/{reference}.fa.sa",
    conda: "../envs/bwa-0.7.17.yml"
    shell:
        """
	bwa index {input.refgen}
	"""

rule split_refgen:
    """
    Split the reference genome according to chromosome.
    """
    input:
        refgen   = expand(rules.decompress_refgen.output.refgen, reference=reference)
    output:
        splitted = expand("data/refgen/splitted/{chr}.fasta", chr=range(1,23))
    shell:
        """
        curr_wd=`pwd`
        cd $(dirname {output.splitted[0]})
        csplit -s -z $curr_wd/{input.refgen} '/>/' '{{*}}'
        for i in xx* ; do \
            n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ; \
            mv "$i" "$n.fasta" ; \
        done
        """
