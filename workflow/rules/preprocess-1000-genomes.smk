from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider(retry=config['FTP']['retries']) # Anonymous 

# ---- Set config variables
configfile: "./config/config.yml"
# ------------------------------------------------------------------------------------------------------------------- #

rule download_1000_genomes:
    """
    Download 1000genomes phase 3 SNPs
    http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

    """
    input:
        vcf = FTP.remote(expand("{url}/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
            url=config["FTP"]["1000g"],
            chrom = "{chrom}"
        ))
    output:
        vcf = "data/g1k-phase3-callset/00-original/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    log: "logs/1000g/download_1000_genomes/download_1000_genomes-chr{chrom}.log"
    shell: """
        mv {input.vcf} {output.vcf} 2> {log}
    """

rule filter_1000_genomes:
    """
    Filter each 1000g variant callset to merely preserve biallelic SNPs with maf >= 0.05
    """
    input:
        vcf = rules.download_1000_genomes.output.vcf
    output:
        vcf = "data/g1k-phase3-callset/01-filtered/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.m2.M2.maf05.vcf.gz"
    threads: 16
    conda: "../envs/bcftools-1.15.yml"
    shell: """
        bcftools view -Oz --threads {threads} --min-alleles 2 --max-alleles 2 -v snps --min-af 0.05 {input.vcf} -o {output.vcf}
    """


rule fetch_samples_panel:
    """
    Download samples metadata from the 1000g FTP website
    """
    input:
        panel = FTP.remote("ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel")
    output:
        panel = "data/g1k-phase3-callset/samples-list/integrated_call_samples_v3.20130502.ALL.panel"
    shell: """
        mv {input.panel} {output.panel}
    """


rule get_target_pop_samples:
    """
    Create a subset list of samples ID for a specific population or superpopulation.
    """
    input:
        panel       = rules.fetch_samples_panel.output.panel
    output:
        target_list = "data/g1k-phase3-callset/samples-list/integrated_call_samples_v3.20130502.{POP}.panel"
    shell: """
        grep {wildcards.POP} {input.panel} | cut -f1 > {output.target_list}
    """


rule subset_1000_genomes:
    """
    Subset the 1000g-phase3 dataset, keeping only the specified population.
    """
    input:
        vcf     = rules.filter_1000_genomes.output.vcf,
	    samples = rules.get_target_pop_samples.output.target_list
    output:
        vcf     = "data/g1k-phase3-callset/01-filtered/{POP}.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.m2.M2.maf05.vcf.gz"
    log: "logs/1000g/subset_1000_genomes/subset_1000_genomes-{POP}-chr{chrom}.log"
    threads: 16
    conda: "../envs/bcftools-1.15.yml"
    shell: """
        bcftools view --threads {threads}            \
                      --samples-file {input.samples} \
                      -Oz {input.vcf}                \
                      -o {output.vcf} > {log} 2>&1
    """


rule concat_1000_genomes:
    """
    Concatenate chromosome file into a single, population-specific VCF. 
    """
    input:
        split_vcfs = expand(rules.subset_1000_genomes.output.vcf, chrom=range(1,23), POP="{POP}")
    output:
        merged_vcf = "data/g1k-phase3-callset/02-merged/{POP}.merged.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.m2.M2.maf05.vcf.gz"
    threads: 16
    conda: "../envs/bcftools-1.15.yml"
    log: "logs/1000g/concat_1000_genomes/concat_1000_genomes-{POP}.log"
    shell: """
        bcftools concat --threads {threads} -Oz -o {output.merged_vcf} {input.split_vcfs} > {log} 2>&1
    """


rule tabix_vcf:
    """
    Index a generic .vcf[.gz] file.
    """
    input:
        vcf = "{vcf}"
    output:
        tbi = "{vcf}.tbi"
    conda: "../envs/bcftools-1.15.yml"
    shell: """
        tabix {input.vcf}
    """