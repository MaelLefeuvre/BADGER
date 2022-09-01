from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider(retry=config['FTP']['retries']) # Anonymous 

# ---- Set config variables
configfile: "./config/config.yml"
# ------------------------------------------------------------------------------------------------------------------- #

rule download_1000_genomes:
    """
    Download 1000genomes phase 3 SNPs
    http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

    """
    input:
        vcf = FTP.remote(expand("{url}/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
            url=config["FTP"]["1000g"],
            chr = "{chr}"
        ))
    output:
        vcf = "data/g1k-phase3-callset/00-original/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    log: "logs/1000g/download_1000_genomes/download_1000_genomes-chr{chr}.log"
    shell: """
        mv {input.vcf} {output.vcf} 2> {log}
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

rule filter_indels_1000_genomes:
    """
    Filter each 1000g variant callset to merely preserve biallelic SNPs
    """
    input:
        vcf = rules.download_1000_genomes.output.vcf
    output:
        vcf = "data/g1k-phase3-callset/00-original/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.snps.vcf.gz"
    conda: "../envs/bcftools-1.15.yml"
    threads: 16
    shell: """
        bcftools view -Oz -v snps < {input.vcf} > {output.vcf}
    """


rule normalize_1000_genomes:
    """
    Normalizes a 1000g vcf file: This step is required to ensure there are no 
    false negatives when filtering for bi-allelic variants, since the 1000G 
    project may contain duplicate records for a given position.
    """
    input:
        vcf = "data/g1k-phase3-callset/00-original/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.{steps}.vcf.gz"
    output:
        vcf = temp("data/g1k-phase3-callset/01-biallelic/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.{steps}.norm.vcf.gz")
    conda: "../envs/bcftools-1.15.yml"
    threads: 16
    shell: """
        bcftools norm --threads {threads} --multiallelics +both -Oz {input.vcf} > {output.vcf}
    """

rule filter_biallelic_1000_genomes:
    """
    Filter-out all multi-allelic positions within a 1000g vcf file.
    @TODO: Getting an error if I create a named pipe output. Snakemake considers this rule as having two consumers, 
           despite it not being the case (?) : the only rule with several consumers seem to be `concat_1000_genomes`
    """
    input:
        vcf = rules.normalize_1000_genomes.output.vcf
    output:
        vcf = "data/g1k-phase3-callset/01-biallelic/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.{steps}.norm.m2.M2.vcf.gz"
    conda: "../envs/bcftools-1.15.yml"
    threads: 16
    shell: """
        bcftools view --threads {threads} --min-alleles 2 --max-alleles 2 -Oz < {input.vcf} > {output.vcf}
    """

rule subset_1000_genomes:
    """
    Subset the 1000g-phase3 dataset, keeping only the specified population.
    """
    input:
        vcf     = rules.filter_biallelic_1000_genomes.output.vcf,
	    samples = rules.get_target_pop_samples.output.target_list
    output:
        vcf     = "data/g1k-phase3-callset/02-subset/{POP}/{POP}.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.{steps}.norm.m2.M2.vcf.gz"
    conda: "../envs/bcftools-1.15.yml"
    threads: 16
    shell: """
        bcftools view --threads {threads} --samples-file {input.samples} -Oz < {input.vcf} > {output.vcf} 
    """

rule update_allele_frequencies:
    """
    Use bcftools +fill-tags to update the "AF" tag contained within the 'INFO' column.
    This ensures proper minor allele frequency removal.
    """
    input:
        vcf = rules.subset_1000_genomes.output.vcf
    output:
        vcf = temp("data/g1k-phase3-callset/02-subset/{POP}/{POP}.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.{steps}.norm.m2.M2.retagged.vcf.gz")
    conda: "../envs/bcftools-1.15.yml"
    threads: 16
    shell: """
        bcftools +fill-tags -Oz --threads {threads} -- -t AF < {input.vcf} > {output.vcf}
    """

rule filter_maf_1000_genomes:
    """
    filter out fixed alleles
    """
    input: 
        vcf = rules.update_allele_frequencies.output.vcf
    output:
        vcf = "data/g1k-phase3-callset/02-subset/{POP}/{POP}.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.{steps}.norm.m2.M2.maf05.vcf.gz"
    params:
        maf = 0.05
    threads: 16
    conda: "../envs/bcftools-1.15.yml"
    shell: """
        bcftools view --threads {threads} --min-af {params.maf}:minor -Ou -o {output.vcf} < {input.vcf} 
    """

def check_maf_option(wildcards):
    if config["ped-sim"]["filter-maf"]:
        output = rules.filter_maf_1000_genomes.output.vcf
    else:
        output = rules.subset_1000_genomes.output.vcf
    return expand(output, chr=range(1,23), POP="{POP}", steps="{steps}")


rule concat_1000_genomes:
    """
    Concatenate chromosome file into a single, population-specific VCF. 
    """
    input:
        split_vcfs = check_maf_option
    output:
        merged_vcf = "data/g1k-phase3-callset/03-merged/{POP}/{POP}.merged.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.{steps}.vcf.gz"
    threads: 16
    conda: "../envs/bcftools-1.15.yml"
    log: "logs/1000g/concat_1000_genomes/concat_1000_genomes-{POP}-{steps}.log"
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