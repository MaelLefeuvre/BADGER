# ---- Set config variables
configfile: "./config/config.yml"

localrules: symlink_original_data


# ------------------------------------------------------------------------------------------------------------------- #
# ---- Generalist / Utility functions. 

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


# ------------------------------------------------------------------------------------------------------------------- #
# ---- 00. Download the original dataset.

#module netrules:
#    snakefile: "00-netrules.smk"
#    config: config
#
#use rule download_1000_genomes from netrules
#use rule fetch_samples_panel from netrules

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 00. B symlink out of the original data dir if required.

rule symlink_original_data:
    """
    Symlink out of the original data dir if required. Protects the raw input from any tampering + useful if there are
    no requested filtering steps.
    """
    input:
        #vcf = rules.download_1000_genomes.output.vcf
        vcf = "data/vcf/1000g-phase3/00-original/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    output:
        link = temp("data/vcf/1000g-phase3/01-filter/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
    threads: 1
    log: "logs/00-preprocess-1000g/symlink_original_data/symlink_original_data-chr-{chr}.log"
    shell: """
        ln -srt data/vcf/1000g-phase3/01-filter/ {input.vcf} > {log} 2>&1
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 01. Filter out INDELS and CNVs.

rule filter_indels_1000_genomes:
    """
    Filter each 1000g variant callset to merely preserve biallelic SNPs
    """
    input:
        link = rules.symlink_original_data.output.link
    output:
        bcf = pipe("data/vcf/1000g-phase3/01-filter/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.snps.bcf")
    log: "logs/00-preprocess-1000g/filter_indels_1000_genomes/filter_indels_1000_genomes-chr{chr}.log"
    conda: "../envs/bcftools-1.15.yml"
    threads: workflow.cores / 22
    shell: """
        bcftools view --threads {threads} -Ob -v snps < {input.link} > {output.bcf} 2> {log}
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 02. Filter out multiallelic variants.

def choose_filter_indel(wildcards):
    if config["ped-sim"]["filter-indels"]:
        return rules.filter_indels_1000_genomes.output.bcf
    else:
        return rules.symlink_original_data.output.link


rule normalize_1000_genomes:
    """
    Normalizes a 1000g vcf file: This step is required to ensure there are no 
    false negatives when filtering for bi-allelic variants, since the 1000G 
    project may contain duplicate records for a given position.
    """
    input:
        bcf = choose_filter_indel
    output:
        bcf = temp("data/vcf/1000g-phase3/01-filter/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.{steps}.norm.bcf")
    log: "logs/00-preprocess-1000g/normalize_1000_genomes/normalize_1000_genomes.{steps}.chr{chr}.log"
    conda: "../envs/bcftools-1.15.yml"
    threads: workflow.cores / 22
    shell: """
        bcftools norm --threads {threads} --multiallelics +both -Ob {input.bcf} > {output.bcf} 2> {log}
    """

rule filter_biallelic_1000_genomes:
    """
    Filter-out : - all multi-allelic positions within a 1000g vcf file.
                 - unphased positions. Yes. some positions are unphased (e.g.: 12:6608369) and can break ped-sim.
    @TODO: Getting an error if I create a named pipe output. Snakemake considers this rule as having two consumers, 
           despite it not being the case (?) : the only rule with several consumers seem to be `concat_1000_genomes`
    """
    input:
        bcf = rules.normalize_1000_genomes.output.bcf
    output:
        bcf = temp("data/vcf/1000g-phase3/01-filter/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.{steps}.norm.m2.M2.bcf")
    log: "logs/00-preprocess-1000g/filter_biallelic_1000_genomes/filter_biallelic_1000_genomes.{steps}.chr{chr}.log"
    conda: "../envs/bcftools-1.15.yml"
    threads: workflow.cores / 22
    shell: """
        bcftools view --threads {threads} --phased --min-alleles 2 --max-alleles 2 -Ob < {input.bcf} > {output.bcf} 2> {log}
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 03. Extract individuals from a given (super)-population and subset our dataset.


rule get_target_pop_samples:
    """
    Create a subset list of samples ID for a specific population or superpopulation.
    """
    input:
        # panel       = rules.fetch_samples_panel.output.panel
        panel = "data/vcf/1000g-phase3/samples-list/integrated_call_samples_v3.20130502.ALL.panel"
    output:
        target_list = "data/vcf/1000g-phase3/samples-list/integrated_call_samples_v3.20130502.{POP}.panel"
    log: "logs/00-preprocess-1000g/get_target_pop_samples/get_target_pop_samples-{POP}.log"
    shell: """
        grep {wildcards.POP} {input.panel} | cut -f1 > {output.target_list} 2> {log}
    """


rule subset_1000_genomes:
    """
    Subset the 1000g-phase3 dataset, keeping only the specified population.
    """
    input:
        bcf     = rules.filter_biallelic_1000_genomes.output.bcf,
	    samples = rules.get_target_pop_samples.output.target_list
    output:
        bcf     = "data/vcf/1000g-phase3/02-subset/{POP}/{POP}.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.{steps}.norm.m2.M2.bcf"
    log: "logs/00-preprocess-1000g/subset_1000_genomes/{POP}/subset_1000_genomes.{steps}.chr{chr}.log"
    conda: "../envs/bcftools-1.15.yml"
    group: "1000g-preprocess"
    threads: 1
    shell: """
        bcftools view --threads {threads} --samples-file {input.samples} -Ob < {input.bcf} > {output.bcf} 2> {log}
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 04. Perform minor allele frequency filtration

rule update_allele_frequencies:
    """
    Use bcftools +fill-tags to update the "AF" tag contained within the 'INFO' column.
    This ensures proper minor allele frequency removal. 
    Using a -e predicate to match "{EUR}_AF" INFO tags would be a simpler alternative,
    but would prohibit us from targeting sub-populations
    """
    input:
        bcf = rules.subset_1000_genomes.output.bcf
    output:
        bcf = pipe("data/vcf/1000g-phase3/02-subset/{POP}/{POP}.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.{steps}.norm.m2.M2.retagged.bcf")
    log: "logs/00-preprocess-1000g/update_allele_frequencies/{POP}/update_allele_frequencies.{steps}.chr{chr}.log"
    conda: "../envs/bcftools-1.15.yml"
    group: "1000g-preprocess"
    threads: 1
    shell: """
        bcftools +fill-tags -Ob --threads {threads} -- -t AF < {input.bcf} > {output.bcf} 2> {log}
    """

rule filter_maf_1000_genomes:
    """
    filter out alleles that are below the requested minor frequency.
    """
    input: 
        bcf = rules.update_allele_frequencies.output.bcf
    output:
        bcf = pipe("data/vcf/1000g-phase3/02-subset/{POP}/{POP}.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.{steps}.norm.m2.M2.maf05.bcf")
    params:
        maf = config['ped-sim']['filter-maf-threshold']
    log: "logs/00-preprocess-1000g/filter_maf_1000_genomes/{POP}/filter_maf_1000_genomes.{steps}.chr{chr}.log"
    conda: "../envs/bcftools-1.15.yml"
    group: "1000g-preprocess"
    threads: 1
    shell: """
        bcftools view --threads {threads} --min-af {params.maf}:minor -Ob -o {output.bcf} < {input.bcf} 2> {log}
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 05. Concatenate our dataset across chromosomes.


def check_maf_option(wildcards):
    """
    Decine on the appropriate input for concat_1000_genome, based whether or
    not the user asked for maf filtration.
    """
    if config["ped-sim"]["filter-maf"]:
        output = rules.filter_maf_1000_genomes.output.bcf
    else:
        output = rules.subset_1000_genomes.output.bcf

    ## Redirect to the appropriate path if filter-indels was requested
    steps ="genotypes.snps" if config["ped-sim"]["filter-indels"] else "genotypes"
    return expand(output, chr=range(1,23), POP="{POP}", steps=steps)


rule concat_1000_genomes:
    """
    Concatenate chromosome file into a single, population-specific VCF. 
    """
    input:
        split_vcfs = check_maf_option
    output:
        merged_vcf = expand(
            "data/vcf/1000g-phase3/03-merged/{POP}/{POP}.merged.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes{indel_tag}.norm.m2.M2{maf_tag}.vcf.gz",
            POP="{POP}", 
            indel_tag = ".snps" if config["ped-sim"]["filter-indels"] else "",
            maf_tag   = ".maf" + str(config["ped-sim"]["filter-maf-threshold"]).split(".")[1] if config["ped-sim"]["filter-maf"] else ""
        )
    threads: 1
    conda: "../envs/bcftools-1.15.yml"
    log: "logs/00-preprocess-1000g/concat_1000_genomes/concat_1000_genomes-{POP}.log"
    group: "1000g-preprocess"
    shell: """
        bcftools concat --threads {threads} -Oz -o {output.merged_vcf} {input.split_vcfs} > {log} 2>&1
    """