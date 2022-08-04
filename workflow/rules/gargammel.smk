import random
import sys
configfile: "./config/config.yml"

wildcard_constraints:
    chr="\d+"

# ---- Utility functions
def get_contaminants(wildcards):
    """
    Returns a list of random ID(s) from a list of potential individuals to use as contamination.
    This will give out a single contaminating individual for each replicate, or generation.
    """
    # Open the pedigree generation file and 1000g panel definition file.
    with open(rules.run_ped_sim.input.definition) as f, open(rules.fetch_samples_panel.output.panel) as samples:
        # get the number of replicates / generations.
        gen_no = int(list(dropwhile(lambda x: x.startswith('#'), [line for line in f]))[0].split(' ')[2])

        # Fetch sample-IDs matching the user-provided contaminating population tag ('EUR', 'AFR', 'YRI', etc...)
        cont_pop_tag = config["gargammel"]["params"]["contam-pop"]
        contaminants = [sample.strip("\n").split("\t")[0] for sample in samples.readlines() if cont_pop_tag in sample]

        # Return a random list of size gen_no (one contaminating individual per replicate).
        return random.sample(contaminants, gen_no)

# ------------------------------------------------------------------------------------------------------------------- #

checkpoint get_contamination_table:
    """
    Call `get_contaminants` and output its result into a tsv file.
    """
    output:
        cont_table = "results/01-gargammel/contaminants/contaminants.tsv"
    params:
        cont      = get_contaminants
    priority: 99
    shell: """
        echo {params.cont} | awk 'BEGIN{{RS=\" \"}}{{print \"ped\"NR, $1}}' > {output.cont_table}
    """

rule create_human_contamination:
    input:
        vcf = expand(rules.concat_1000_genomes.output.merged_vcf, POP=config["gargammel"]["params"]["contam-pop"]),
        tbi = expand(f"{rules.concat_1000_genomes.output.merged_vcf}.tbi", POP=config["gargammel"]["params"]["contam-pop"]),
        chr_ref = "data/refgen/splitted/{chr}.fasta"
    output:
        hap1 = "results/01-gargammel/contaminants/{sample}/{chr}/{sample}_chr{chr}_haplo1.fasta",
        hap2 = "results/01-gargammel/contaminants/{sample}/{chr}/{sample}_chr{chr}_haplo2.fasta"
    group: "contaminate"
    resources: scatter=2
    threads: 2
    priority: 99
    conda: "../envs/bcftools-1.15.yml"
    shell: """
    bcftools consensus -H 1 -f {input.chr_ref} --sample {wildcards.sample} {input.vcf} > {output.hap1} \
    & \
    bcftools consensus -H 2 -f {input.chr_ref} --sample {wildcards.sample} {input.vcf} > {output.hap2} \
    """


def find_contaminant(wildcards):
    """
    Assign the correct contaminating individual to a given pedigree sample, according to its generation #.
    """

    # Get the current pedigree sampleID, chromosome and generation.
    sample = wildcards.sample
    chromo = wildcards.chr
    gen    = sample.split("_")[0] 
    # @TODO: what is up with these dummy files ??! I don't think we need a second context manager for that.
    # Run the checkpoint, read the contamination table and return the sampleID corresponding to the generation number.
    # return the path leading to the appropriate fasta files..
    with checkpoints.get_contamination_table.get().output.cont_table.open() as dummy, open(rules.get_contamination_table.output.cont_table) as f:
        for line in f.readlines():
            if line.startswith(f"{gen} "):
                contaminant = line.strip("\n").split(" ")[1]
                break
    #return expand("results/01-gargammel/contaminants/{sample}/{chr}", sample=contaminant, chr=chromo)
    return expand("results/01-gargammel/contaminants/{sample}/{{chr}}/{sample}_chr{{chr}}_haplo{haplo}.fasta", sample=contaminant, chr=chromo, haplo=[1,2])


rule get_consensus:
    input:
        vcf     = expand(rules.extract_twins.output.merged_vcf, POP=config["ped-sim"]["params"]["pop"]),
        chr_ref = "data/refgen/splitted/{chr}.fasta"
    output:
        hap1=temp("results/01-gargammel/{sample}/{chr}/endo/{sample}_chr{chr}_haplo1.fasta"),
        hap2=temp("results/01-gargammel/{sample}/{chr}/endo/{sample}_chr{chr}_haplo2.fasta"),
    group: "scatter"
    threads: 4
    shell:
        """
        bcftools consensus -H 1 -f {input.chr_ref} --sample {wildcards.sample} {input.vcf} > {output.hap1} \
        & \
        bcftools consensus -H 2 -f {input.chr_ref} --sample {wildcards.sample} {input.vcf} > {output.hap2} \
        """


rule fetch_bacterial_contamination:
    """
    Placeholder rule. This is the step where we'll download contaminating bacterial fastas
    """
    output:
        directory("data/gargammel/bact")
    shell: """
        mkdir -p {output}
    """


rule get_fastq:
    input:
        hap1                    = rules.get_consensus.output.hap1,
        hap2                    = rules.get_consensus.output.hap2,
        cont_table              = rules.get_contamination_table.output.cont_table,
        human_contamination     = find_contaminant,
        bacterial_contamination = rules.fetch_bacterial_contamination.output,
    output:
        forw  = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}_s1.fq.gz"),
        rev   = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}_s2.fq.gz"),
        a     = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}_a.fa.gz"),
        b     = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}.b.fa.gz"),
        c     = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}.c.fa.gz"),
        d     = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}_d.fa.gz"),
        e     = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}.e.fa.gz")
    params:
        coverage         = config['gargammel']['coverage'],
        size_freq        = config['gargammel']['sizefreq'],
        comp_endo        = config['gargammel']['comp_endo'],
        comp_cont        = config['gargammel']['comp_cont'],
        comp_bact        = config['gargammel']['comp_bact'],
        misincorporation = config['gargammel']['misincorporation'],
        output_base_name = "results/01-gargammel/{sample}/{chr}",
        input_directory  = directory("results/01-gargammel/{sample}/{chr}")
    group: "scatter"
    conda: "../envs/gargammel-1.1.2.yml"
    threads: 4
    priority: 2
    shell: """
        mkdir -p {params.output_base_name}/cont
        ln -sfrt {params.output_base_name}/cont {input.human_contamination}
        ln -sfrt {params.output_base_name} {input.bacterial_contamination}
        gargammel --comp {params.comp_bact},{params.comp_cont},{params.comp_endo}    \
                  -mapdamage {params.misincorporation} double                        \
                  -c {params.coverage}                                               \
                  -f {params.size_freq}                                              \
                  -o {params.output_base_name}/{wildcards.sample}_chr{wildcards.chr} \
                  {params.input_directory}    
    """


rule merge_chromosomes:
    input:
        forw=expand(rules.get_fastq.output.forw, chr=range(1,23), sample="{sample}"),
        rev =expand(rules.get_fastq.output.rev,  chr=range(1,23), sample="{sample}")
    output:
        forw="results/01-gargammel/fastqs/{sample}_s1.fq.gz",
        rev ="results/01-gargammel/fastqs/{sample}_s2.fq.gz"
    group: "scatter"
    priority: 3
    threads: 4
    shell: """
        zcat {input.forw} | gzip > {output.forw}
        zcat {input.rev} | gzip > {output.rev} 
    """
