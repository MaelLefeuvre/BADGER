import random
import sys
from os.path import dirname
configfile: "./config/config.yml"

# ---- Utility functions
def get_contaminants(wildcards):
    """
    Returns a list of random ID(s) from a list of potential individuals to use as contamination.
    This will give out a single contaminating individual for each replicate, or generation.
    """

    # Prevent unneccessary rule re-run trigger events: If the contamination table already exists,
    # simply print out its contents and leave...
    if Path(rules.get_contamination_table.output.cont_table).exists():
        with open(rules.get_contamination_table.output.cont_table, "r") as contamination_table:
            samples = [line.strip("\n").split()[1] for line in contamination_table.readlines()]
            return samples
    # else, open the pedigree generation file and 1000g panel definition file and return a list
    # of random samples.
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
    input:
        samples_panel = rules.fetch_samples_panel.output.panel
    output:
        cont_table    = "results/01-gargammel/contaminants/contaminants.tsv",
    params:
        cont          = get_contaminants
    priority: 99
    shell: """
        echo {params.cont} | awk 'BEGIN{{RS=\" \"}}{{print \"ped\"NR, $1}}' > {output.cont_table}
    """

rule create_human_contamination:
    input:
        vcf = lambda wildcards: expand(rules.concat_1000_genomes.output.merged_vcf, 
            steps=check_ped_sim_input_vcf(wildcards), 
            POP=config["gargammel"]["params"]["contam-pop"]
        ),
        tbi = lambda wildcards: expand(f"{rules.concat_1000_genomes.output.merged_vcf}.tbi",
            steps=check_ped_sim_input_vcf(wildcards),
            POP=config["gargammel"]["params"]["contam-pop"]
        ),
        chr_ref = lambda wildcards: dirname(config["refgen"]) + "/splitted/{chr}.fasta"
    output:
        hap1 = "results/01-gargammel/contaminants/{cont}/{chr}/{cont}_chr{chr}_haplo1.fasta",
        hap2 = "results/01-gargammel/contaminants/{cont}/{chr}/{cont}_chr{chr}_haplo2.fasta"
    log: 
        hap1 = "logs/01-gargammel/create_human_contamination/{cont}_chr{chr}_haplo1.log",
        hap2 = "logs/01-gargammel/create_human_contamination/{cont}_chr{chr}_haplo2.log"
    threads: 2
    priority: 99
    conda: "../envs/bcftools-1.15.yml"
    shell: """
    bcftools consensus -e 'ALT~"<.*>"' -H 1 -f {input.chr_ref} --sample {wildcards.cont} {input.vcf} -o {output.hap1} 2> {log.hap1} \
    & \
    bcftools consensus -e 'ALT~"<.*>"' -H 2 -f {input.chr_ref} --sample {wildcards.cont} {input.vcf} -o {output.hap2} 2> {log.hap2} \
    """


def find_contaminant(wildcards):
    """
    Assign the correct contaminating individual to a given pedigree sample, according to its generation #.
    """
    # Get the current pedigree sampleID, chromosome and generation.
    sample = wildcards.sample
    chromo = wildcards.chr
    gen    = sample.split("_")[0] 
    
    # ---- Do not trigger the checkpoint if the corresponding file already exists...
    if Path(rules.get_contamination_table.output.cont_table).exists():
        contamination_table = rules.get_contamination_table.output.cont_table
    else:
        contamination_table = checkpoints.get_contamination_table.get().output.cont_table

    with open(contamination_table) as f:
        for line in f.readlines():
            if line.startswith(f"{gen} "):
                contaminant = line.strip("\n").split(" ")[1]
                break
    
    out = "results/01-gargammel/contaminants/{cont}/{{chr}}/{cont}_chr{{chr}}_haplo{haplo}.fasta"
    return expand(out, cont=contaminant, chr=chromo, haplo=[1,2])


rule get_consensus:
    """
    Generate phased haplotypes for pedigree individuals using ped-sim's output VCF and a provided reference fasta file.
    (Copy Number Variants, i.e. 'CN[0-9]+' are unsupported by bcftools, thus they must be excluded.)
    """
    input:
        vcf     = expand(rules.extract_twins.output.merged_vcf, POP=config["ped-sim"]["params"]["pop"]),
        chr_ref = lambda wildcards: dirname(config["refgen"]) + "/splitted/{chr}.fasta"
    output:
        hap1=temp("results/01-gargammel/{sample}/{chr}/endo/{sample}_chr{chr}_haplo1.fasta"),
        hap2=temp("results/01-gargammel/{sample}/{chr}/endo/{sample}_chr{chr}_haplo2.fasta"),
    log: 
        hap1 = "logs/01-gargammel/get_consensus/{sample}_chr{chr}_haplo1.log",
        hap2 = "logs/01-gargammel/get_consensus/{sample}_chr{chr}_haplo2.log"
    group: "scatter"
    threads: 2
    shell: """
        bcftools consensus -e 'ALT~"<CN[0-9]+>"' -H 1 -f {input.chr_ref} --sample {wildcards.sample} {input.vcf} -o {output.hap1} 2> {log.hap1} \
        & \
        bcftools consensus -e 'ALT~"<CN[0-9]+>"' -H 2 -f {input.chr_ref} --sample {wildcards.sample} {input.vcf} -o {output.hap2} 2> {log.hap2} \
    """


rule fetch_bacterial_contamination:
    """
    Placeholder rule. This is the step where we'll download contaminating bacterial fastas
    As of now, stolen from gargammel's Makefile 'makebacterialex' target
    """
    output:
        out_dir        = directory("data/gargammel/bact"),
        fasta_dir      = directory("data/gargammel/bact/bact"),
        abundance_list = "data/gargammel/bact/bact/list",
        taxa_list      = "data/gargammel/bact/all_taxa.tsv"
    params:
        url = "https://www.dropbox.com/s/obmr48d72ahjvhp/clovis.tar.gz?dl=1"
        #url = "https://www.dropbox.com/s/1pdbqbguw0jfzib/k14.tar.gz?dl=1"
    shell: """
        cd {output.out_dir}

        # Download bacterial data
        wget -O- {params.url} | tar xvzf - --strip-components=1

        # Rename "fasta" directory to "bact"
        mv fasta/* bact/ 
    """


rule run_gargammel:
    input:
        hap1                    = rules.get_consensus.output.hap1,
        hap2                    = rules.get_consensus.output.hap2,
        cont_table              = rules.get_contamination_table.output.cont_table,
        human_contamination     = find_contaminant,
        bacterial_contamination = rules.fetch_bacterial_contamination.output.fasta_dir,
    output:
        forwd = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}_s1.fq.gz"),
        revrs = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}_s2.fq.gz"),
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
        input_directory  = directory("results/01-gargammel/{sample}/{chr}"),
    log: "logs/01-gargammel/run_gargammel/{sample}_chr{chr}.log"
    group: "scatter"
    conda: "../envs/gargammel-1.1.2.yml"
    threads: 1
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
                  {params.input_directory} > {log} 2>&1
    """


rule merge_chromosomes:
    input:
        forwd = expand(rules.run_gargammel.output.forwd, chr=range(1,23), sample="{sample}"),
        revrs = expand(rules.run_gargammel.output.revrs, chr=range(1,23), sample="{sample}")
    output:
        forwd = "results/02-preprocess/00-raw/{sample}_s1.fq.gz",
        revrs = "results/02-preprocess/00-raw/{sample}_s2.fq.gz"
    log: 
        forwd = "logs/01-gargammel/merge_chromosomes/{sample}_s1.log",
        revrs = "logs/01-gargammel/merge_chromosomes/{sample}_s2.log"
    group: "scatter"
    priority: 3
    threads: 2
    shell: """
        zcat {input.forwd} | gzip > {output.forwd} 2> {log.forwd}
        zcat {input.revrs} | gzip > {output.revrs} 2> {log.revrs}
    """
