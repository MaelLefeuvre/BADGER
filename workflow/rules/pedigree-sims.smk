from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from itertools import dropwhile
import os 

HTTP = HTTPRemoteProvider()

# ---- Set config variables

configfile: "./config/config.yml"

# ---- Utility functions
def get_twins(wildcards):
    """
    Runs through the pedigree definition file, finds which individuals are supposed to become twins and outputs their IDs
    """
    with open(rules.run_ped_sim.input.definition) as f:
        # Find how many pedigree replicates are being run
        gen_no = int(list(dropwhile(lambda x: x.startswith('#'), [line for line in f]))[0].split(' ')[2])
    with open(config["ped-sim"]["data"]["codes"], "r") as f:
        # Run through the pedigree_codes and find potential twins relationship
        lines = [line for line in f.read().split('\n') if line.strip() != ''] 
        twins = []
        for line in lines[:len(lines)]:
            cells = line.split("\t")
            twins.append(cells[1]) if cells[1] == cells[2].lower() else None
    if twins:
        return expand("ped{gen}_{id}", gen=range(1,gen_no+1), id=twins)
    else:
        raise RuntimeError(f"Failed to extract any twins within pedigree codes file: '{f.name}'")

# ------------------------------------------------------------------------------------------------------------------- #

rule fetch_sex_specific_gen_map:
    """
    See: Bhérer, C., Campbell, C. & Auton, A. Refined genetic maps reveal sexual dimorphism in human meiotic
         recombination at multiple scales. Nat Commun 8, 14994 (2017). https://doi.org/10.1038/ncomms14994

     - https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination
    """
    input:
        gen_map = HTTP.remote(config["ped-sim"]["input"]["refined-genetic-map-url"])
    output:
        map_dir  = directory("data/ped-sim/Refined_genetic_map_b37"),
        gen_maps = temp(expand("data/ped-sim/Refined_genetic_map_b37/{sex}_chr{chrom}.txt", chrom=range(1,23), sex=["female", "male", "sexavg"]))
    shell: """
        tar --strip-components=1 -xvzf {input.gen_map} -C {output.map_dir} && rm -rf {input.gen_map}
    """


rule format_sex_specific_gen_map:
    """
    Download a sex-specific genetic map. Required by ped-sim
    """
    input:
        gen_map = rules.fetch_sex_specific_gen_map.output.gen_maps
    output:
        sim_map = config["ped-sim"]['data']["map"]
    params:
        map_dir = rules.fetch_sex_specific_gen_map.output.map_dir
    shell: """
        printf "#chr\tpos\tmale_cM\tfemale_cM\n" > {output.sim_map};                    \
        for chr in {{1..22}}; do                                                        \
            paste {params.map_dir}/male_chr$chr.txt {params.map_dir}/female_chr$chr.txt \
                | awk -v OFS="\t" 'NR > 1 && $2 == $6 {{print $1,$2,$4,$8}}'            \
                | sed 's/^chr//' >> {output.sim_map};                                   \
        done
    """

rule fetch_interference_map:
    """
    Download a sex-refined genetic interference map. Required for ped-sim
    See: Campbell, C., Furlotte, N., Eriksson, N. et al. Escape from crossover interference increases with maternal age.
         Nat Commun 6, 6260 (2015). https://doi.org/10.1038/ncomms7260
     - https://github.com/williamslab/ped-sim/blob/master/interfere/nu_p_campbell.tsv
    """
    input:
        intf_map = HTTP.remote(config["ped-sim"]["input"]["interference-map-url"])
    output:
        intf_map = config['ped-sim']['data']['interference']
    shell: """
        mv {input.intf_map} {output.intf_map}
    """

def check_ped_sim_input_vcf(wildcards):
    filter_indels = config["ped-sim"]["filter-indels"]
    if filter_indels:
        return "genotypes.snps" 
    else:
        return "genotypes"


rule run_ped_sim:
    """
    Run pedigree-simulator, using individuals from a specified population of the 1000g project as founder individuals.
    """
    input:
        vcf          = lambda wildcards: expand(rules.concat_1000_genomes.output.merged_vcf, steps=check_ped_sim_input_vcf(wildcards), POP="{POP}") ,
        definition   = config['ped-sim']['data']['definition'],
        map          = config['ped-sim']['data']['map'],
        interference = config['ped-sim']['data']['interference']
    output:
        fam = "results/00-ped-sim/{POP}-pedigrees-everyone.fam",
        log = "results/00-ped-sim/{POP}-pedigrees.log",
        seg = "results/00-ped-sim/{POP}-pedigrees.seg",
        vcf = "results/00-ped-sim/{POP}-pedigrees.vcf.gz",
    params:
        output_basename  = "results/00-ped-sim/{POP}-pedigrees",
        error_rate       = config['ped-sim']['params']['error_rate'],
        missingness      = config['ped-sim']['params']['missingness']
    log: "logs/00-ped-sim/run_ped_sim/{POP}_ped-sim.log"
    conda: "../envs/ped-sim-1.3.yml"
    shell: """
        ped-sim -d {input.definition}            \
                -m {input.map}                   \
                -i {input.vcf}                   \
                -o {params.output_basename}      \
                --intf {input.interference}      \
                --fam                            \
                --keep_phase                     \
                --miss_rate {params.missingness} \
                --err_rate {params.error_rate} > {log} 2>&1 
    """

rule filter_snps:
    """
    Keep only biallelic SNPs from the definition file and convert to BGZIP
    """
    input:
        rules.run_ped_sim.output.vcf
    output:
        vcf = "results/00-ped-sim/{POP}-pedigrees-M2-m2-snps.vcf.gz"
    conda: "../envs/bcftools-1.15.yml"
    threads: 16
    shell: """
        bcftools view --threads {threads} -M2 -m2 -v snps {input} | bgzip -c > {output} && tabix {output}
    """

rule extract_twins:
    """
    Pick a selected individual within ped-sim's vcf output (we run through our pedigree_codes file using get_twins()
    and search for any self-comparison).
    Copy the genotype information as a separate sample and create a merged vcf.
    """
    input:
        vcf        = rules.filter_snps.output.vcf
    output:
        twins_vcf  = temp("results/00-ped-sim/{POP}-pedigrees-M2-m2-snps-twins.vcf.gz"),
        twin_codes = "results/00-ped-sim/{POP}-twin_codes.txt",
        merged_vcf = "results/00-ped-sim/{POP}-pedigrees-M2-m2-snps-merged.vcf.gz"
    params:
        twins = get_twins
    conda: "../envs/bcftools-1.15.yml"
    shell: """
        echo -e {params.twins} | tr '[:lower:]' '[:upper:]' | sed 's/PED/ped/g' | awk 'BEGIN{{RS=" "}}{{print}}' | head -n -1 > {output.twin_codes}
        bcftools view -s $(echo {params.twins} | sed 's/ /,/g') {input.vcf} | bcftools reheader -s {output.twin_codes} | bgzip > {output.twins_vcf}
        tabix {output.twins_vcf}
        bcftools merge -Oz {input.vcf} {output.twins_vcf} > {output.merged_vcf}
        tabix {output.merged_vcf}
    """


checkpoint get_samples:
    """
    Extract the ped-sim samples identifiers contained within the VCF. These will become our wildcards
    for the remainder of this simulations pipeline.
    """
    input:
        expand(rules.extract_twins.output.merged_vcf, POP=config["ped-sim"]["params"]["pop"])
    output:
        "results/00-ped-sim/sample_names.tsv"
    priority: 50
    shell: """
        zcat {input} | grep '#CHROM' | cut -f 10- > {output}
    """
