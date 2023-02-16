from itertools import dropwhile
import os 


#module netrules:
#    snakefile: "00-netrules.smk"
#    config: config
#
#use rule fetch_sex_specific_recombination_map from netrules
#use rule fetch_interference_map from netrules

# ---- Set config variables
configfile: "./config/config.yml"

localrules: get_samples, format_sex_specific_gen_map

# ---- Utility functions
def get_twins(wildcards):
    """
    Runs through the pedigree definition file, finds which individuals are supposed to become twins and outputs their IDs
    """
    #with open(rules.run_ped_sim.input.definition) as f:
    #    # Find how many pedigree replicates are being run
    #    gen_no = int(list(dropwhile(lambda x: x.startswith('#'), [line for line in f]))[0].split(' ')[2])

    gen_no = config['ped-sim']['replicates']

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

rule format_sex_specific_gen_map:
    """
    Concatenate a set of sex-specific chromosome genetic map into a single file.
    """
    input:
        #gen_map = rules.fetch_sex_specific_recombination_map.output.gen_maps
        gen_map = expand("data/recombination-maps/Refined_genetic_map_b37/{sex}_chr{chrom}.txt", chrom=range(1,23), sex=["female", "male", "sexavg"])
    output:
        sim_map = config["ped-sim"]['data']["map"]
    params:
        #map_dir = rules.fetch_sex_specific_recombination_map.output.map_dir
        map_dir = directory("data/recombination-maps/Refined_genetic_map_b37"),
    log:   "logs/00-ped-sim/format_sex_specific_gen_map.log"
    shell: """
        printf "#chr\tpos\tmale_cM\tfemale_cM\n" > {output.sim_map}; \
        for chr in {{1..22}}; do \
            paste {params.map_dir}/male_chr$chr.txt {params.map_dir}/female_chr$chr.txt \
                | awk -v OFS="\t" 'NR > 1 && $2 == $6 {{print $1,$2,$4,$8}}' \
                | sed 's/^chr//' >> {output.sim_map}; \
        done 2> {log}
    """


def check_ped_sim_input_vcf(wildcards):
    """
    This is the part were we resolve the input to ped-sim according to whether or not the user
    requested INDEL filtration.
    """
    filter_indels = config["ped-sim"]["filter-indels"]
    if filter_indels:
        return "genotypes.snps" 
    else:
        return "genotypes"

def set_ped_sim_seed(wildcards):
    seed              = config['ped-sim']['params']['seed']
    chained_error_msg = RuntimeError(f"Invalid seed value for ped-sim: '{seed}'. Check your config file.")
    
    if seed is None:
        with open(rules.meta.output.metadata) as f:
            metadata = yaml.load(f, Loader=yaml.loader.SafeLoader)
            seed     = metadata['seed']
        
    return seed


rule format_pedigree_definition:
    """
    @ TODO use jinja / YTE to template.
    # See: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#template-rendering-integration
    """
    input:
        template = config['ped-sim']['data']['definition']
    output:
        definition = "results/00-ped-sim/pedigree.def"
    params:
        replicates = config['ped-sim']['replicates']
    log: "logs/00-ped-sim/format_pedigree_definition.log"
    shell: """
        REPLICATES={params.replicates} envsubst < {input.template} > {output.definition} 2> {log}
    """


rule run_ped_sim:
    """
    Run pedigree-simulator, using individuals from a specified population of the 1000g project as founder individuals.
    @ TODO: Maybe the fact that we're fetching the interference map is the cause of the FTP os error. ?
    """
    input:
        vcf          = rules.concat_1000_genomes.output.merged_vcf,
        definition   = rules.format_pedigree_definition.output.definition,
        map          = config['ped-sim']['data']['map'],
        #interference = rules.fetch_interference_map.output.intf_map,
        interference = "data/ped-sim/interference_maps/nu_p_campbell.tsv",
        metadata     = "results/meta/pipeline-metadata.yml",
    output:
        fam = "results/00-ped-sim/{POP}-pedigrees-everyone.fam",
        log = "results/00-ped-sim/{POP}-pedigrees.log",
        seg = "results/00-ped-sim/{POP}-pedigrees.seg",
        vcf = temp(pipe("results/00-ped-sim/{POP}-pedigrees.vcf")),
    params:
        output_basename  = "results/00-ped-sim/{POP}-pedigrees",
        error_rate       = config['ped-sim']['params']['error_rate'],
        missingness      = config['ped-sim']['params']['missingness'],
        retain_extra     = config['ped-sim']['params']['retain-extra'],
        seed             = set_ped_sim_seed
    log:   "logs/00-ped-sim/run_ped_sim/{POP}_ped-sim.log"
    conda: "../envs/ped-sim-1.3.yml"
    shell: """
        ped-sim \
        -d {input.definition} \
        -m {input.map} \
        -i {input.vcf} \
        -o {params.output_basename} \
        --intf {input.interference} \
        --miss_rate {params.missingness} \
        --err_rate {params.error_rate} \
        --pseudo_hap 0 \
        --err_hom_rate 0 \
        --fam \
        --keep_phase \
        --founder_ids \
        --mrca \
        --bp \
        --nogz \
        --retain_extra {params.retain_extra} \
        --seed {params.seed} > {log} 2>&1 
    """


rule dopplegang_twins:
    """
    Pick a selected individual within ped-sim's vcf output (we run through our pedigree_codes file using get_twins()
    and search for any self-comparison). Copy the genotype information as a separate sample and create a merged vcf
    To imitate monozygotic twins downstream.

    NOTE: At least 6 bgzip threads is required for it to keep up.

    NOTE: dopplegang-samples is an absolutely barbaric script which merely pastes the requested column(s) at the 
          end of each line. This is MUCH faster than using bcftools merge, but keep in mind that this program:
           - does not update the INFO tags (e.g.: AC, AF, NS, etc.)
           - does not check for genomic coordinate correspondance.
          
          We know this is absolutely fine, because we're only feeding this file to bcftools consensus, wich only
          requires the raw allele information, coordinate, and type of mutation. But keep in mind that this
          output file is 'corrupted' in the sense that the information contained within the INFO column now makes
          absolutely no sense.
    """
    input:
        vcf        = rules.run_ped_sim.output.vcf
    output:
        merged_vcf = "results/00-ped-sim/{POP}-pedigrees-twins-merged.vcf.gz"
    params:
        twins      = get_twins
    log:       "logs/00-ped-sim/dopplegang_twins/dopplegang_twins-{POP}.log"
    benchmark: "benchmarks/{POP}.dopplegang_twins.txt"
    conda:     "../envs/dopplegang-samples.yml"
    threads:   8
    shell: """
        cat {input.vcf} | dopplegang-samples {params.twins} | bgzip --compress-level 9 --threads 6 > {output.merged_vcf}
    """

checkpoint get_samples:
    """
    Extract the ped-sim samples identifiers contained within the VCF. These will become our wildcards
    for the remainder of this simulations pipeline.
    """
    input:
        vcf = expand(rules.dopplegang_twins.output.merged_vcf, POP=config["ped-sim"]["params"]["pop"]),
    output:
        "results/00-ped-sim/sample_names.tsv"
    log: "logs/00-ped-sim/get_samples.log"
    priority: 50
    shell: """
        zcat {input.vcf} | grep '#CHROM' | cut -f 10- > {output} 2> {log}
    """