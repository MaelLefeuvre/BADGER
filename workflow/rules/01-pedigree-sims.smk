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
configfile: "./config/netrules.yml"

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
        sim_map = config["netrules"]["ped-sim"]['data']["default-map"]
    params:
        #map_dir = rules.fetch_sex_specific_recombination_map.output.map_dir
        map_dir = directory("data/recombination-maps/Refined_genetic_map_b37")
    resources:
        cores = lambda w, threads: threads
    log:     "logs/00-ped-sim/format_sex_specific_gen_map.log"
    conda:   "../envs/coreutils-9.1.yml"
    threads: 1
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
        
    return seed + int(wildcards.rep.replace("ped", ""))

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
        #replicates = config['ped-sim']['replicates']
        replicates = 1
    resources:
        cores = lambda w, threads: threads
    threads: 1
    template_engine: "jinja2"

def get_ped_sim_map_file(wildcards):
    user_defined = config['ped-sim']['data']['map']
    default      = config['netrules']['ped-sim']['data']['default-map']
    return user_defined or default

def get_ped_sim_interference_file(wildcards):
    user_defined = config['ped-sim']['data']['interference'] 
    default      = config['netrules']['ped-sim']['data']['default-interference']
    return user_defined or default

rule run_ped_sim:
    """
    Run pedigree-simulator, using individuals from a specified population of the 1000g project as founder individuals.
    @ TODO: Maybe the fact that we're fetching the interference map is the cause of the FTP os error. ?
    """
    input:
        vcf          = rules.concat_1000_genomes.output.merged_vcf,
        definition   = rules.format_pedigree_definition.output.definition,
        map          = get_ped_sim_map_file,
        interference = get_ped_sim_interference_file,
        metadata     = "results/meta/pipeline-metadata.yml",
    output:
        fam  = temp("results/00-ped-sim/splitted/{rep}/{POP}-pedigrees-everyone.fam"),
        log  = temp("results/00-ped-sim/splitted/{rep}/{POP}-pedigrees.log"),
        seg  = temp("results/00-ped-sim/splitted/{rep}/{POP}-pedigrees.seg"),
        ids  = temp("results/00-ped-sim/splitted/{rep}/{POP}-pedigrees.ids"),
        bp   = temp("results/00-ped-sim/splitted/{rep}/{POP}-pedigrees.bp"),
        mrca = temp("results/00-ped-sim/splitted/{rep}/{POP}-pedigrees.mrca"),
        vcf = temp(pipe("results/00-ped-sim/splitted/{rep}/{POP}-pedigrees.vcf")),
    params:
        output_basename  = "results/00-ped-sim/splitted/{rep}/{POP}-pedigrees",
        error_rate       = config['ped-sim']['params']['error_rate'],
        missingness      = config['ped-sim']['params']['missingness'],
        retain_extra     = config['ped-sim']['params']['retain-extra'],
        seed             = set_ped_sim_seed
    resources:
        runtime = 90,
        cores   = lambda w, threads: threads
    log:       "logs/00-ped-sim/{POP}/splitted/{rep}/run_ped_sim.log"
    benchmark: "benchmarks/00-ped-sim/{POP}/splitted/{rep}/run_ped_sim.tsv"
    conda:     "../envs/ped-sim-gsl-1.4.yml"
    threads: 1
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

    # Replace pedigree replicate in every output file.
    for file in {output.fam} {output.log} {output.seg} {output.seg} {output.ids} {output.bp} {output.mrca}; do
        sed -i "s/ped1/{wildcards.rep}/g" $file > {log} 2>&1
    done
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
        merged_vcf = temp("results/00-ped-sim/splitted/{rep}/{POP}-pedigrees-doppleganged.vcf.gz")
    params:
        twins      = get_twins
    resources:
        runtime = 90,
        mem_mb  = 128,
        cores   = lambda w, threads: threads
    log:       "logs/00-ped-sim/{POP}/splitted/{rep}/dopplegang_twins.log"
    benchmark: "benchmarks/00-ped-sim/{POP}/splitted/{rep}/dopplegang_twins.tsv"
    conda:     "../envs/dopplegang-samples.yml"
    threads:   7
    shell: """
        cat {input.vcf} | sed '/#CHROM/{{s/ped1/{wildcards.rep}/g}}' | dopplegang-samples {params.twins} | bgzip --compress-level 9 --threads 6 > {output.merged_vcf}
    """

rule merge_ped_sim:
    input:
        fam  = expand(rules.run_ped_sim.output.fam,  rep = get_generations(), POP = "{POP}"),
        log  = expand(rules.run_ped_sim.output.log,  rep = get_generations(), POP = "{POP}"),
        seg  = expand(rules.run_ped_sim.output.seg,  rep = get_generations(), POP = "{POP}"),
        ids  = expand(rules.run_ped_sim.output.ids,  rep = get_generations(), POP = "{POP}"),
        bp   = expand(rules.run_ped_sim.output.bp,   rep = get_generations(), POP = "{POP}"),
        mrca = expand(rules.run_ped_sim.output.mrca, rep = get_generations(), POP = "{POP}"),
        vcf  = expand(rules.dopplegang_twins.output.merged_vcf, rep = get_generations(), POP = "{POP}"),
        tbi  = expand(rules.dopplegang_twins.output.merged_vcf + ".tbi", rep = get_generations(), POP = "{POP}")

    output:
        fam  = "results/00-ped-sim/{POP}-pedigrees-everyone.fam",
        log  = "results/00-ped-sim/{POP}-pedigrees.log",
        seg  = "results/00-ped-sim/{POP}-pedigrees.seg",
        ids  = "results/00-ped-sim/{POP}-pedigrees.ids",
        bp   = "results/00-ped-sim/{POP}-pedigrees.bp",
        mrca = "results/00-ped-sim/{POP}-pedigrees.mrca",
        vcf  = "results/00-ped-sim/{POP}-pedigrees-doppleganged-merged.vcf.gz"
    params:
        replicates = config['ped-sim']['replicates']
    wildcard_constraints:
        POP = "[A-Z]{3}"
    log:       "logs/00-ped-sim/{POP}/merge_ped_sim.log"
    benchmark: "benchmarks/00-ped-sim/{POP}/merge_ped_sim.tsv"
    conda: "../envs/bcftools-1.15.yml"
    threads: 16
    shell: """
        cat {input.fam} >> {output.fam}   2>  {log}
        cat {input.log} >> {output.log}   2>> {log}
        cat {input.seg} >> {output.seg}   2>> {log}
        cat {input.ids} >> {output.ids}   2>> {log}
        cat {input.bp}  >> {output.bp}    2>> {log}
        cat {input.mrca} >> {output.mrca} 2>> {log}
        ([[ {params.replicates} > 1 ]] \
          && bcftools merge --threads {threads} -Oz -o {output.vcf} {input.vcf} \
          || cp {input.vcf} {output.vcf}) >> {log} 2>&1
    """

rule plot_ped_sim:
    input:
        fam  = rules.merge_ped_sim.output.fam
    output:
        pdf = report("results/00-ped-sim/{POP}-pedigree-simulations.pdf",
            caption= "../report/01-pedigree-sims/plot-ped-sim.rst",
            category = "00. Pedigree simulations"
        )
    log:       "logs/00-ped-sim/{POP}/plot_ped_sim.log"
    benchmark: "benchmarks/00-ped-sim/{POP}/plot_ped_sim.tsv"
    conda: "../envs/plot-pedigree-sims.yml"
    threads: 1
    shell: """
        Rscript workflow/scripts/plot-pedigree.R {input.fam} {output.pdf} > {log} 2>&1
    """

# ---- Temporary workaround to Snakemake's issue #2668
# See https://github.com/snakemake/snakemake/issues/2668
include: "02-adna-simulations.smk"

checkpoint get_samples:
    """
    Extract the ped-sim samples identifiers contained within the VCF. These will become our wildcards
    for the remainder of this simulations pipeline.
    """
    input:
        _dummy = rules.get_contamination_table.output, # ---- Temporary workaround to Snakemake's issue #2668
        vcf = expand(rules.merge_ped_sim.output.vcf, POP=config["ped-sim"]["params"]["pop"]),
    output:
        "results/00-ped-sim/sample_names.tsv"
    resources:
        cores = lambda w, threads: threads
    log:      "logs/00-ped-sim/get_samples.log"
    conda:    "../envs/bcftools-1.15.yml"
    threads:  3
    priority: 50
    shell: """
        bcftools view -h {input.vcf} | grep '#CHROM' | cut -f 10- > {output} 2> {log}
    """
