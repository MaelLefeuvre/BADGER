import yaml

#module netrules:
#    snakefile: "00-netrules.smk"
#    config: config
#
#use rule download_HapMapII_recombination_map from netrules
#use rule download_1000_genomes from netrules



# ---- Utility functions
def GRUPS_output(wildcards):
    """
    Returns the list of required generation_wise comparisons for GRUPS-rs
    @ TODO: This checkpoint is not needed anymore
    """
    with checkpoints.get_samples.get().output[0].open() as f:
        samples = str.split(f.readline().replace('\n', ''), '\t')
        generations = set([sample.split('_')[0] for sample in samples])
        return directory(expand("results/04-kinship/GRUPS/{generation}", generation=generations))
    #generations = range(1, config['ped-sim']['replicates'])
    #return directory(expand("results/04-kinship/GRUPS/{generation}", generation=generations))

# ------------------------------------------------------------------------------------------------------------------- #

rule GRUPS_generate_fst_set:
    """
    Generate a .fst and .fst.frq dataset from a set of VCF files. Allows for
    generally faster IO processing when performing multiple runs.
    @ TODO: The download directive might be the reason why we have an FTP OS error
    """
    input:
        #data    = expand(rules.download_1000_genomes.output.vcf, chr=range(1,23)),
        #panel   = rules.fetch_samples_panel.output.panel
        data     = expand(
            "data/vcf/1000g-phase3/00-original/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
            chr=range(1,23)
        ),
        
        panel = "data/vcf/1000g-phase3/samples-list/integrated_call_samples_v3.20130502.ALL.panel",
    output:
        fst     = protected(expand(
            "data/grups/fst/g1k-phase3-v5/{ped_pop}-{cont_pop}/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes-{ped_pop}-{cont_pop}.{ext}",
            ped_pop  = config["kinship"]["GRUPS"]["pedigree-pop"],
            cont_pop = config["kinship"]["GRUPS"]["contam-pop"],
            chrom    = range(1,23),
            ext      = ["fst", "fst.frq"]
        ))
    params:
        pedigree_pop = config["kinship"]["GRUPS"]["pedigree-pop"],
        contam_pop   = config["kinship"]["GRUPS"]["contam-pop"]
    log: "logs/04-kinship/GRUPS/GRUPS_generate_fst_set/{params.pedigree_pop}-{params.contam_pop}-GRUPS_generate_fst_set.log"
    conda: "../envs/grups-rs.yml"
    threads: 22
    shell: """
        grups fst \
        --vcf-dir $(dirname {input.data} | uniq) \
        --output-dir $(dirname {output.fst} | uniq) \
        --pop-subset {params.pedigree_pop} {params.contam_pop} \
        --panel {input.panel} \
        --threads {threads} \
        --verbose > {log} 2>&1
    """

def format_seed(wildcards):
    """
    Return a user-defined seed from the config file if it was set. Else, fetch
    the randomly generated backup seed from our metadata file.
    """
    seed = config['kinship']['GRUPS']['seed']
    if seed is None:
        with open(rules.meta.output.metadata) as f:
            metadata = yaml.load(f, Loader=yaml.loader.SafeLoader)
            seed     = metadata['seed']
    
    return {seed}


rule run_GRUPS:
    """
    Run grups on an entire pedigree generation.
    @ TODO: VCF mode is not yet implemented within the pipeline.
            A simple function should do the trick.
    @ TODO: Fetch sample panel might be the cause of our download error.
    """
    input:
        pileup       = rules.samtools_pileup.output.pileup,
        data         = rules.GRUPS_generate_fst_set.output.fst,
        #panel        = rules.fetch_samples_panel.output.panel,
        #recomb_map   = rules.download_HapMapII_recombination_map.output.map,
        recomb_map   = expand("data/recombination-maps/HapMapII_GRCh37/genetic_map_GRCh37_chr{chr}.txt", chr=range(1, 23)),
        targets      = config["kinship"]["targets"],
        metadata     = "results/meta/pipeline-metadata.yml"
    output:
        output_dir   = directory("results/04-kinship/GRUPS/{generation}/"),
        results      = multiext("results/04-kinship/GRUPS/{generation}/{generation}", ".pwd", ".result")
    params:
        sample_names = lambda wildcards: expand(get_samples_ids(wildcards), generation = wildcards.generation),
        samples      = lambda wildcards: len(get_samples_ids(wildcards)) - 1,
        data_dir     = lambda wildcards, input: dirname(input.data[0]),
        recomb_dir   = lambda wildcards, input: dirname(input.recomb_map[0]),
        pedigree     = config["kinship"]["GRUPS"]["pedigree"],
        pedigree_pop = config["kinship"]["GRUPS"]["pedigree-pop"],
        contam_pop   = config["kinship"]["GRUPS"]["contam-pop"],
        reps         = config["kinship"]["GRUPS"]["reps"],
        mode         = config["kinship"]["GRUPS"]["mode"],
        min_depth    = config["kinship"]["GRUPS"]["min-depth"],
        min_quality  = config["kinship"]["GRUPS"]["min-qual"], 
        maf          = config["kinship"]["GRUPS"]["maf"],
        q_error_rate = config["kinship"]["GRUPS"]["seq-error-rate"],
        seed         = format_seed
    log: "logs/04-kinship/GRUPS/run_GRUPS/{generation}-run_GRUPS.log"
    conda: "../envs/grups-rs.yml"
    shell: """
        grups pedigree-sims \
        --pileup {input.pileup} \
        --data-dir {params.data_dir} \
        --recomb-dir {params.recomb_dir} \
        --pedigree {params.pedigree} \
        --pedigree-pop {params.pedigree_pop} \
        --contam-pop {params.contam_pop} \
        --min-depth {params.min_depth} \
        --samples 0-{params.samples} \
        --sample-names {params.sample_names} \
        --reps {params.reps} \
        --mode {params.mode} \
        --output-dir {output.output_dir} \
        --maf {params.maf} \
        --min-qual {params.min_quality} \
        --seq-error-rate {params.q_error_rate} \
        --seed {params.seed} \
        --print-blocks \
        --ignore-dels \
        --verbose > {log} 2>&1
    """
