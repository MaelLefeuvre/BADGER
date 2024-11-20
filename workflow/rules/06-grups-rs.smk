import yaml

# ---- Utility functions
def GRUPS_output(wildcards):
    """
    Returns the list of required generation_wise comparisons for GRUPS-rs
    """
    # Get the number of expected generations.
    gen_no   = config['ped-sim']['replicates']
    template = expand("results/04-kinship/GRUPS/{{generation}}/{{generation}}.{ext}",
        ext=["pwd", "result", "probs"]
    )
    return expand(template, generation = get_generations())


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

# ------------------------------------------------------------------------------------------------ #

rule GRUPS_generate_fst_set:
    """
    Generate a .fst and .fst.frq dataset from a set of VCF files. Allows for
    generally faster IO processing when performing multiple runs.
    # @TODO:
    - Fetch sample panel might be the cause of our FTP error when fetching data (?)
    """
    input:
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
    resources:
        cores        = lambda w, threads: threads
    log:       "logs/04-kinship/GRUPS/GRUPS_generate_fst_set/GRUPS_generate_fst_set.log"
    benchmark: "benchmarks/04-kinship/GRUPS/GRUPS_generate_fst_set/GRUPS_generate_fst_set.tsv"
    conda:     "../envs/grups-rs-0.3.2.yml"
    threads:   22
    shell: """
        grups-rs fst \
        --vcf-dir $(dirname {input.data} | uniq) \
        --output-dir $(dirname {output.fst} | uniq) \
        --pop-subset {params.pedigree_pop} {params.contam_pop} \
        --panel {input.panel} \
        --threads {threads} \
        --verbose > {log} 2>&1
    """


rule run_GRUPS:
    """
    Run grups-rs on an entire pedigree generation to perform kinship estimation.

    # See:
    - Martin MD, Jay F, Castellano S, Slatkin M. Determination of genetic relatedness from 
      low-coverage human genome sequences using pedigree simulations. 
      Mol Ecol. 2017;26(16):4145-4157. https://doi.org/10.1111/mec.14188
    - Lefeuvre M, Martin M, Jay F, Marsolier M, Bon C. GRUPS-rs, a high-performance ancient
      DNA genetic relatedness estimation software relying on pedigree simulations.
      Hum Popul Genet Genom 2024; 4(1):0001. https://doi.org/10.47248/hpgg2404010001
         
    # Github: 
    - (py-grups) https://github.com/sameoldmike/grups/
    - (grups-rs) https://github.com/MaelLefeuvre/grups-rs.git

    # @TODO: 
    - VCF mode is not yet implemented within the pipeline. A simple function should do the trick.

    # Benchmarks: 
    | depth | max h:m:s | max_rss |
    | ----- | --------- | ------- |
    | 0.01X |           |         |
    | 0.05X | 0:09:42   | 2013    |
    | 0.10X | 0:44:15   | 2157    |
    """
    input:
        samples_def  = rules.get_samples.output, 
        pileup       = rules.samtools_pileup.output.pileup,
        data         = rules.GRUPS_generate_fst_set.output.fst,
        recomb_map   = expand("data/recombination-maps/HapMapII_GRCh37/genetic_map_GRCh37_chr{chr}.txt", chr=range(1, 23)),
        targets      = get_snp_targets(ext=".snp"),
        metadata     = "results/meta/pipeline-metadata.yml"
    output:
        results      = multiext("results/04-kinship/GRUPS/{generation}/{generation}", ".pwd", ".result", ".probs")
    params:
        output_dir     = "results/04-kinship/GRUPS/{generation}/",
        sample_names   = lambda wildcards: expand(get_samples_ids_filtered(wildcards), generation = wildcards.generation),
        samples        = lambda wildcards: len(get_samples_ids_filtered(wildcards)) - 1,
        data_dir       = lambda wildcards, input: dirname(input.data[0]),
        recomb_dir     = lambda wildcards, input: dirname(input.recomb_map[0]),
        pedigree       = config["kinship"]["GRUPS"]["pedigree"],
        pedigree_pop   = config["kinship"]["GRUPS"]["pedigree-pop"],
        contam_pop     = config["kinship"]["GRUPS"]["contam-pop"],
        reps           = config["kinship"]["GRUPS"]["reps"],
        mode           = config["kinship"]["GRUPS"]["mode"],
        min_depth      = config["kinship"]["GRUPS"]["min-depth"],
        min_quality    = config["kinship"]["GRUPS"]["min-qual"], 
        maf            = config["kinship"]["GRUPS"]["maf"],
        seq_error_rate = config["kinship"]["GRUPS"]["seq-error-rate"],
        seed           = format_seed
    resources:
        runtime = 60,
        mem_mb  = 2048,
        cores   = lambda w, threads: threads
    log:       "logs/04-kinship/GRUPS/run_GRUPS/{generation}-run_GRUPS.log"
    benchmark: "benchmarks/04-kinship/GRUPS/run_GRUPS/{generation}-run_GRUPS.tsv"
    conda:     "../envs/grups-rs-0.3.2.yml"
    threads:   1
    shell: """
        grups-rs pedigree-sims \
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
        --output-dir {params.output_dir} \
        --maf {params.maf} \
        --min-qual {params.min_quality} \
        --seq-error-rate {params.seq_error_rate} \
        --seed {params.seed} \
        --overwrite \
        --verbose > {log} 2>&1
    """
