from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider() # Anonymous 
# https://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz

# ---- Utility functions
def GRUPS_output(wildcards):
    """
    Returns the list of required generation_wise comparisons for GRUPS-rs
    """
    with checkpoints.get_samples.get().output[0].open() as f:
        samples = str.split(f.readline().replace('\n', ''), '\t')
        generations = set([sample.split('_')[0] for sample in samples])
        return directory(expand("results/04-kinship/GRUPS/{generation}", generation=generations))

# ------------------------------------------------------------------------------------------------------------------- #

rule GRUPS_fetch_recombination_map:
    """
    Download samples metadata from the 1000g FTP website
    """
    input:
        tarball = HTTP.remote("http://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz")
    output:
        map     = expand("data/recombination-maps/HapMapII_GRCh37/genetic_map_GRCh37_chr{chr}.txt", chr=range(1, 23)),
        exclude = temp(expand("data/recombination-maps/HapMapII_GRCh37/genetic_map_GRCh37_chr{chr}.txt", chr=["X", "X_par1", "X_par2"])),
        readme  = temp("data/recombination-maps/HapMapII_GRCh37/README.txt")
    params:
        output_dir = lambda wildcards, output: dirname(output.map[0])
    shell: """
        tar -xvzf {input.tarball} -C {params.output_dir}
        rm {input.tarball}
    """

rule GRUPS_generate_fst_set:
    input:
        data    = expand(rules.download_1000_genomes.output.vcf, chr=range(1,23)),
        panel   = rules.fetch_samples_panel.output.panel
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
    log: "logs/04-kinship/GRUPS_generate_fst_set/{params.pedigree_pop}-{params.contam_pop}-GRUPS_generate_fst_set.log"
    conda: "../envs/grups-rs.yml"
    threads: 22
    shell: """
        grups fst --vcf-dir $(dirname {input.data} | uniq)               \
                  --output-dir $(dirname {output.fst} | uniq)            \
                  --pop-subset {params.pedigree_pop} {params.contam_pop} \
                  --panel {input.panel}                                  \
                  --threads {threads}                                    \
                  --verbose > {log} 2>&1
    """

rule run_GRUPS:
    input:
        pileup       = rules.samtools_pileup.output.pileup,
        data         = rules.GRUPS_generate_fst_set.output.fst,
        #panel        = rules.fetch_samples_panel.output.panel,
        recomb_map   = rules.GRUPS_fetch_recombination_map.output.map,
        targets      = config["kinship"]["targets"],
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
        maf          = config["kinship"]["GRUPS"]["maf"]
    log: "logs/04-kinship/run_GRUPS/{generation}-run_GRUPS.log"
    conda: "../envs/grups-rs.yml"
    shell: """
        grups pedigree-sims --pileup {input.pileup}                                           \
                            --data-dir {params.data_dir}                                      \
                            --recomb-dir {params.recomb_dir}                                  \
                            --pedigree {params.pedigree}                                      \
                            --pedigree-pop {params.pedigree_pop}                              \
                            --contam-pop {params.contam_pop}                                  \
                            --min-depth {params.min_depth}                                    \
                            --samples 0-{params.samples}                                      \
                            --sample-names {params.sample_names}                              \
                            --reps {params.reps}                                              \
                            --mode {params.mode}                                              \
                            --output-dir {output.output_dir}                                  \
                            --maf {params.maf}                                                \
                            --min-qual {params.min_quality}                                   \
                            --print-blocks                                                    \
                            --verbose > {log} 2>&1
    """