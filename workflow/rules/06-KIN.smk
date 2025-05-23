localrules: symlink_KINgaroo_input_bams

# ---- Utility functions
def KIN_output(wildcards):
    """
    Returns the list of required generation_wise comparisons from KIN
    """
    # Get the number of expected generations.
    template = "results/04-kinship/KIN/{generation}/{generation}-KIN-results/KIN_results.csv"
    return expand(template, generation = get_generations())

# ------------------------------------------------------------------------------------------------ #

rule symlink_KINgaroo_input_bams:
    input:
        bamlist = get_pileup_input_bams,
    output:
        bamlist  = "results/04-kinship/KIN/{generation}/{generation}.bamlist",
        linkdir  = directory("results/04-kinship/KIN/{generation}/symbams/")
    resources:
        cores    = lambda w, threads: threads
    log:       "logs/04-kinship/KIN/symlink_KINgaroo_input_bams/{generation}.log"
    benchmark: "benchmarks/04-kinship/KIN/symlink_KINgaroo_input_bams/{generation}.tsv"
    conda:     "../envs/coreutils-9.1.yml"
    threads: 1
    shell: """
        mkdir -p {output.linkdir} && ln -srft {output.linkdir} {input.bamlist} 2>  {log}
        basename -a {input.bamlist} | sed 's/.bam$//' > {output.bamlist}       2>> {log}
    """


def parse_KINgaroo_optargs(wildcards):
    correct_contamination = False
    diversity_parameter = config['kinship']['KIN']['diversity-parameter']
    noisy_windows       = config['kinship']['KIN']['noisy-windows']
    optargs = ""
    if diversity_parameter is not None:
        optargs += f"--diversity_parameter_p_0 {diversity_parameter} "
    if noisy_windows is not None:
        optargs += f"--noisy_wins {noisy_windows} "

    return optargs

rule run_KINgaroo:
    """
    Run KINgaroo on a set of pedigree individuals (generation-wise.)
    # @ TODO: Add contamination estimate. 

    # Benchmarks: 
    | depth | max h:m:s | max_rss |
    | ----- | --------- | ------- |
    | 0.01X |           |         |
    | 0.05X | 0:11:35   | 4922.52 |
    | 0.10X | 0:45:00   | 3539    |
    """
    input:
        bamlist        = rules.symlink_KINgaroo_input_bams.output.bamlist,
        linkdir        = rules.symlink_KINgaroo_input_bams.output.linkdir,
        targets        = expand(rules.get_target_panel_intersect.output.targets,
            maf      = config['variant-calling']['maf'],
            superpop = config['variant-calling']['maf-superpop'],
        ),
    output:
        kingaroo_dir      = directory("results/04-kinship/KIN/{generation}/kingaroo"),
        bedfiles          = directory("results/04-kinship/KIN/{generation}/kingaroo/bedfiles"),
        hapProbs          = directory("results/04-kinship/KIN/{generation}/kingaroo/hapProbs"),
        hbd_results       = directory("results/04-kinship/KIN/{generation}/kingaroo/hbd_results"),
        hmm_parameters    = directory("results/04-kinship/KIN/{generation}/kingaroo/hmm_parameters"),
        splitbams         = directory("results/04-kinship/KIN/{generation}/kingaroo/splitbams"),
        goodpairs         = "results/04-kinship/KIN/{generation}/kingaroo/goodpairs.csv",
        overlap           = "results/04-kinship/KIN/{generation}/kingaroo/overlap.csv",
        identical_overlap = "results/04-kinship/KIN/{generation}/kingaroo/identical_overlap.csv",
        diffs_hmm         = "results/04-kinship/KIN/{generation}/kingaroo/input_diffs_hmm.csv",
        hbd_hmm_diffs     = "results/04-kinship/KIN/{generation}/kingaroo/input_hbd_hmm_diffs.csv",
        hbd_hmm_total     = "results/04-kinship/KIN/{generation}/kingaroo/input_hbd_hmm_total.csv",
        total_hmm         = "results/04-kinship/KIN/{generation}/kingaroo/input_total_hmm.csv"
    params:
        interval       = config['kinship']['KIN']['interval'],
        threshold      = config['kinship']['KIN']['p0-threshold'],
        contam_param   = config['kinship']['KIN']['contam-parameter'],
        optargs        = parse_KINgaroo_optargs
    resources:
        runtime        = 600,
        mem_mb         = 6000,
        cores          = lambda w, threads: threads
    retries:   3
    log:       "logs/04-kinship/KIN/run_KINgaroo/{generation}.log"
    benchmark: "benchmarks/04-kinship/KIN/run_KINgaroo/{generation}.tsv"
    conda:     "../envs/kin-3.1.3.yml"
    threads:   8
    shell: """
        CWD=`pwd`
        cd {output.kingaroo_dir}
        KINgaroo \
        --cores {threads} \
        --bamfiles_location $CWD/{input.linkdir} \
        --target_location $CWD/{input.bamlist} \
        --bedfile $CWD/{input.targets} \
        --contam_parameter {params.contam_param} \
        --interval {params.interval} \
        --threshold {params.threshold} \
        {params.optargs} \
        > $CWD/{log} 2>&1
    """


def parse_KIN_optargs(wildcards):
    optargs = ""
    roh_threshold       = config['kinship']['KIN']['roh-threshold']
    diversity_parameter = config['kinship']['KIN']['diversity-parameter']
    if roh_threshold is not None:
        optargs += f"--threshold {roh_threshold} "
    if diversity_parameter is not None:
        optargs += f"--diversity_parameter_p_0 {diversity_parameter} "

    return optargs

rule run_KIN:
    """
    Perform generation-wise kinship estimation using KIN, and the output of KINgaroo.
    See: Popli, D., Peyrégne, S. & Peter, B.M. KIN: a method to infer relatedness from low-coverage
         ancient DNA. Genome Biol 24, 10 (2023). https://doi.org/10.1186/s13059-023-02847-7

    Repo: https://github.com/DivyaratanPopli/Kinship_Inference.git

    # @TODO: Add contamination estimate. 

    # Benchmarks: 
    | depth | max h:m:s | max_rss |
    | ----- | --------- | ------- |
    | 0.01X |           |         |
    | 0.05X | 0:00:08   | 2464    |
    | 0.10X | 0:00:17   | 2121    |
    """
    input:
        bamlist      = rules.symlink_KINgaroo_input_bams.output.bamlist,
        linkdir      = rules.symlink_KINgaroo_input_bams.output.linkdir,
        kingaroo_dir = rules.run_KINgaroo.output.kingaroo_dir
    output:
        kin_results  = directory("results/04-kinship/KIN/{generation}/{generation}-KIN-results"),
        gamma_files  = directory("results/04-kinship/KIN/{generation}/{generation}-KIN-results/gammafiles"),
        lik_files    = directory("results/04-kinship/KIN/{generation}/{generation}-KIN-results/likfiles"),
        res_files    = directory("results/04-kinship/KIN/{generation}/{generation}-KIN-results/resfiles"),
        IBDadded     = "results/04-kinship/KIN/{generation}/{generation}-KIN-results/IBDadded.csv",
        KIN_results  = "results/04-kinship/KIN/{generation}/{generation}-KIN-results/KIN_results.csv",
        likelihoods  = "results/04-kinship/KIN/{generation}/{generation}-KIN-results/relatable_allLikelihoods_fil0.csv"
    params:
        workdir      = lambda w, output: os.path.dirname(output.kin_results),
        interval     = config['kinship']['KIN']['interval'],
        optargs      = parse_KIN_optargs
    resources:
        runtime      = 60,
        mem_mb       = 3000,
        cores        = lambda w, threads: threads
    log:       "logs/04-kinship/KIN/run_KIN/{generation}.log"
    benchmark: "benchmarks/04-kinship/KIN/run_KIN/{generation}.tsv"
    conda:     "../envs/kin-3.1.3.yml"
    threads:   16
    shell: """
        CWD=`pwd`; cd {params.workdir}
        KIN \
        --cores {threads} \
        --input_location $CWD/{input.kingaroo_dir}/ \
        --output_location $CWD/{output.kin_results}/ \
        --target_location $CWD/{input.bamlist} \
        --interval {params.interval} \
        {params.optargs} \
        > $CWD/{log} 2>&1
    """
