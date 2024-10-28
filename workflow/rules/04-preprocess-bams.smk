configfile: "./config/config.yml"

# ------------------------------------------------------------------------------------------------------------------- #
# ---- Generic / Utility rules
rule samtools_index:
    """
    Index a generic .bam file
    """
    input:
        bam   = "{directory}/{bam}"
    output:
        bai   = "{directory}/{bam}.bai"
    resources:
        cores = lambda w, threads: threads
    log:       "logs/generics/{directory}/samtools_index-{bam}.log"
    benchmark: "benchmarks/generics/{directory}/samtools_index-{bam}.log"
    conda:     "../envs/samtools-1.15.yml"
    threads:   4
    shell: """
        samtools index -@ {threads} {input.bam}
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 01. Perform read-length and BQ quality filtration.

rule samtools_filter_unmapped:
    """
    Apply QC filtration after raw alignment. Three filters are applied concurrently:
    - (user-defined) minimum base quality
    - (user-defined) minimum length
    - (constant)     remove unmapped

    | depth | max_rss |
    | ----- | ------- | 
    | 0.01X |         |
    | 0.05X | 7.59    |
    | 0.10X | 30.89   |
    """
    input:
        sam        = assign_aligner_algorithm,
        reference  = config["reference"],
        bwt        = rules.index_reference_genome.output.bwt
    output:
        bam        = temp("results/02-preprocess/03-filter/{sample}/{sample}.bam")
    params:
        min_MQ     = config["preprocess"]["filter"]["min-MQ"],
        min_length = config["preprocess"]["filter"]["min-length"],
    resources:
        runtime       = 10,
        mem_mb        = 128,
        cores         = lambda w, threads: threads
    log:       "logs/02-preprocess/03-filter/samtools_filter_unmapped/{sample}.log"
    benchmark: "benchmarks/02-preprocess/03-filter/samtools_filter_unmapped/{sample}.tsv"
    conda:     "../envs/samtools-1.15.yml"
    priority:  6
    threads:   1
    shell: """
        samtools view \
        --threads {threads} \
        --reference {input.reference} \
        -q {params.min_MQ} \
        -e 'length(seq)>{params.min_length}' \
        -F4 -Sb \
        {input.sam} > {output.bam} 2> {log}
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 02. Sort BAM file. 


rule samtools_sort:
    """
    Apply coordinate sorting on a BAM file.

    # Benchmarks:
    | depth | max_rss |
    | ----- | ------- | 
    | 0.01X |         |
    | 0.05X | 898     |
    | 0.10X | 894     |
    """
    
    input:
        bam       = rules.samtools_filter_unmapped.output.bam,
        reference = config["reference"]
    output:
        bam       = temp("results/02-preprocess/04-sort/{sample}/{sample}.srt.bam")
    resources:
        runtime       = 10,
        mem_mb        = 1024,
        cores         = lambda w, threads: threads
    log:       "logs/02-preprocess/04-sort/samtools_sort/{sample}.log"
    benchmark: "benchmarks/02-preprocess/04-sort/samtools_sort/{sample}.tsv"
    conda:     "../envs/samtools-1.15.yml"
    threads:   1
    priority:  7
    shell: """
        samtools sort -@ {threads} --reference {input.reference} --output-fmt BAM -o {output.bam} {input.bam} 2> {log}
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 03. Remove duplicates using Picard, Dedup or Samtools.

rule picard_rmdup:
    """
    Remove PCR duplicates from a BAM file using picard.

    # Benchmarks: 
    | depth | max_rss |
    | ----- | ------- |
    | 0.01X |         |
    | 0.05X | 1612    |
    | 0.10X | 1828    |
    """
    input:
        bam     = rules.samtools_sort.output.bam,
    output:
        bam     = "results/02-preprocess/05-dedup/picard/{sample}/{sample}.srt.rmdup.bam",
        metrics = "results/02-preprocess/05-dedup/picard/{sample}/{sample}.rmdup.metrics.txt"
    resources:
        runtime = 10,
        mem_mb  = 2048,
        tmpdir  = config["tempdir"],
        cores   = lambda w, threads: threads
    log:       "logs/02-preprocess/05-dedup/picard/picard_rmdup/{sample}.log"
    benchmark: "benchmarks/02-preprocess/05-dedup/picard/picard_rmdup/{sample}.tsv"
    conda:     "../envs/picard-2.27.4.yml"
    threads:   1
    shell: """
        picard MarkDuplicates \
        -I {input.bam} \
        -O {output.bam} \
        -M {output.metrics} \
        --ASSUME_SORT_ORDER coordinate \
        --REMOVE_DUPLICATES true \
        --VALIDATION_STRINGENCY LENIENT \
        --TMP_DIR {resources.tmpdir} 2> {log}
    """

rule apeltzer_dedup:
    """
    Remove PCR duplicates from a BAM file using Alexander Peltzer's 'dedup'
    tool (i.e. from nf-core EAGER pipeline)
    @ TODO: FIXME: not yet fully implemented: AdapterRemovalFixPrefix
            should be applied on fastq files, right after AdapterRemoval not bams.
    """
    input:
        bam   = rules.samtools_sort.output.bam
    output:
        bam   = temp("results/02-preprocess/05-dedup/dedup/{sample}/{sample}.srt.rmdup.bam")
    resources:
        cores = lambda w, threads: threads
    log:       "logs/02-preprocess/05-dedup/dedup/apeltzer_dedup/{sample}.log"
    benchmark: "benchmarks/02-preprocess/05-dedup/dedup/apeltzer_dedup/{sample}.tsv"
    conda:     "../envs/dedup-0.12.8.yml"
    threads:   1
    shell: """
        AdapterRemovalFixPrefix {input.bam} | dedup --output {output.bam} 2> {log}
    """


rule samtools_rmdup:
    """
    Remove PCR duplicates from a BAM file using samtool's rmdup
    """
    input:
        bam    = rules.samtools_sort.output.bam
    output:
        bam    = temp("results/02-preprocess/05-dedup/samtools/{sample}/{sample}.srt.rmdup.bam")
    params:
        tmpdir = lambda wildcards, output: splitext(output.bam)[0]
    resources:
        cores  = lambda w, threads: threads
    log:       "logs/02-preprocess/05-dedup/samtools/samtools_rmdup/{sample}.log"
    benchmark: "benchmarks/02-preprocess/05-dedup/samtools/samtools_rmdup/{sample}.tsv"
    conda:     "../envs/samtools-1.15.yml"
    threads:   4
    shell: """
    samtools index -@ {threads} {input.bam}
    samtools collate -O {input.bam} \
    | samtools fixmate -m - - \
    | samtools sort - \
    | samtools markdup -r -T {params.tmpdir} -s - {output.bam} 2> {log}
    """


# ------------------------------------------------------------------------------------------------------------------- #
# ---- 04. Perform PMD Base Quality rescaling using MapDamage 

def define_dedup_input_bam(wildcards):
    """
    Choose the appropriate input for run_pmdtools / run_mapdamage, based on
    which duplicate removal tool was requested by the user.
    """
    match config['preprocess']['dedup']['method']:
        case "picard":
            return rules.picard_rmdup.output.bam
        case "dedup":
            return rules.apeltzer_dedup.output.bam
        case "samtools":
            return rules.samtools_rmdup.output.bam

    raise RuntimeError(f"Invalid config value for dedup method '{config['preprocess']['dedup']['method']}'.")


rule run_pmdtools:
    input:
        bam        = define_dedup_input_bam,
        bai        = lambda wildcards: define_dedup_input_bam(wildcards) + ".bai",
        reference  = config['reference'],
    output:
        bam        = "results/02-preprocess/06-pmdtools/{sample}/{sample}.srt.rmdup.filtercontam.bam"
    params:
        threshold  = config['preprocess']['pmd-rescaling']['pmdtools']['threshold'],
        mask_deams = config['preprocess']['pmd-rescaling']['pmdtools']['mask-terminal-deams'],
    resources:
        cores      = lambda w, threads: threads
    log:       "logs/02-preprocess/06-pmdtools/run_pmdtools/{sample}.log"
    benchmark: "benchmarks/02-preprocess/06-pmdtools/run_pmdtools/{sample}.tsv"
    conda:     "../envs/pmdtools-0.60.yml"
    threads:   3
    shell: """
        ### --threshold {params.threshold} 
        samtools calmd {input.bam} {input.reference} \
        | pmdtools \
          --header \
          --adjustbaseq \
          --stats \
          --threshold -9999 \
          --upperthreshold 999999 \
          --maskterminaldeaminations {params.mask_deams} \
        | samtools view -OBAM > {output.bam} 2> {log}
    """


def define_mapdamage_seed(wildcards):
    """
    Add additional seeding arguments for MapDamage if the user-requested for 
    downsampling. If the downsample seed was provided within the config file,
    provide the program with said value. If not, fetch the randomly generated
    backup seed from our metadata file.
    """
    downsample_n = config['preprocess']['pmd-rescaling']['map-damage']['downsample']

    if downsample_n is None:
        return ""

    seed = config['preprocess']['pmd-rescaling']['map-damage']['downsample-seed']
    if seed is None:
        with open(rules.meta.output.metadata) as f:
            metadata = yaml.load(f, Loader=yaml.loader.SafeLoader)
            seed     = metadata['seed']

    return f"--downsample {downsample_n} --downsample-seed {seed}"


rule run_mapdamage:
    """
    Apply PMD base quality score recalibration on a bam file using MapDamage.

    # Benchmarks: 
    | depth | max h:m:s | max_rss |
    | ----- | --------- | ------- |
    | 0.01X |           |         |
    | 0.05X | 0:20:21   | 537     |
    | 0.10X |           | 500     |
    """
    input:
        bam       = define_dedup_input_bam,
        bai       = lambda wildcards: define_dedup_input_bam(wildcards) + ".bai",
        reference = config["reference"],
        metadata  = "results/meta/pipeline-metadata.yml"
    output:
        bam              = [] if config['preprocess']['pmd-rescaling']['apply-masking'] else "results/02-preprocess/06-mapdamage/{sample}/{sample}.srt.rmdup.rescaled.bam",
        g2a_freq         = "results/02-preprocess/06-mapdamage/{sample}/3pGtoA_freq.txt",
        c2t_freq         = "results/02-preprocess/06-mapdamage/{sample}/5pCtoT_freq.txt",
        dnacomp_genome   = "results/02-preprocess/06-mapdamage/{sample}/dnacomp_genome.csv",
        dnacomp          = "results/02-preprocess/06-mapdamage/{sample}/dnacomp.txt",
        lgdistribution   = "results/02-preprocess/06-mapdamage/{sample}/lgdistribution.txt",
        misincorporation = "results/02-preprocess/06-mapdamage/{sample}/misincorporation.txt",
        runtime_log      = "results/02-preprocess/06-mapdamage/{sample}/Runtime_log.txt",
        mcmc_prob        = "results/02-preprocess/06-mapdamage/{sample}/Stats_out_MCMC_correct_prob.csv",
        mcmc_iter        = "results/02-preprocess/06-mapdamage/{sample}/Stats_out_MCMC_iter.csv",
        mcmc_iter_summ   = "results/02-preprocess/06-mapdamage/{sample}/Stats_out_MCMC_iter_summ_stat.csv",
        frag_plot        = report("results/02-preprocess/06-mapdamage/{sample}/Fragmisincorporation_plot.pdf",
            caption     = "../report/04-preprocess-bams/run_mapdamage/frag-plot.rst",
            category    = "03. BAM Preprocessing",
            subcategory = "06. MapDamage",
            labels      = {"sample": "{sample}", "figure": "Fragment misincorporation"}
        ),
        length_plot      = report("results/02-preprocess/06-mapdamage/{sample}/Length_plot.pdf",
            caption     = "../report/04-preprocess-bams/run_mapdamage/length-plot.rst",
            category    = "03. BAM Preprocessing",
            subcategory = "06. MapDamage",
            labels      = {"sample": "{sample}", "figure": "Fragment Length distribution"}
        ),
        mcmc_post_pred   = report("results/02-preprocess/06-mapdamage/{sample}/Stats_out_MCMC_post_pred.pdf",
            caption     = "../report/04-preprocess-bams/run_mapdamage/mcmc-post-pred.rst",
            category    = "03. BAM Preprocessing",
            subcategory = "06. MapDamage",
            labels      = {"sample": "{sample}", "figure": "MCMC Posterior Prediction intervals"}
        ),
        mcmc_hist        = report("results/02-preprocess/06-mapdamage/{sample}/Stats_out_MCMC_hist.pdf",
            caption     = "../report/04-preprocess-bams/run_mapdamage/mcmc-hist.rst",
            category    = "03. BAM Preprocessing",
            subcategory = "06. MapDamage",
            labels      = {"sample": "{sample}", "figure": "MCMC posterior distributions histogram"}
        ),
        mcmc_trace       = report("results/02-preprocess/06-mapdamage/{sample}/Stats_out_MCMC_trace.pdf",
            caption     = "../report/04-preprocess-bams/run_mapdamage/mcmc-trace.rst",
            category    = "03. BAM Preprocessing",
            subcategory = "06. MapDamage",
            labels      = {"sample": "{sample}", "figure": "MCMC parameter trace plots"}
        )
    params:
        downsample_seed  = define_mapdamage_seed,
        rescale          = "" if config['preprocess']['pmd-rescaling']['apply-masking'] else "--rescale",
        rescale_out_flag = "" if config['preprocess']['pmd-rescaling']['apply-masking'] else "--rescale-out"
    resources:
        runtime = 60,
        mem_mb  = 1024,
        cores   = lambda w, threads: threads
    log:       "logs/02-preprocess/06-mapdamage/run_mapdamage/{sample}.log"
    benchmark: "benchmarks/02-preprocess/06-mapdamage/run_mapdamage/{sample}.tsv"
    conda:     "../envs/mapdamage-2.2.1.yml"
    threads: 1
    priority: 10
    shell: """
        mapDamage \
        -i {input.bam} \
        -r {input.reference} \
        --folder $(dirname {output.misincorporation}) \
        {params.downsample_seed} \
        {params.rescale} \
        {params.rescale_out_flag} {output.bam} \
        --verbose > {log} 2>&1 
    """


def define_masking_input_bam(wildcards):
    rescaler = config['preprocess']['pmd-rescaling']['rescaler']
    if rescaler is None:
        return define_dedup_input_bam(wildcards)
    else:
        return define_rescale_input_bam(wildcards)
    raise RuntimeError("Failed to define a proper input bam file for pmd-mask")


rule run_pmd_mask:
    input:
        bam              = define_masking_input_bam,
        bai              = lambda wildcards: define_masking_input_bam(wildcards) + ".bai",
        reference        = config["reference"],
        misincorporation = rules.run_mapdamage.output.misincorporation
    output:
        bam     = "results/02-preprocess/06-pmd-mask/{sample}/{sample}.pmdmasked.bam",
        metrics = "results/02-preprocess/06-pmd-mask/{sample}/{sample}.pmdmasked.metrics",
    params:
        threshold = config['preprocess']['pmd-rescaling']['pmd-mask']['threshold']
    log:       "logs/02-preprocess/06-pmd-mask/run_pmd_mask/{sample}.log"
    benchmark: "benchmarks/02-preprocess/06-pmd-mask/run_pmd_mask/{sample}.tsv" 
    conda:     "../envs/pmd-mask-0.3.2.yml"
    threads:   8
    shell: """
        pmd-mask -@ {threads} -b {input.bam} -f {input.reference} -m {input.misincorporation} --threshold {params.threshold} -M {output.metrics} -Ob -o {output.bam} --verbose > {log} 2>&1
    """




def get_contaminants(wildcards):
    """
    Returns a list of random ID(s) from a list of potential individuals to use as contamination.
    This will give out a single contaminating individual for each pedigree replicate.
    @ TODO: Maybe adding this info to the global metadata.yml file would be a good idea ?
    """

    # Prevent unneccessary rule re-run trigger events: If the contamination table already exists,
    # simply print out its contents and leave...
    if Path(rules.get_contamination_table.output.cont_table).exists():
        with open(rules.get_contamination_table.output.cont_table, "r") as contamination_table:
            samples = [line.strip("\n").split()[1] for line in contamination_table.readlines()]
            return samples

    # else, open the pedigree generation file and 1000g panel definition file and return a list
    # of random samples.
    # @ TODO: This might be where the FTP OS error is located ?
    #with open(rules.run_ped_sim.input.definition) as f, open(rules.fetch_samples_panel.output.panel) as samples:

    gen_no              = config['ped-sim']['replicates']
    contam_samples_file = rules.get_target_pop_samples.output.target_list.format(
        POP=config["gargammel"]["params"]["contam-pop"]
    )
    with open(contam_samples_file) as samples:
        # Fetch sample-IDs matching the user-provided contaminating population tag ('EUR', 'AFR', 'YRI', etc...)
        cont_pop_tag = config["gargammel"]["params"]["contam-pop"]
        contaminants = [sample.strip("\n").split("\t")[0] for sample in samples.readlines()]
        # Return a random list of size gen_no (one contaminating individual per replicate).
        return random.sample(contaminants, gen_no)


# ------------------------------------------------------------------------------------------------ #
# ---- Define the output of this smk. Optionally output QC statistics with FastQC + MultiQC

def define_rescale_input_bam(wildcards):
    rescaler = config['preprocess']['pmd-rescaling']['rescaler']
    match rescaler:
        case "mapdamage":
            return "results/02-preprocess/06-mapdamage/{sample}/{sample}.srt.rmdup.rescaled.bam"
        case "pmdtools":
            return "results/02-preprocess/06-pmdtools/{sample}/{sample}.srt.rmdup.filtercontam.bam"
        case other:
            raise RuntimeError(f'Invalid rescaler value "{rescaler}"')


def get_pileup_input_bams(wildcards, print_log = False, logfile = sys.stderr, filter_samples = True):
    """
    Define the appropriate input bam for the variant caller, based on which 
    PMD-rescaling method was requested by the user.
    """
    # Run through the initial samples files and extract pedigree ids 
    samples_ids = get_samples_ids_filtered(wildcards) if filter_samples else get_samples_ids(wildcards)

    # If masking is required, delegate input definition to the appropriate rule.
    apply_masking = config['preprocess']['pmd-rescaling']['apply-masking']
    if apply_masking:
        if print_log:
            print("Applying pmd-mask for variant calling.", file=logfile)
        return expand(rules.run_pmd_mask.output.bam, sample = samples_ids)


    # Return a list of input bam files for pileup
    rescaler = config['preprocess']['pmd-rescaling']['rescaler']
    if rescaler is None:
        if print_log:
            print("WARNING: Skipping PMD Rescaling for variant calling!", file=logfile)
        return expand(define_dedup_input_bam(wildcards), sample = samples_ids)
    else:
        if print_log:
            print(f"NOTE: Applying {rescaler} for variant calling.", file=logfile)
        return expand(define_rescale_input_bam(wildcards), sample = samples_ids)        
    
    raise RuntimeError(f'Invalid rescaler value "{rescaler}"')


ruleorder: multiqc_preprocess_bams > multiqc
use rule multiqc as multiqc_preprocess_bams with:
    input:
        required_files = lambda w: temp(expand(
            expand(rules.run_fastqc_bam.output.data, zip,
                directory=[dirname(file) for file in get_pileup_input_bams(w)],
                file=[basename(splitext(file)[0]) for file in get_pileup_input_bams(w)]
            ),
            generation = get_generations()
        )),
    output: 
        html = report("results/00-qc/04-preprocess-bams/multiqc-report.html",
            caption     = "../report/04-preprocess-bams/multiqc-preprocess.rst",
            category    = "03. BAM Preprocessing",
            subcategory = "MultiQC"
        )
    params:
        extra_dirs = "results/02-preprocess",
        extra_args = expand("--ignore results/02-preprocess/{subdir}", subdir = ["00-raw", "01-adapter_removal", "02-align"])


rule goleft_covstats:
    input:
        bam = lambda w: expand(get_pileup_input_bams(w),
            sample = get_samples_ids(w),
            generation="{generation}"
        ),
        bai = lambda w: expand([bam + ".bai" for bam in get_pileup_input_bams(w)],
            sample = get_samples_ids(w),
            generation="{generation}"
        )
    output:
        covstats = "results/00-qc/04-preprocess-bams/coveragestats/{generation}-covstats.tsv"
    log:       "logs/00-qc/04-preprocess-bams/{generation}-goleft_covstats.log"
    benchmark: "benchmarks/00-qc/04-preprocess-bams/{generation}-goleft_covstats.tsv"
    conda: "../envs/goleft-0.2.4.yml"
    shell: """
        goleft covstats {input.bam} > {output.covstats} 2> {log}
    """

rule plot_covstats:
    input:
        covstats = rules.goleft_covstats.output.covstats
    output:
        plot = report("results/00-qc/04-preprocess-bams/coveragestats/{generation}-covstats.svg",
            category    = "03. BAM Preprocessing",
            subcategory = "Coverage",
            labels      = {"replicate": "{generation}", "figure": "goleft-covstats-plot"}

        )
    params:
        expected = config['gargammel']['coverage']
    conda: "../envs/plot-covstats.yml"
    shell: """
        workflow/scripts/plot-covstats.py --input {input.covstats} --output {output.plot} --expect {params.expected}
    """

