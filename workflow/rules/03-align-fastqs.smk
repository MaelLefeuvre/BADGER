import yaml

configfile: "./config/config.yml"


# ----------------------------------------------------------------------------------------------- #
# ---- 00. Clip Adapters.

def assign_adapter_removal_seed(wildcards):
    """
    Return a user-defined seed from the config file if it was set. Else, fetch
    the randomly generated backup seed from our metadata file.
    """
    seed = config['preprocess']['trimming']['seed']
    if seed is None:
        with open(rules.meta.output.metadata) as f:
            metadata = yaml.load(f, Loader=yaml.loader.SafeLoader)
            seed     = metadata['seed']

    return seed


rule adapter_removal_pe:
    """
    Perform Adapter Trimming for Illumina Paired-End sequuencing data. Contrary to some other
    workflows, we don't output a combined fq file, to allow specific mapping using bwa aln.

    # Benchmarks:
    | depth | h:m:s   | max_rss |
    | ----- | ------- | ------- | 
    | 0.01X |         |         |
    | 0.05X | 0:05:32 | 41.30   |
    | 0.10X | 0:03:31 | 23      |
    """
    input:
        forwd       = "results/02-preprocess/00-raw/{sample}_s1.fq.gz",
        revrs       = "results/02-preprocess/00-raw/{sample}_s2.fq.gz",
        metadata    = "results/meta/pipeline-metadata.yml",
    output:
        trimmed     = temp("results/02-preprocess/01-adapter_removal/{sample}/{sample}.collapsed.gz"),
        truncated   = temp("results/02-preprocess/01-adapter_removal/{sample}/{sample}.collapsed.truncated.gz"),
        discarded   = temp("results/02-preprocess/01-adapter_removal/{sample}/{sample}.discarded.gz"),
        pair1       = temp("results/02-preprocess/01-adapter_removal/{sample}/{sample}.pair1.truncated.gz"),
        pair2       = temp("results/02-preprocess/01-adapter_removal/{sample}/{sample}.pair2.truncated.gz"),
        singleton   = temp("results/02-preprocess/01-adapter_removal/{sample}/{sample}.singleton.truncated.gz")
    params:
        base_name   = "results/02-preprocess/01-adapter_removal/{sample}/{sample}",
        min_overlap = config["preprocess"]["trimming"]["min-overlap"],
        min_length  = config["preprocess"]["trimming"]["min-overlap"],
        min_quality = config["preprocess"]["trimming"]["min-quality"],
        quality_max = config["preprocess"]["trimming"]["qualitymax"],
        seed        = assign_adapter_removal_seed
    resources:
        runtime     = 10,
        mem_mb      = 128,
        cores       = lambda w, threads: threads
    log:       "logs/02-preprocess/01-adapter_removal/adapter_removal_pe/{sample}.log"
    benchmark: "benchmarks/02-preprocess/01-adapter_removal/adapter_removal_pe/{sample}-bench.tsv"
    conda:     "../envs/adapterremoval-2.3.3.yml"
    priority: 4
    threads:  4
    shell: """
        AdapterRemoval \
        --threads {threads} \
        --file1 {input.forwd} \
        --file2 {input.revrs} \
        --basename {params.base_name} \
        --minlength {params.min_length} \
        --minquality {params.min_quality} \
        --qualitymax {params.quality_max} \
        --minadapteroverlap {params.min_overlap} \
        --seed {params.seed} \
        --collapse \
        --gzip 2> {log}
    """


# ----------------------------------------------------------------------------------------------- #
# ---- 01-A. Align to genome using BWA aln

def align_thread_assign(wildcards):
    """
    We know that most paired alignments will be collapsed => More threads should be allocated for
    that file.
    """
    if wildcards.extension == "collapsed":
        return 16
    else:
        return 1


rule sam_to_tmp_bam:
    """
    Convert a samfile to a **temporary** bam file.
    """
    input:
        sam = "{directory}/{file}.sam"
    output:
        bam = temp("{directory}/{file}.tmp-bam")
    resources:
        cores = lambda w, threads: threads
    log:       "logs/generics/{directory}/sam_to_tmp_bam-{file}.log"
    benchmark: "benchmarks/generics/{directory}/sam_to_tmp_bam-{file}.tsv"
    conda:     "../envs/samtools-1.15.yml"
    group:     f"bwa{config['preprocess']['bwa']['aligner']}"
    threads: 1
    shell: """
        samtools view -@ {threads} -OBAM {input.sam} > {output.bam}
    """

def assign_bwa_aln_mem(wildcards):
    match wildcards.extension:
        case "collapsed":
            return 7000
        case "pair1.truncated" | "pair2.truncated":
            return 3400
        case "collapsed.truncated" | "singleton.truncated":
            return 3000
    raise RuntimeError(f"Invalid wildcard for bwa aln extension {wildcards.extension}")
        

rule bwa_aln:
    """
    Use bwa aln to prepare alignment and generate suffix array. Gives off very robust results, at 
    the cost of an impoverished runtime performance compared to bwa mem.
    See:    Oliva, A., Tobler, R., Llamas, B. and Souilmi, Y. (2021), Additional evaluations show that specific
            BWA-aln settings still outperform BWA-mem for ancient DNA data alignment. Ecol Evol, 11: 18743-18748. 
            https://doi.org/10.1002/ece3.8297

    extension: '.collapsed' || '.collapsed.truncated' || '.pair1.truncated' || '.pair2.truncated' || '.singleton.truncated'

    # Benchmarks max RSS:
    | depth | .collapsed | .collapsed.truncated | .pair{1,2} | .singleton.truncated |
    | ----- | ---------- | -------------------- | ---------- | -------------------- |
    | 0.01X | 6911       | 2964                 | 3211       | 2964                 |
    | 0.05X | 7000       | 2964                 | 3360       | 2964                 |
    | 0.10X | 7200       | 2960                 | 3370       | 2960                 |

    # Benchmarks h:m:s:
    | depth | .collapsed | .collapsed.truncated | .pair{1,2} | .singleton.truncated |
    | ----- | ---------- | -------------------- | ---------- | -------------------- |
    | 0.01X |            |                      |            |                      |
    | 0.05X | 0:25:22    | 0:02:53              | 0.29:02    | 0:02:48              |
    | 0.10X |            |                      |            |                      |
    """
    input:
        trimmed       = "results/02-preprocess/01-adapter_removal/{sample}/{sample}.{extension}.gz",
        reference     = ReferenceGenome.get_path(),
        bwt           = rules.index_reference_genome.output.bwt
    output:
        sai           = temp(pipe("results/02-preprocess/02-align/{sample}/{sample}.bwaaln.{extension}.sai"))
    params:
        seed_length   = config['preprocess']['bwa']['bwa-aln']['seed-length'],
        max_open_gap  = config['preprocess']['bwa']['bwa-aln']['max-open-gap'],
        missing_prob  = config['preprocess']['bwa']['bwa-aln']['max-miss-prob'],
        max_seed_diff = config['preprocess']['bwa']['bwa-aln']['max-seed-diff']
    resources:
        runtime       = 60,
        mem_mb        = assign_bwa_aln_mem,
        cores         = lambda w, threads: threads
    log:       "logs/02-preprocess/02-align/bwa_aln/{extension}/{sample}.log"
    benchmark: "benchmarks/02-preprocess/02-align/bwa_aln/{extension}/{sample}-bench.tsv"
    conda:     "../envs/bwa-0.7.17.yml"
    #group:     "bwaaln"
    threads:   align_thread_assign
    priority:  50
    shell: """
        bwa aln \
        {input.reference} \
        {input.trimmed} \
        -t {threads} \
        -l {params.seed_length} \
        -n {params.missing_prob} \
        -k {params.max_seed_diff} \
        -o {params.max_open_gap}  > {output.sai} 2> {log}
    """

rule bwa_samse:
    """
    Perform single-end read mapping on a 'bwa aln' suffix array.
    extension: '.collapsed' || '.collapsed.truncated' || '.singleton.truncated'

    # Benchmarks max RSS:
    | depth | .collapsed | .collapsed.truncated | .singleton.truncated |
    | ----- | ---------- | -------------------- | -------------------- |
    | 0.01X | 4588       | 7.43                 | 7.35                 |
    | 0.05X | 4989       | 7.44                 | 18.23                |
    | 0.10X | 4583       | 21.67                | 29.45                |

    # Benchmarks max h:m:s:
    | depth | .collapsed | .collapsed.truncated | .singleton.truncated |
    | ----- | ---------- | -------------------- | -------------------- |
    | 0.01X |            |                      |                      |
    | 0.05X | 0:25:30    | 0:02:53              | 0:02:48              |
    """
    input:
        trimmed   = "results/02-preprocess/01-adapter_removal/{sample}/{sample}.{extension}.gz",
        sai       = rules.bwa_aln.output.sai,
        reference = ReferenceGenome.get_path()
    output:
        sam       = temp(pipe("results/02-preprocess/02-align/{sample}/{sample}.bwaaln.{extension}.sam"))
    params:
        RG        = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina' ## ADD IT LATER (-r argument)
    resources:
        runtime   = 60,
        mem_mb    = 5000,
        cores     = lambda w, threads: threads
    log:       "logs/02-preprocess/02-align/bwa_samse/{sample}.{extension}.log"
    benchmark: "benchmarks/02-preprocess/02-align/bwa_samse/{extension}/{sample}-bench.tsv"
    conda:     "../envs/bwa-0.7.17.yml"
    #group:     "bwaaln"
    priority:  50
    threads:   1
    shell: """
        bwa samse -r '{params.RG}' {input.reference} {input.sai} {input.trimmed} > {output.sam} 2> {log}
    """

rule bwa_sampe:
    """
    Perform paired-end read mapping on a 'bwa aln' suffix array.
    extension: '.pair1.truncated' || '.pair2.truncated'

    # Benchmarks :
    | depth | max h:m:s | max_rss | 
    | ----- | --------- | ------- | 
    | 0.01X |           | 4478    |
    | 0.05X |  0:29:33  | 4618    |
    | 0.10X |           | 4755    |
    """
    input:
        pair1     = "results/02-preprocess/01-adapter_removal/{sample}/{sample}.pair1.truncated.gz",
        pair2     = "results/02-preprocess/01-adapter_removal/{sample}/{sample}.pair2.truncated.gz",
        sai1      = "results/02-preprocess/02-align/{sample}/{sample}.bwaaln.pair1.truncated.sai",
        sai2      = "results/02-preprocess/02-align/{sample}/{sample}.bwaaln.pair2.truncated.sai",
        reference = ReferenceGenome.get_path()
    output:
        sam       = temp(pipe("results/02-preprocess/02-align/{sample}/{sample}.bwaaln.paired.sam"))
    params:
        RG        = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina' # @TODO ADD it LATER
    resources:
        runtime   = 60,
        mem_mb    = 5000,
        cores     = lambda w, threads: threads
    log:       "logs/02-preprocess/02-align/bwa_sampe/{sample}.log"
    benchmark: "benchmarks/02-preprocess/02-align/bwa_sampe/{sample}.paired-bench.tsv"
    conda:     "../envs/bwa-0.7.17.yml"
    #group:     "bwaaln"
    priority: 50
    threads: 1
    shell: """
        bwa sampe -r '{params.RG}' {input.reference} {input.sai1} {input.sai2} {input.pair1} {input.pair2} > {output.sam} 2> {log}
    """


rule samtools_merge_aln:
    """
    Merge the outputs of bwa_samse & bwa_sampe for a given sample, and output a single merged bam file.

    # Benchmarks:
    | depth | max h:m:s | max_rss |
    | ----- | --------- | ------- | 
    | 0.01X |           |         |
    | 0.05X | 0:1:15    | 12.81   |
    | 0.10X |           | 30.12   |
    """
    input:
        paired_end = "results/02-preprocess/02-align/{sample}/{sample}.bwaaln.paired.tmp-bam",
        single_end = expand(
            "results/02-preprocess/02-align/{{sample}}/{{sample}}.bwaaln.{extension}.tmp-bam",
            extension=["collapsed", "collapsed.truncated", "singleton.truncated"]
        )
    output:
        merged    = "results/02-preprocess/02-align/{sample}/{sample}.bwaaln.merged.bam"
    params:
        RG        = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina' # @TODO ADD it LATER
    resources:
        runtime   = 60,
        mem_mb    = 128,
        cores     = lambda w, threads: threads
    log:       "logs/02-preprocess/02-align/samtools_merge/{sample}.log"
    benchmark: "benchmarks/02-preprocess/02-align/samtools_merge/{sample}-bench.tsv"
    conda:     "../envs/samtools-1.15.yml"
    #group:     "bwaaln"
    threads:   1
    priority:  60
    shell: """
        samtools merge -@ {threads} -o - {input.paired_end} {input.single_end} \
        | samtools addreplacerg -r '{params.RG}' -w -OBAM -o {output.merged} - 2> {log}
    """


# ----------------------------------------------------------------------------------------------- #
# ---- 01-B. Align to genome using BWA mem

"""
    Use bwa mem to perform read alignment and generate suffix array. This is much much faster than
    the traditional bwa aln + offers multithreading capability. but gives off less reliable results.
    Useful for prototyping, debugging, or when exactness isn't paramount.
    
    See:  Xu, W, Lin, Y, Zhao, K, et al. An efficient pipeline for ancient DNA mapping and recovery
          of endogenous ancient DNA from whole-genome sequencing data.
          https://doi.org/10.1002/ece3.7056 

    See:  Oliva, A., Tobler, R., Llamas, B. and Souilmi, Y. (2021), Additional evaluations show that specific
          BWA-aln settings still outperform BWA-mem for ancient DNA data alignment. 
          https://doi.org/10.1002/ece3.8297
"""

rule bwa_mem_se:
    """
    Perform a one pass single-end read mapping using bwa mem
    extension: '.collapsed' || '.collapsed.truncated' || '.singleton.truncated'
    """
    input:
        trimmed   = "results/02-preprocess/01-adapter_removal/{sample}/{sample}.{extension}.gz",
        reference = ReferenceGenome.get_path(),
        bwt       = rules.index_reference_genome.output.bwt
    output:
        sam       = temp(pipe("results/02-preprocess/02-align/{sample}/{sample}.bwamem.{extension}.sam"))
    params:
        RG        = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina'
    resources:
        cores     = lambda w, threads: threads
    log:       "logs/02-preprocess/02-align/bwa_mem_se/{sample}.bwamem.{extension}.log"
    benchmark: "benchmarks/02-preprocess/02-align/bwa_mem_se/{extension}/{sample}.bwamem-bench.tsv"
    conda:     "../envs/bwa-0.7.17.yml"
    group:     "bwamem"
    priority:  50
    threads:   4
    shell: """
        bwa mem -M -t {threads} -R \'{params.RG}\' -p {input.reference} {input.trimmed} > {output.sam} 2> {log}
    """


rule bwa_mem_pe:
    """
    Perform a one pass paired-end read mapping using bwa mem
    extension: '.pair1.truncated' || '.pair2.truncated'
    """
    input:
        pair1     = "results/02-preprocess/01-adapter_removal/{sample}/{sample}.pair1.truncated.gz",
        pair2     = "results/02-preprocess/01-adapter_removal/{sample}/{sample}.pair2.truncated.gz",
        reference = ReferenceGenome.get_path(),
        bwt       = rules.index_reference_genome.output.bwt
    output:
        sam       = temp(pipe("results/02-preprocess/02-align/{sample}/{sample}.bwamem.paired.sam"))
    params:
        RG        = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina'
    resources:
        cores     = lambda w, threads: threads
    log:       "logs/02-preprocess/02-align/bwa_mem_pe/{sample}.bwamem.paired.log"
    benchmark: "benchmarks/02-preprocess/02-align/bwa_mem_pe/{sample}.paired-bench.tsv"
    conda:     "../envs/bwa-0.7.17.yml"
    group:     "bwamem"
    priority:  50
    threads:   4
    shell: """
        bwa mem -M -t {threads} -R \'{params.RG}\' {input.reference} {input.pair1} {input.pair2} > {output.sam} 2> {log}
    """


rule samtools_merge_mem:
    """
    Merge the outputs of bwa_mem_se & bwa_mem_pe for a given sample, and output a single merged bam file.
    """
    input:
        paired_end = "results/02-preprocess/02-align/{sample}/{sample}.bwamem.paired.tmp-bam",
        single_end = expand(
            "results/02-preprocess/02-align/{{sample}}/{{sample}}.bwamem.{extension}.tmp-bam",
            extension=["collapsed", "collapsed.truncated", "singleton.truncated"]
        )
    output:
        merged     = "results/02-preprocess/02-align/{sample}/{sample}.bwamem.merged.bam"
    params:
        RG        = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina'
    resources:
        cores     = lambda w, threads: threads
    log:       "logs/02-preprocess/02-align/samtools_merge/{sample}.log"
    benchmark: "benchmarks/02-preprocess/02-align/samtools_merge/{sample}.bwamem-bench.tsv"
    conda:     "../envs/samtools-1.15.yml"
    group:     "bwamem"
    priority:  60
    threads:   4
    shell: """
        samtools merge -@ {threads} -o - {input.paired_end} {input.single_end} \
        | samtools addreplacerg -r '{params.RG}' -w -OBAM -o {output.merged} - 2> {log}
    """


# ------------------------------------------------------------------------------------------------------------------- #
# ---- 00. Resolve the required aligner algorithm.

def assign_aligner_algorithm(wildcards):
    """
    Decide on the appropriate bwa algorithm (aln or mem), based on the user input.
    """

    if config["preprocess"]["bwa"]["aligner"] == "mem":
        if config["preprocess"]["bwa"]["collapsed-only"]:
            return expand(rules.bwa_mem_se.output.sam, sample="{sample}", extension="collapsed")
        else:
            return rules.samtools_merge_mem.output.merged
    elif config["preprocess"]["bwa"]["aligner"] == "aln":
        if config["preprocess"]["bwa"]["collapsed-only"]:
            return expand(rules.bwa_samse.output.sam, sample="{sample}", extension="collapsed")
        else:
            return rules.samtools_merge_aln.output.merged


ruleorder: multiqc_align_fastqs > multiqc
use rule multiqc as multiqc_align_fastqs with:
    input:
        required_files = lambda w: temp(expand(expand(
                rules.run_fastqc_bam.output.data,
                directory=dirname(assign_aligner_algorithm(w)),
                file=basename(splitext(assign_aligner_algorithm(w))[0])
            ),
            sample = get_all_samples_ids(w)
        )),
    output: 
        html = report("results/00-qc/03-align-fastqs/multiqc-report.html",
            caption     = "../report/03-align-fastqs/multiqc-bwa.rst",
            category    = "02. Align"
        )
    params:
        extra_args = rules.multiqc.params.extra_args,
        extra_dirs = expand("results/02-preprocess/{subdir}", subdir=["01-adapter_removal", "02-align"]),