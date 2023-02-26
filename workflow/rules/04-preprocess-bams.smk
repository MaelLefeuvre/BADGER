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
    """
    
    input:
        bam       = rules.samtools_filter_unmapped.output.bam,
        reference = config["reference"]
    output:
        bam       = "results/02-preprocess/04-sort/{sample}/{sample}.srt.bam"
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
    """
    input:
        tmpdir  = config["tempdir"]
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
    """
    input:
        bam       = define_dedup_input_bam,
        bai       = lambda wildcards: define_dedup_input_bam(wildcards) + ".bai",
        reference = config["reference"],
        metadata     = "results/meta/pipeline-metadata.yml"
    output:
        bam = "results/02-preprocess/06-mapdamage/{sample}/{sample}.srt.rmdup.rescaled.bam"
    params:
        downsample_seed = define_mapdamage_seed
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
        -i {input.bam}                   \
        -r {input.reference}             \
        --rescale                        \
        --folder $(dirname {output.bam}) \
        --rescale-out {output.bam}       \
        {params.downsample_seed}         \
        --verbose > {log} 2>&1 
    """


