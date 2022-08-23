configfile: "./config/config.yml"

rule adapter_removal_pe:
    input:
        forwd = "{preprocess_dir}/00-raw/{sample}_s1.fq.gz",
        revrs = "{preprocess_dir}/00-raw/{sample}_s2.fq.gz"
    output:
        trimmed   = "{preprocess_dir}/01-adapter_removal/{sample}/{sample}.collapsed.gz",
        truncated = "{preprocess_dir}/01-adapter_removal/{sample}/{sample}.collapsed.truncated.gz",
        discarded = "{preprocess_dir}/01-adapter_removal/{sample}/{sample}.discarded.gz",
        pair1     = "{preprocess_dir}/01-adapter_removal/{sample}/{sample}.pair1.truncated.gz",
        pair2     = "{preprocess_dir}/01-adapter_removal/{sample}/{sample}.pair2.truncated.gz",
        singleton = "{preprocess_dir}/01-adapter_removal/{sample}/{sample}.singleton.truncated.gz"
    params:
        base_name   = "{preprocess_dir}/01-adapter_removal/{sample}/{sample}",
        min_overlap = config["preprocess"]["trimming"]["min-overlap"],
        min_length  = config["preprocess"]["trimming"]["min-overlap"]
    log: "logs/{preprocess_dir}/adapter_removal_pe/{sample}.log"

    conda: "../envs/adapterremoval-2.3.3.yml"
    priority: 4
    threads:  4
    shell: """
        AdapterRemoval --threads {threads}                       \
                       --gzip                                    \
                       --file1 {input.forwd}                     \
                       --file2 {input.revrs}                     \
                       --basename {params.base_name}             \
                       --collapse                                \
                       --minlength {params.min_length}           \
                       --minadapteroverlap {params.min_overlap} 2> {log}
    """

rule bwa_mem:
    input:
        trimmed   = rules.adapter_removal_pe.output.trimmed,
        reference = config["refgen"],
        bwt       = rules.index_reference_genome.output.bwt
    output:
        sam = "{preprocess_dir}/02-align/{sample}.bwamem.sam"
    params:
        RG = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina'
    log: "logs/{preprocess_dir}/bwa_mem/{sample}.log"
    conda: "../envs/bwa-0.7.17.yml"
    priority: 5
    threads: 4
    shell: """
        bwa mem -t {threads} -R \'{params.RG}\' -p {input.reference} {input.trimmed} > {output.sam} 2> {log}
    """

rule bwa_aln:
    """
    Use bwa aln to perform alignment. Gives off very robust results, at the cost of an impoverished runtime performance
    compared to bwa mem.
    See:    Oliva, A., Tobler, R., Llamas, B. and Souilmi, Y. (2021), Additional evaluations show that specific
            BWA-aln settings still outperform BWA-mem for ancient DNA data alignment. Ecol Evol, 11: 18743-18748. 
            https://doi.org/10.1002/ece3.8297

    extension: '.collapsed' || '.collapsed.truncated' || '.pair1.truncated' || '.pair1.truncated' || '.singleton.truncated'
    """
    input:
        trimmed = "{preprocess_dir}/01-adapter_removal/{sample}/{sample}.{extension}.gz",
        reference = config["refgen"],
        bwt       = rules.index_reference_genome.output.bwt
    output:
        sai     = "{preprocess_dir}/02-align/{sample}/{sample}.{extension}.sai"
    params:
        seed_length = 1024,
        missing_prob = 0.01,
        max_open_gap = 2
    log: "logs/{preprocess_dir}/bwa_aln/{sample}.{extension}.log"
    conda: "../envs/bwa-0.7.17.yml"
    threads: 16
    shell: """
        bwa aln {input.reference}        \
                {input.trimmed}          \
                -t {threads}             \
                -l {params.seed_length}  \
                -n {params.missing_prob} \
                -o {params.max_open_gap} > {output.sai} 2> {log}
    """

rule bwa_samse:
    input:
        trimmed   = "{preprocess_dir}/01-adapter_removal/{sample}/{sample}.{extension}.gz",
        sai       = rules.bwa_aln.output.sai,
        reference = config["refgen"]
    output:
        sam       = "{preprocess_dir}/02-align/{sample}/{sample}.{extension}.sam"
    #params:
    #    RG        = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina' ## ADD IT LATER (-r argument)
    log: "logs/{preprocess_dir}/bwa_samse/{sample}.{extension}.log"
    conda: "../envs/bwa-0.7.17.yml"
    threads: 1
    shell: """
        bwa samse {input.reference} {input.sai} {input.trimmed} > {output.sam} 2> {log}
    """

rule bwa_sampe:
    input:
        pair1 = "{preprocess_dir}/01-adapter_removal/{sample}/{sample}.pair1.truncated.gz",
        pair2 = "{preprocess_dir}/01-adapter_removal/{sample}/{sample}.pair2.truncated.gz",
        sai1  = "{preprocess_dir}/02-align/{sample}/{sample}.pair1.truncated.sai",
        sai2  = "{preprocess_dir}/02-align/{sample}/{sample}.pair2.truncated.sai",
        reference = config["refgen"]
    output:
        sam       = temp("{preprocess_dir}/02-align/{sample}/{sample}.paired.sam")
    #params:
    #    RG        = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina' # ADD it LATER
    log: "logs/{preprocess_dir}/bwa_sampe/{sample}.log"
    conda: "../envs/bwa-0.7.17.yml"
    threads: 1
    shell: """
        bwa sampe {input.reference} {input.sai1} {input.sai2} {input.pair1} {input.pair2} > {output.sam} 2> {log}
    """

rule samtools_merge:
    input:
        paired_end = "{preprocess_dir}/02-align/{sample}/{sample}.paired.sam",
        single_end = expand("{{preprocess_dir}}/02-align/{{sample}}/{{sample}}.{extension}.sam", extension=["collapsed", "collapsed.truncated", "singleton.truncated"],)
    output:
        merged = "{preprocess_dir}/02-align/{sample}/{sample}.merged.sam"
    log: "logs/{preprocess_dir}/samtools_merge/{sample}.log"
    threads: 1
    shell: """
        samtools merge -@ {threads} {output.merged} {input.paired_end} {input.single_end} 2> {log}
    """

def assign_aligner_algorithm(wildcards):
    if config["preprocess"]["bwa"]["aligner"] == "mem":
        return rules.bwa_mem.output.sam
    elif config["preprocess"]["bwa"]["aligner"] == "aln":
        return rules.samtools_merge.output.merged

rule samtools_filter_unmapped:
    input:
        sam = assign_aligner_algorithm,
        reference = config["refgen"],
        bwt       = rules.index_reference_genome.output.bwt
    output:
        bam = "{preprocess_dir}/03-filter/{sample}.bam"
    conda: "../envs/samtools-1.15.yml"
    params:
        min_MQ = config["preprocess"]["filter"]["min-MQ"]
    log: "logs/{preprocess_dir}/samtools_filter_unmapped/{sample}.log"
    priority: 6
    threads: 4
    shell: """
        samtools view --threads {threads}           \
                      --reference {input.reference} \
                      -q {params.min_MQ}            \
                      -F4 -Sb                       \
                      {input.sam} > {output.bam} 2> {log}
    """

rule samtools_sort:
    input:
        bam       = rules.samtools_filter_unmapped.output.bam,
        reference = config["refgen"]
    output:
        bam = "{preprocess_dir}/04-sort/{sample}.srt.bam"
    log: "logs/{preprocess_dir}/samtools_sort/{sample}.log"
    conda: "../envs/samtools-1.15.yml"
    priority: 7
    threads: 4
    shell: """
        samtools sort -@ {threads} --reference {input.reference} --output-fmt BAM -o {output.bam} {input.bam} 2> {log}
    """

rule samtools_rmdup:
    input:
        bam       = rules.samtools_sort.output.bam,
        reference = config["refgen"],
    output:
        bam = "{preprocess_dir}/05-dedup/{sample}.srt.rmdup.bam"
    params:
        max_length = config["preprocess"]["dedup"]["max-length"]
    log: "logs/{preprocess_dir}/samtools_rmdup/{sample}.log"
    conda: "../envs/samtools-1.15.yml"
    threads: 4
    priority: 8
    shell: """
        samtools markdup -rs --threads {threads}                   \
                         --reference {input.reference}             \
                         -l {params.max_length}                    \
                         --output-fmt BAM {input.bam} {output.bam} 2> {log}
    """
    
rule samtools_index:
    input:
        bam = "{bam}"
    output:
        bai = "{bam}.bai"
    threads: 4
    conda: "../envs/samtools-1.15.yml"
    shell: """
        samtools index -@ {threads} {input.bam}
    """

rule map_damage:
    input:
        bam       = rules.samtools_rmdup.output.bam,
        bai       = expand(rules.samtools_index.output.bai, bam=rules.samtools_rmdup.output.bam),
        reference = config["refgen"]
    output:
        bam = "{preprocess_dir}/06-mapdamage/{sample}/{sample}.srt.rmdup.rescaled.bam"
    log: "logs/{preprocess_dir}/samtools_mapdamage/{sample}.log"
    conda: "../envs/mapdamage-2.2.1.yml"
    threads: 1
    priority: 10
    resources: scatter=16
    shell: """
        mapDamage -i {input.bam}                   \
                  -r {input.reference}             \
                  --rescale                        \
                  --folder $(dirname {output.bam}) \
                  --rescale-out {output.bam}       \
                  --verbose > {log} 2>&1 
    """


