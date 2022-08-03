configfile: "./config/config.yml"

output_dir = "results/01-preprocess/"

rule adapter_removal:
    input:
        forwd = rules.merge_chromosomes.output.forw,
        revrs = rules.merge_chromosomes.output.rev
    output:
        trimmed = temp("{output_dir}/00-adapter_removal/{sample}_trimmed.fq.gz")
    params:
        min_quality = 30
    conda: "envs/adapterremoval-2.3.3"
    priority: 4
    threads:  4
    shell: """
        AdapterRemoval --threads {threads}               \
                       --gzip                            \
                       --file1 {input.forwd}             \
                       --file2 {input.revrs}             \
                       --basename {wildcards.sample}.fq  \
                       --collapse                        \
                       --combined-output                 \
                       --output1 {output.trimmed}        \
                       --minquality {params.min_quality}
    """

rule bwa_mem:
    input:
        trimmed = rules.adapter_removal.output.trimmed,
        reference = rules.index_refgen.output.bwt
    output:
        sam = temp("{output_dir}/01-bwa-mem/{sample}.sam")
    params:
        RG='@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina\\t'
    conda: "envs/bwa-0.7.17.yml"
    priority: 5
    threads: 4
    shell: """
        bwa mem -t {threads}
                -R \'{params.RG}\'
                -p ${{input.reference}%.bwt}
                {input.trimmed} > {output.sam}"
    """