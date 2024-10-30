# ---- Set config variables
configfile: "./config/config.yml"

rule run_fastqc_bam:
    input: "{directory}/{file}.bam"
    output:
        data = temp("{directory}/{file}_fastqc.zip"),
        html = temp("{directory}/{file}_fastqc.html")
    params:
        extra   = "--svg",
        tmpdir = config['tempdir'],
    resources:
        mem_mb = lambda w, threads: max(512, threads * 250)
    log:       "logs/00-qc/{directory}/{file}/run_fastqc.log"
    benchmark: "benchmarks/00-qc/{directory}/{file}/run_fastqc.tsv"
    conda:     "../envs/fastqc-0.12.1.yml"
    threads:   1
    shell: """
        fastqc {input} --dir {params.tmpdir} --outdir {wildcards.directory} {params.extra} > {log} 2>&1
    """

use rule run_fastqc_bam as run_fastqc_fq with:
    input: "{directory}/{file}.fq.gz"

use rule run_fastqc_bam as run_fastqc_sam with:
    input: "{directory}/{file}.sam"

rule multiqc:
    output:
        html = "{file}.html",
        data = directory("{file}_data")
    params:
        #basedir  = lambda w, output: dirname(output.html),
        #basename = lambda w, output: basename(splitext(output.html)[0]),
        extra_args    = "",
        extra_dirs    = ""
    conda: "../envs/multiqc-1.21.yml"
    shell: """
        BASEDIR=$(dirname {output.html})
        BASENAME=$(basename {output.html})
        multiqc {input} {params.extra_dirs} --outdir $BASEDIR --filename ${{BASENAME%.html}} --force --no-ansi {params.extra_args}
    """
