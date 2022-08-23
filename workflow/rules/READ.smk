# ---- Utility functions
def READ_output(wildcards):
    """
    Returns the list of required generation_wise comparisons from READ
    """
    with checkpoints.get_samples.get().output[0].open() as f:
        samples = str.split(f.readline().replace('\n', ''), '\t')
        generations = set([sample.split('_')[0] for sample in samples])
        return expand("results/03-kinship/READ/{generation}/READ_results", generation=generations)

def READ_bam_input_id(wildcards):
    # Run through the initial samples files and extract pedigree ids 
    with checkpoints.get_samples.get().output[0].open() as f:
        samples = str.split(f.readline().replace('\n', ''), '\t')
        ids     = set([sample.split('_')[1] for sample in samples])
        return expand(
            "results/02-preprocess/06-mapdamage/{generation}_{ids}/{generation}_{ids}.srt.rmdup.rescaled.bam",
            ids=ids,
            generation="{generation}",
        )

def READ_bam_samples_id(wildcards):
    # Run through the initial samples files and extract pedigree ids 
    with checkpoints.get_samples.get().output[0].open() as f:
        samples = str.split(f.readline().replace('\n', ''), '\t')
        ids     = set([sample.split('_')[1] for sample in samples])
        return expand(
            "{generation}_{ids}",
            ids=ids,
            generation="{generation}",
        )

rule READ_bam_list:
    """
    Create an ordered list of bam files using the original samples file
    """
    input:
        bamlist = READ_bam_input_id
    output:
        bamlist = "results/03-kinship/READ/{generation}/{generation}.bam.list"
    priority: 15
    shell: """
        ls {input.bamlist} > {output.bamlist}
    """

rule pileup_READ:
    """
    Perform batch mpileup
    """
    input:
        bams      = READ_bam_input_id,
        bamlist   = rules.READ_bam_list.output.bamlist,
        targets   = config["kinship"]["targets"],
        reference = config["refgen"],
        fai       = config["refgen"] + ".fai",
    output:
        pileup    = "results/03-kinship/READ/{generation}/{generation}.pileup"
    params:
        min_MQ = 30,
        min_BQ = 30
    priority: 15
    conda: "../envs/samtools-1.15.yml"
    shell: """
        samtools mpileup -B -q {params.min_MQ} -Q {params.min_BQ} -l {input.targets} -f {input.reference} -b {input.bamlist} > {output}
    """

rule ANGSD_random_haploid:
    """
    Perform random pseudo-haploid variant calling with ANGSD
    """
    input:
        pileup    = rules.pileup_READ.output.pileup,
        bamlist   = rules.READ_bam_list.output.bamlist,
        reference = config['refgen'],
        fai       = config["refgen"] + ".fai",
    output:
        haplos    = "results/03-kinship/READ/{generation}/{generation}.haplo.gz"
    params:
        out       = "results/03-kinship/READ/{generation}/{generation}"
    threads: 4
    priority: 15
    conda: "../envs/angsd-0.939.yml"
    shell: """
        angsd -pileup {input.pileup}         \
        -fai {input.fai}                     \
        -nInd $(cat {input.bamlist} | wc -l) \
        -rf <(seq 1 22)                      \
        -dohaplocall 1                       \
        -doCounts 1                          \
        -nthreads {threads}                  \
        -out {params.out}
    """

rule ANGSD_haplo_to_plink:
    """
    Convert a .haplo.gz ANGSD file to a set of PLINK .tped / .tfam files 
    """
    input:
        haplos  = rules.ANGSD_random_haploid.output.haplos,
        bamlist = rules.READ_bam_list.output.bamlist 
    output:
        tped = "results/03-kinship/READ/{generation}/{generation}.tped",
        tfam = "results/03-kinship/READ/{generation}/{generation}.tfam",
    params:
        outputname = "results/03-kinship/READ/{generation}/{generation}"
    conda: "../envs/angsd-0.939.yml"
    priority: 15
    shell: """
        haploToPlink {input.haplos} {params.outputname}
        sed -i 's/N/0/g' {output.tped}
        cat {input.bamlist} | xargs basename -a \
                            | grep -oP '^[^.]+(?=(\.[^.]+)*(\.bam$))' \
                            | awk 'BEGIN{{OFS="\t"}}{{print "{wildcards.generation}", $1, 0, 0, 0, 0}}' \
                            > {output.tfam}
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- Serial READ
rule run_READ:
    """
    Run the original, serialized version of READ
    See: Monroy Kuhn JM, Jakobsson M, Günther T (2018) Estimating genetic kin relationships in prehistoric populations.
         PLoS ONE 13(4): e0195491. https://doi.org/10.1371/journal.pone.0195491
         https://bitbucket.org/tguenther/read.git
    """
    input:
        tped = rules.ANGSD_haplo_to_plink.output.tped,
        tfam = rules.ANGSD_haplo_to_plink.output.tfam
    output:
        results = "results/03-kinship/READ/{generation}/READ_results",
        plot    = "results/03-kinship/READ/{generation}/READ_results_plot.pdf"
    params:
        window_size = config["kinship"]["READ"]["window-size"]
    conda: "../envs/READ.yml"
    shell: """
        cwd=$(pwd)
        cd $(dirname {output.results})
        ln -srf $(which READscript.R) READscript.R
        touch meansP0_AncientDNA_normalized
        python2 $(which READ.py)  {wildcards.generation} median - --window_size {params.window_size}
        cd -
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- Parallel READ
rule READ_split_plink:
    """
    Split a genome wide set of .tped/.tfam files into a set of per-chromosome .tped/.tfam files.
    """
    input:
        tped = rules.ANGSD_haplo_to_plink.output.tped,
        tfam = rules.ANGSD_haplo_to_plink.output.tfam
    output:
        tped = "results/03-kinship/READ/{generation}/{chr}/{generation}_{chr}.tped",
        tfam = "results/03-kinship/READ/{generation}/{chr}/{generation}_{chr}.tfam"
    params:
        out  = "results/03-kinship/READ/{generation}/{chr}/{generation}_{chr}"
    conda: "../envs/plink-1.9.yml"
    group: "READ"
    resources: scatter = 1
    priority: 15
    shell: """
        plink --tped {input.tped}   \
              --tfam {input.tfam}   \
              --chr {wildcards.chr} \
              --out {params.out}    \
              --recode transpose    \
    """

rule run_READ_parallel:
    """
    Run a modified, parallelized version of READ
    See: Monroy Kuhn JM, Jakobsson M, Günther T (2018) Estimating genetic kin relationships in prehistoric populations.
         PLoS ONE 13(4): e0195491. https://doi.org/10.1371/journal.pone.0195491
         https://bitbucket.org/tguenther/read.git
    """
    input:
        tped = "results/03-kinship/READ/{generation}/{chr}/{generation}_{chr}.tped",
        tfam = "results/03-kinship/READ/{generation}/{chr}/{generation}_{chr}.tfam"
    output:
        results = "results/03-kinship/READ/{generation}/{chr}/Read_intermediate_output"
    priority: 15
    group: "READ"
    resources: scatter = 1
    shell: """
        cwd=$(pwd)
        cd $(dirname {output.results})
        python2 ${{cwd}}/methods/read/READ_parallel.py $(ls *.tped | sed 's/.tped//g') median -
        cd -
    """

rule merge_READ_parallel:
    """
    Merge the intermediary results of parallel-READ into a single output.
    """
    input:
        results=expand(rules.run_READ_parallel.output.results, generation="{generation}", chr=range(1,23))
    output:
        "results/03-kinship/READ/{generation}/Read_intermediate_output"
    priority: 15
    group: "READ"
    resources: scatter = 1
    shell: """
        cat {input.results} > {output}
    """

rule plot_READ_parallel:
    """
    Compute kinship and plot the output of parallel-read using the merged intermediary results.
    """
    input:
        results = rules.merge_READ_parallel.output,
        tped = rules.ANGSD_haplo_to_plink.output.tped,
        tfam = rules.ANGSD_haplo_to_plink.output.tfam
    output:
        results = "results/03-kinship/READ/{generation}/READ_results_parallel"
    priority: 15
    shell: """
        cwd=$(pwd)
        cd $(dirname {output.results})
        ln -srf ${{cwd}}/methods/read/READscript.R READscript.R
        touch meansP0_AncientDNA_normalized
        python2 ${{cwd}}/methods/read/READ_merge.py  {wildcards.generation} median -
        cd -
    """