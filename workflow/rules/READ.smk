# ---- Utility functions
def READ_output(wildcards):
    """
    Returns the list of required generation_wise comparisons from READ
    """
    with checkpoints.get_samples.get().output[0].open() as f:
        samples = str.split(f.readline().replace('\n', ''), '\t')
        generations = set([sample.split('_')[0] for sample in samples])
        return expand("results/04-kinship/READ/{generation}/READ_results", generation=generations)


def READ_define_random_haploid_caller(wildcards):
    match config['variant-calling']['caller']:
        case "ANGSD":
            basedir = "results/03-variant-calling/02-ANGSD/"
        case "pileupCaller":
            basedir = "results/03-variant-calling/02-pileupCaller/"
        case _:
            raise RuntimeError("Incorrect pileupCaller mode selected.")    

    return multiext(basedir + "{generation}/{generation}", ".tped", ".tfam")


def get_READ_norm_value(wildcards):
    if config["kinship"]["READ"]["norm-method"] == "value":
        return config["kinship"]["READ"]["norm-value"]
    else:
        return "-"

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
        tplink = READ_define_random_haploid_caller
    output:
        results = "results/04-kinship/READ/{generation}/READ_results",
        plot    = "results/04-kinship/READ/{generation}/READ_results_plot.pdf"
    params:
        window_size = config["kinship"]["READ"]["window-size"],
        norm_method = config["kinship"]["READ"]["norm-method"],
        norm_value  = get_READ_norm_value,
        basename    = lambda wildcards, input: splitext(input.tplink[0])[0],
    conda: "../envs/READ.yml"
    log: "logs/04-kinship/run_READ/{generation}.log"
    shell: """
        cwd=$(pwd)
        cd $(dirname {output.results})
        ln -srf $(which READscript.R) READscript.R
        touch meansP0_AncientDNA_normalized
        python2 $(which READ.py) $cwd/{params.basename} {params.norm_method} {params.norm_value} --window_size {params.window_size} > $cwd/{log} 2>&1 
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
        tped = "results/04-kinship/READ/{generation}/{chr}/{generation}_{chr}.tped",
        tfam = "results/04-kinship/READ/{generation}/{chr}/{generation}_{chr}.tfam"
    params:
        out  = "results/04-kinship/READ/{generation}/{chr}/{generation}_{chr}"
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
        tped = "results/04-kinship/READ/{generation}/{chr}/{generation}_{chr}.tped",
        tfam = "results/04-kinship/READ/{generation}/{chr}/{generation}_{chr}.tfam"
    output:
        results = "results/04-kinship/READ/{generation}/{chr}/Read_intermediate_output"
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
        "results/04-kinship/READ/{generation}/Read_intermediate_output"
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
        results = "results/04-kinship/READ/{generation}/READ_results_parallel"
    priority: 15
    shell: """
        cwd=$(pwd)
        cd $(dirname {output.results})
        ln -srf ${{cwd}}/methods/read/READscript.R READscript.R
        touch meansP0_AncientDNA_normalized
        python2 ${{cwd}}/methods/read/READ_merge.py  {wildcards.generation} median -
        cd -
    """