# ---- Utility functions
def READ_output(wildcards):
    """
    Returns the list of required generation_wise comparisons from READ
    """
    # Get the number of expected generations.
    template = "results/04-kinship/READ/{generation}/READ_results"
    return expand(template, generation = get_generations())


def READ_define_random_haploid_caller(wildcards):
    """
    Define the input for READ, based on which variant-calling method was 
    requested by the user.
    """
    match config['variant-calling']['caller']:
        case "ANGSD":
            basedir = "results/03-variant-calling/02-ANGSD/"
        case "pileupCaller":
            basedir = "results/03-variant-calling/02-pileupCaller/"
        case _:
            raise RuntimeError("Incorrect pileupCaller mode selected.")    

    return multiext(basedir + "{generation}/{generation}", ".tped", ".tfam")


def get_READ_norm_value(wildcards):
    """
    Inject the user-defined normalization value within READ's command line 
    arguments (if the user requested it). Or add a placeholder "-" character
    in place if the user wishes to compute that value from the cohort.
    """
    if config["kinship"]["READ"]["norm-method"] == "value":
        return config["kinship"]["READ"]["norm-value"]
    else:
        return "-"

# ------------------------------------------------------------------------------------------------ #

rule run_READ:
    """
    Run the original, serialized version of READ
    See: Monroy Kuhn JM, Jakobsson M, Günther T (2018) Estimating genetic kin relationships in prehistoric populations.
         PLoS ONE 13(4): e0195491. https://doi.org/10.1371/journal.pone.0195491
         
    Repo: https://bitbucket.org/tguenther/read.git

    # Benchmarks: 
    | depth | max h:m:s | max_rss |
    | ----- | --------- | ------- |
    | 0.01X |           |         |
    | 0.05X | 0:00:49   | 47.87   |
    | 0.10X | 0:04:10   | 30.97   |
    """
    input:
        tplink = READ_define_random_haploid_caller
    output:
        results = "results/04-kinship/READ/{generation}/READ_results",
        means   = "results/04-kinship/READ/{generation}/meansP0_AncientDNA_normalized",
        raw     = "results/04-kinship/READ/{generation}/READ_output_ordered",
        plot    = report("results/04-kinship/READ/{generation}/READ_results_plot.pdf",
            caption     = "../report/06-READ/run-READ/plot.rst",
            category    = "04. Kinship",
            subcategory = "READ",
            labels      = {"replicate": "{generation}", "figure": "READ results plot"}
        ),
        tempraw = temp("results/04-kinship/READ/{generation}/Read_intermediate_output"),
    params:
        window_size = config["kinship"]["READ"]["window-size"],
        norm_method = config["kinship"]["READ"]["norm-method"],
        norm_value  = get_READ_norm_value,
        basename    = lambda wildcards, input: splitext(input.tplink[0])[0],
    resources:
        runtime = 60,
        mem_mb  = 128,
        cores   = lambda w, threads: threads
    log:       "logs/04-kinship/READ/run_READ/{generation}.log"
    benchmark: "benchmarks/04-kinship/READ/run_READ/{generation}.tsv"
    conda:     "../envs/READ.yml"
    threads:   1
    shell: r""" # use of raw string to prevent python escape character SyntaxWarning
        cwd=$(pwd)
        cd $(dirname {output.results})               2>  $cwd/{log}
        ln -srf $(which READscript.R) READscript.R   2>> $cwd/{log}
        touch meansP0_AncientDNA_normalized          2>> $cwd/{log}
        python2 $(which READ.py) $cwd/{params.basename} {params.norm_method} {params.norm_value} \
        --window_size {params.window_size} >> $cwd/{log} 2>&1

        # Remove symlink
        find . -type l -exec rm {{}} \;
    """
