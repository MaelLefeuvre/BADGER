
# ---- Utility functions
def READ_output(wildcards):
    """
    Returns the list of required generation_wise comparisons from READ
    @ TODO: This checkpoint is not needed anymore
    """
    with checkpoints.get_samples.get().output[0].open() as f:
        samples     = str.split(f.readline().replace('\n', ''), '\t')
        generations = set([sample.split('_')[0] for sample in samples])
        return expand("results/04-kinship/READ/{generation}/READ_results", generation=generations)
    #generations = range(1, config['ped-sim']['replicates'])
    #return expand("results/04-kinship/READ/{generation}/READ_results", generation=generations)


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
        means   = "results/04-kinship/READ/{generation}/meansP0_AncientDNA_normalized",
        raw     = "results/04-kinship/READ/{generation}/READ_output_ordered",
        plot    = "results/04-kinship/READ/{generation}/READ_results_plot.pdf",
        tempraw = temp("results/04-kinship/READ/{generation}/Read_intermediate_output"),
    params:
        window_size = config["kinship"]["READ"]["window-size"],
        norm_method = config["kinship"]["READ"]["norm-method"],
        norm_value  = get_READ_norm_value,
        basename    = lambda wildcards, input: splitext(input.tplink[0])[0],
    conda: "../envs/READ.yml"
    log: "logs/04-kinship/READ/run_READ/{generation}.log"
    shell: """
        cwd=$(pwd)
        cd $(dirname {output.results})               2>  $cwd/{log}
        ln -srf $(which READscript.R) READscript.R   2>> $cwd/{log}
        touch meansP0_AncientDNA_normalized          2>> $cwd/{log}
        python2 $(which READ.py) $cwd/{params.basename} {params.norm_method} {params.norm_value} \
        --window_size {params.window_size} >> $cwd/{log} 2>&1

        # Remove symlink
        find . -type l -exec rm {{}} \;
    """
