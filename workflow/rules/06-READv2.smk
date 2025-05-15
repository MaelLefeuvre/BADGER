# ---- Utility functions
def READv2_output(wildcards):
    """
    Returns the list of required generation_wise comparisons from READv2
    """
    # Get the number of expected generations.
    template = "results/04-kinship/READv2/{generation}/Read_Results.tsv" 
    return expand(template, generation = get_generations())


def READv2_define_random_haploid_caller(wildcards):
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

    return multiext(basedir + "{generation}/{generation}", ".bed", ".bim", ".fam")


def parse_READv2_norm_method(wildcards):
    """
    Inject the user-defined normalization value within READ's command line argument (if the user requested it).
    """
    norm_method = config['kinship']['READv2']['norm-method']
    norm_value  = config['kinship']['READv2']['norm-value']
    match norm_method:
        case None:
            return ""
        case "median" | "mean" | "max" :
            return f"--norm_method {norm_method} "
        case "value":
            if norm_value is None:
                raise RuntimeError("[READv2] Requested a set value as a normalisation value, but none was provided.")
            return f"--norm_method {norm_method} --norm_value {config['kinship']['READ']['norm-value']}"
        case other:
            raise RuntimeError(f"[READv2] Invalid norm-method: {other}")

def parse_READv2_window_estimate(wildcards):
    """
    Inject the user-defined window estimate strategy within READ's command line argument (if the user requested it)
    """
    window_estimate_requested = config['kinship']['READv2']['window-est']
    window_size               = config['kinship']['READv2']['window-size']
    output_optargs = ""
    if window_estimate_requested:
        output_optargs += "--window_est"
        if window_size is not None:
            try:
                output_optargs += f" --window_size {int(window_size)}"
            except:
                raise RuntimeError(f"[READv2] Invalid window-size provided: {window_size}")

    return output_optargs

def parse_READv2_alternate_thresholds(wildcards):
    """
    Inject the user-requested window threshold stragtegy within READv2's command line argument

    - When False: [0.96875,0.90625,0.8125,0.625]
    - When True:  [1-1/(2**4.5),1-1/(2**3.5),1-1/(2**2.5),1-1/(2**1.5)]
    """
    alternate_threshold_requested = config['kinship']['READv2']['2pow']
    return "--2pow" if alternate_threshold_requested else ""


def count_comparisons(wildcards, include_self = True):
    with open(rules.generate_bam_list.output.bamlist) as f:
        n = len([line for line in f.readlines()])
        combinations = (n * (n-1)) / 2
    return combinations + n if include_self else combinations

# ------------------------------------------------------------------------------------------------ #

rule run_READv2:
    """
    Run the updated version of READ:
    See   : https://doi.org/10.1101/2024.01.23.576660 
    Github: https://github.com/GuntherLab/READv2
    """
    input:
        bfile = READv2_define_random_haploid_caller
    output:
        results = "results/04-kinship/READv2/{generation}/Read_Results.tsv",
        meansP0 = "results/04-kinship/READv2/{generation}/meansP0_AncientDNA_normalized_READv2",
        plot    = report( # This will break in n > 45 samples, as READv2 only plots when under 1000 pairs.
            "results/04-kinship/READv2/{generation}/READ.pdf",
            caption     = "../report/06-READ/run-READ/plot.rst",
            category    = "04. Kinship",
            subcategory = "READv2",
            labels      = {"replicate": "{generation}", "figure": "READv2 results plot"}
        )
    params:
        basename              = lambda wildcards, input: splitext(input.bfile[0])[0],
        norm_method_optargs   = parse_READv2_norm_method,
        window_size_optargs   = parse_READv2_window_estimate,
        alt_threshold_optargs = parse_READv2_alternate_thresholds
    log:       "logs/04-kinship/READv2/run_READv2/{generation}.log"
    benchmark: "benchmarks/04-kinship/READv2/run_READv2/{generation}.tsv"
    conda:     "../envs/READv2.yml"
    threads:   1
    shell: """
        cwd=$(pwd)
        cd $(dirname {output.results}) > $cwd/{log} 2>&1
        READ2 \
        --input $cwd/{params.basename} \
        {params.norm_method_optargs} \
        {params.window_size_optargs} \
        {params.alt_threshold_optargs} >> $cwd/{log} 2>&1
    """
