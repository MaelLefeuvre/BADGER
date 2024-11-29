from os.path import splitext, basename
from functools import partial
localrules: merge_TKGWV2_results

# ---- Utility functions
def TKGWV2_output(wildcards):
    """
    Returns the list of required generation_wise comparisons from TKGWV2
    """
    # Get the number of expected generations.
    gen_no = config['ped-sim']['replicates']
    template = "results/04-kinship/TKGWV2/{generation}/{generation}-TKGWV2_Results.txt"
    return expand(template, generation=["ped" + str(i) for i in range(1, gen_no + 1)])


def get_TKGWV2_input_bams(wildcards, print_log = False, logfile = sys.stderr):
    apply_masking = config['preprocess']['pmd-rescaling']['apply-masking']
    rescaler = config['preprocess']['pmd-rescaling']['rescaler']
    if apply_masking:
        if print_log:
            print("Applying pmd-mask for TKGWV2", file=logfile)
        root = expand(rules.run_pmd_mask.output.bam, sample = "{generation}_{pairs}")


    elif rescaler is None:
        if print_log: 
            print("WARNING: Skipping PMD rescaling for TKGWV2!", file=logfile)
        root = expand(define_dedup_input_bam(wildcards), sample = "{generation}_{pairs}")
    else:
        if print_log:
            print(f"Using: {rescaler} bam files for TKGWV2", file=logfile)
        root = expand(define_rescale_input_bam(wildcards), sample = "{generation}_{pairs}")
    
    return expand(root,
        generation = "{generation}",
        pairs      = ["{pairA}", "{pairB}"],
        ext        = ["bam", "bam.bai"]
    )

# ------------------------------------------------------------------------------------------------ #

def TKGWV2_downsample_seed(wildcards):
    seed = config['kinship']['TKGWV2']['downsample-seed']
    if seed is None:
        with open(rules.meta.output.metadata) as f:
            metadata = yaml.load(f, Loader=yaml.loader.SafeLoader)
            seed     = metadata['seed']
    return seed

rule TKGWV2_downsample_bam:
    """
    Downsample the .bam file if there are more than 1_500_000 reads.
     - As is, this helper script is applied to all files ending with the ".bam" suffix within the working directory.
     - output files are identified with the "_subsampled.bam" suffix.
    """
    input:
        pairs =  get_TKGWV2_input_bams,
        metadata = "results/meta/pipeline-metadata.yml"
    output:
        pairA = "results/04-kinship/TKGWV2/{generation}/{pairA}_{pairB}/{generation}_{pairA}.srt.rmdup.rescaled_subsampled.bam",
        pairB = "results/04-kinship/TKGWV2/{generation}/{pairA}_{pairB}/{generation}_{pairB}.srt.rmdup.rescaled_subsampled.bam"
    params:
        workdir     = lambda wildcards, output: dirname(output.pairA),
        downsampleN = config['kinship']['TKGWV2']['downsample-N'],
        seed        = TKGWV2_downsample_seed,
    resources:
        cores       = lambda w, threads: threads
    log:       "logs/04-kinship/TKGWV2/TKGWV2_downsample_bam/{generation}/{pairA}_{pairB}.log"
    benchmark: "benchmarks/04-kinship/TKGWV2/TKGWV2_downsample_bam/{generation}/{pairA}_{pairB}.tsv"
    conda:     "../envs/TKGWV2.yml"
    threads:   1
    shell: """
        root_dir=`pwd`                                                   # Keep current dir in memory.
        ln -sfrt {params.workdir} {input.pairs} >  $root_dir/{log}       # temporary symlink
        cd {params.workdir}                     >> $root_dir/{log}       # go to output workdir
        TK-helpers.py downsampleBam \
        --downsampleN {params.downsampleN} \
        --downsampleSeed {params.seed}          >> $root_dir/{log} 2>&1  # run TK-helpers
        find . -maxdepth 1 -type l -delete      >> $root_dir/{log}       # delete symlinks.
    """


def get_target_frequencies():
    user_defined = config["kinship"]["TKGWV2"]["target-frequencies"]
    default      = config["netrules"]["TKGWV2"]["default-target-frequencies"]
    return user_defined or default

rule intersect_freq_file:
    input:
        frequencies   = get_target_frequencies(),
        targets       = get_snp_targets(ext=".snp"),
        intersect     = expand(rules.get_target_panel_intersect.output.targets, 
            maf      = config['variant-calling']['maf'],
            superpop = config['variant-calling']['maf-superpop'],
        ),
    output:
        frequencies   = "results/04-kinship/TKGWV2/{file}.frq".format(file = basename(splitext(get_target_frequencies())[0]))
    params:
    conda: "../envs/TKGWV2.yml"
    shell: """
        workflow/scripts/intersect-plink-frq.R {input.frequencies} {input.targets} {input.intersect} > {output.frequencies}
    """


def define_TKGWV2_input(wildcards):
    """
    Define the input for TKGWV2: either a downsampled bam file if the user
    requested it, or skip downsampling entirely.
    """
    if config['kinship']['TKGWV2']['downsample']:
        return rules.TKGWV2_downsample_bam.output
    else:
        return get_TKGWV2_input_bams(wildcards)

rule run_TKGWV2:
    """
    Perform kinship estimation on a single pair of individuals, using TKGWV2.
    See: Fernandes DM, et al. TKGWV2: an ancient DNA relatedness pipeline for 
         ultra-low coverage whole genome shotgun data. Sci Rep 11, 21262 (2021). 
         https://doi.org/10.1038/s41598-021-00581-3

    Repo: https://github.com/danimfernandes/tkgwv2.git
    
    @ TODO: Keep track of the 'frq' and 'tped' output files (TKGWV2 reorders pairA and pairB in lexical order.......)

    # Benchmarks: 
    | depth | max h:m:s | max_rss |
    | ----- | --------- | ------- |
    | 0.01X |           |         |
    | 0.05X | 0:05:33   | 3155    |
    | 0.10X | 0:04:37   | 7325    |
    | 1.00X | 0:25:17   | 19402   |
    """
    input:
        bams          = define_TKGWV2_input,
        reference     = ReferenceGenome.get_path(),
        bed_targets   = "data/TKGWV2/genomeWideVariants_hg19/1000GP3_22M_noFixed_noChr.bed",
        plink_targets = multiext("data/TKGWV2/genomeWideVariants_hg19/DummyDataset_EUR_22M_noFixed", ".bed", ".bim", ".fam"),
        frequencies   = rules.intersect_freq_file.output.frequencies
    output:
        results = "results/04-kinship/TKGWV2/{generation}/{pairA}_{pairB}/TKGWV2_Results.txt",
        peds    = temp([
            "results/04-kinship/TKGWV2/{generation}/{pairA}_{pairB}/{generation}_{pairA}.ped",
            "results/04-kinship/TKGWV2/{generation}/{pairA}_{pairB}/{generation}_{pairB}.ped"
        ]),
        maps    = temp([
            "results/04-kinship/TKGWV2/{generation}/{pairA}_{pairB}/{generation}_{pairA}.map",
            "results/04-kinship/TKGWV2/{generation}/{pairA}_{pairB}/{generation}_{pairB}.map"
        ]),
        pileups = temp([
            "results/04-kinship/TKGWV2/{generation}/{pairA}_{pairB}/{generation}_{pairA}.pileupsamtools.gwv.txt",
            "results/04-kinship/TKGWV2/{generation}/{pairA}_{pairB}/{generation}_{pairB}.pileupsamtools.gwv.txt"
        ]),
    params:
        plink_basename = lambda wildcards, input: splitext(input.plink_targets[0])[0],
        min_MQ     = config["kinship"]["TKGWV2"]["min-MQ"],
        min_BQ     = config["kinship"]["TKGWV2"]["min-BQ"],
        min_depth  = config["kinship"]["TKGWV2"]["min-depth"],
        bam_ext    = lambda wildcards, input: basename(input.bams[0]).split(".",1)[1]
    resources:
        runtime    = 10,
        #mem_mb     = 4000,
        mem_mb     = lambda w: round(2600 + 16900 * float(config['gargammel']['coverage'])),
        cores      = lambda w, threads: threads
    log:       "logs/04-kinship/TKGWV2/run_TKGWV2/{generation}/{pairA}_{pairB}.log"
    benchmark: "benchmarks/04-kinship/TKGWV2/run_TKGWV2/{generation}/{pairA}_{pairB}.tsv"
    conda:     "../envs/TKGWV2.yml"
    threads:   1
    shell: """
        base_dir=`pwd`                                     # Keep a record of the base directory
        cd $(dirname {output.results}) 2> $base_dir/{log}  # Go into the results directory

        # If the file is not present (i.e. no downsample has been made, create symlink.)
        for bam in {input.bams}; do
            [[ -f $(basename $bam) ]] || ln -sr $base_dir/$bam
        done 2>> $base_dir/{log}

        # Run TKGWV2
        TKGWV2.py bam2plink \
        --referenceGenome $base_dir/{input.reference} \
        --gwvList $base_dir/{input.bed_targets} \
        --gwvPlink $base_dir/{params.plink_basename} \
        --minMQ {params.min_MQ} \
        --minBQ {params.min_BQ} \
        --bamExtension .{params.bam_ext} \
        plink2tkrelated \
        --freqFile $base_dir/{input.frequencies} \
        --ignoreThresh {params.min_depth} \
        --verbose >> $base_dir/{log} 2>&1
    """


def define_TKGWV2_requested_dyads(wildcards):
    """
    Returns the list of requested pairwise comjparisons for TKGWV2,
    using the pedigree codes file. This will merge the results of each
    considered dyad into a single results file. one per generated pedigree.
    """

    # Open the pedigree codes definition file, and extract the relevant pedigree comparisons.
    with open(config["ped-sim"]["data"]["codes"]) as f:
        ind1, ind2 = zip(
            *[tuple(str.split(comparison.replace('\n', ''), '\t')[1:3]) for comparison in list(f.readlines())]
        )

        # Guard gen wildcard. We must do this in advance, or zipped expand will bug_out.
        out = expand(rules.run_TKGWV2.output.results,
            pairA="{pairA}",
            pairB="{pairB}",
            generation="{{generation}}"
        )
        relevant_comparisons = expand(
            out,
            zip,
            pairA=ind1,
            pairB=ind2
        )
        return relevant_comparisons

rule merge_TKGWV2_results:
    """
    Merge the pair-specific results of run_TKGWV2 into a single results file.
    (while removing header duplicates.)
    """
    input:
        results = define_TKGWV2_requested_dyads
    output:
        result  = "results/04-kinship/TKGWV2/{generation}/{generation}-TKGWV2_Results.txt"
    resources:
        cores   = lambda w, threads: threads
    log:   "logs/04-kinship/TKGWV2/merge_TKGWV2_results/{generation}.log"
    conda: "../envs/coreutils-9.1.yml"
    threads: 1
    shell: """
        awk 'FNR>1 || NR==1' {input.results} > {output.result} 2> {log}
    """
