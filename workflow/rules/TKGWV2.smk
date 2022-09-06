from os.path import splitext, basename
def TKGWV2_output(wildcards):
    """
    Returns the list of requires pairwise comjparisons for TKGWV2
    """

    # Get the number of expected generations.
    with open(rules.run_ped_sim.input.definition) as f:
        gen_no = int(list(dropwhile(lambda x: x.startswith("#"), [line for line in f]))[0].split(' ')[2])

    # Open the pedigree codes definition file, and extract the relevant pedigree comparisons.
    with open(config["ped-sim"]["data"]["codes"]) as f:
        ind1, ind2 = zip(*[tuple(str.split(comparison.replace('\n', ''), '\t')[1:3]) for comparison in list(f.readlines())])

        out = "results/03-kinship/TKGWV2/ped{gen}/{{pairA}}_{{pairB}}/TKGWV2_Results.txt"

        relevant_comparisons = expand(expand(out, gen=range(1, gen_no+1)), zip, pairA=ind1, pairB=ind2)
        return relevant_comparisons


rule TKGWV2_downsample_bam:
    """
    Downsample the .bam file if there are more than 1_500_000 reads.
     - As is, this helper script is applied to all files ending with the ".bam" suffix within the working directory.
     - output files are identified with the "_subsampled.bam" suffix.
    """
    input:
        pairA = multiext("results/02-preprocess/06-mapdamage/{gen}_{pairA}/{gen}_{pairA}.srt.rmdup.rescaled", ".bam", ".bam.bai"),
        pairB = multiext("results/02-preprocess/06-mapdamage/{gen}_{pairB}/{gen}_{pairB}.srt.rmdup.rescaled", ".bam", ".bam.bai")
    output:
        pairA = "results/03-kinship/TKGWV2/{gen}/{pairA}_{pairB}/{gen}_{pairA}.srt.rmdup.rescaled_subsampled.bam",
        pairB = "results/03-kinship/TKGWV2/{gen}/{pairA}_{pairB}/{gen}_{pairB}.srt.rmdup.rescaled_subsampled.bam"
    params:
        workdir = lambda wildcards, output: dirname(output.pairA)
    conda: "../envs/TKGWV2.yml"
    shell: """
        ln -sfrt {params.workdir} {input.pairA} {input.pairB} # temporary symlink
        pushd {params.workdir}                                # go to output workdir
        Rscript $(which downsampleBam.R)                      # run TKGWV2 helper-script for downsampling
        find . -maxdepth 1 -type l -delete                    # delete symlinks.
    """


rule run_TKGWV2:
    input:
        bams           = rules.TKGWV2_downsample_bam.output,
        reference     = config["refgen"],
        bed_targets   = lambda wildcards: splitext(config["kinship"]["targets"])[0] + ".ucscbed",
        plink_targets = lambda wildcards: multiext(splitext(config["kinship"]["targets"])[0] + "_onesample", ".bed", ".bim", ".fam"),
        frequencies   = lambda wildcards: splitext(config["kinship"]["targets"])[0] + ".frq"
    output:
        results       = "results/03-kinship/TKGWV2/{gen}/{pairA}_{pairB}/TKGWV2_Results.txt"
    conda: "../envs/TKGWV2.yml"
    params:
        plink_basename = lambda wildcards: splitext(config["kinship"]["targets"])[0],
        min_MQ     = config["kinship"]["TKGWV2"]["min-MQ"],
        min_BQ     = config["kinship"]["TKGWV2"]["min-BQ"],
        min_depth  = config["kinship"]["TKGWV2"]["min-depth"],
        bam_ext    = ".srt.rmdup.rescaled_subsampled.bam"
    log: "logs/03-kinship/TKGWV2/run_TKGWV2/{gen}/run_TKGWV2-{pairA}_{pairB}.log"
    shell: """
        base_dir=`pwd`                                                         # Keep a record of the base directory
        cd $(dirname {output.results})                                         # Go into the results directory
        TKGWV2.py bam2plink --referenceGenome $base_dir/{input.reference}    \
                            --gwvList $base_dir/{input.bed_targets}          \
                            --gwvPlink $base_dir/{params.plink_basename}     \
                            --minMQ {params.min_MQ}                          \
                            --minBQ {params.min_BQ}                          \
                            --bamExtension {params.bam_ext}                  \
                  plink2tkrelated --freqFile $base_dir/{input.frequencies}   \
                                  --ignoreThresh {params.min_depth}          \
                                  --verbose 2>&1 {log}
    """