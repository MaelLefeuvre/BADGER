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

        out = "results/04-kinship/TKGWV2/ped{gen}/{{pairA}}_{{pairB}}/TKGWV2_Results.txt"

        relevant_comparisons = expand(expand(out, gen=range(1, gen_no+1)), zip, pairA=ind1, pairB=ind2)
        return relevant_comparisons


rule TKGWV2_download_support_files:
    output:
        support_files = expand("{directory}/{dataset}", 
            directory = config['kinship']['TKGWV2']['support-files-dir'],
            dataset = [
                "1000GP3_EUR_1240K.frq",
                "1000GP3_22M_noFixed_noChr.bed",
                "DummyDataset_EUR_22M_noFixed.bed",
                "DummyDataset_EUR_22M_noFixed.bim",
                "DummyDataset_EUR_22M_noFixed.fam"
            ]
        )
    params:
        url = config['kinship']['TKGWV2']['support-files-url'],
        output_dir = config['kinship']['TKGWV2']['support-files-dir']
    conda: "../envs/gdown-4.4.0.yml"
    shell: """
        gdown "{params.url}" -O {params.output_dir} --folder
    """


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
        pairA = "results/04-kinship/TKGWV2/{gen}/{pairA}_{pairB}/{gen}_{pairA}.srt.rmdup.rescaled_subsampled.bam",
        pairB = "results/04-kinship/TKGWV2/{gen}/{pairA}_{pairB}/{gen}_{pairB}.srt.rmdup.rescaled_subsampled.bam"
    params:
        workdir     = lambda wildcards, output: dirname(output.pairA),
        downsampleN = config['kinship']['TKGWV2']['downsample-N']
    log: "logs/04-kinship/TKGWV2/{gen}/{pairA}_{pairB}-downsample_bam.log"
    conda: "../envs/TKGWV2.yml"
    shell: """
        root_dir=`pwd`                                                                         # Keep current dir in memory.
        ln -sfrt {params.workdir} {input.pairA} {input.pairB}                                  # temporary symlink
        cd {params.workdir}                                                                    # go to output workdir
        TK-helpers.py downsampleBam --downsampleN {params.downsampleN} > $root_dir/{log} 2>&1  # run TK-helpers
        find . -maxdepth 1 -type l -delete                                                     # delete symlinks.
    """


def define_TKGWV2_input(wildcards):
    if config['kinship']['TKGWV2']['downsample']:
        return rules.TKGWV2_downsample_bam.output
    else:
        return rules.TKGWV2_downsample_bam.input

rule run_TKGWV2:
    input:
        bams          = define_TKGWV2_input,
        reference     = config["reference"],
        bed_targets   = "data/TKGWV2/1000GP3_22M_noFixed_noChr.bed",
        plink_targets = multiext("data/TKGWV2/DummyDataset_EUR_22M_noFixed", ".bed", ".bim", ".fam"),
        frequencies   = config['kinship']['TKGWV2']['target-frequencies']
    output:
        results       = "results/04-kinship/TKGWV2/{gen}/{pairA}_{pairB}/TKGWV2_Results.txt"
    conda: "../envs/TKGWV2.yml"
    params:
        plink_basename = lambda wildcards: splitext(config["kinship"]["targets"])[0],
        min_MQ     = config["kinship"]["TKGWV2"]["min-MQ"],
        min_BQ     = config["kinship"]["TKGWV2"]["min-BQ"],
        min_depth  = config["kinship"]["TKGWV2"]["min-depth"],
        bam_ext    = lambda wildcards, input: basename(input.bams[0]).split(".",1)[1]
    log: "logs/04-kinship/TKGWV2/{gen}/{pairA}_{pairB}-run_TKGWV2.log"
    shell: """
        base_dir=`pwd`                                                         # Keep a record of the base directory
        cd $(dirname {output.results})                                         # Go into the results directory

        # If the file is not present (i.e. no downsample has been made, create symlink.)
        for bam in {input.bams}; do
            [[ -f $(basename $bam) ]] || ln -s $bam
        done

        # Run TKGWV2
        TKGWV2.py bam2plink --referenceGenome $base_dir/{input.reference}    \
                            --gwvList $base_dir/{input.bed_targets}          \
                            --gwvPlink $base_dir/{params.plink_basename}     \
                            --minMQ {params.min_MQ}                          \
                            --minBQ {params.min_BQ}                          \
                            --bamExtension .{params.bam_ext}                 \
                  plink2tkrelated --freqFile $base_dir/{input.frequencies}   \
                                  --ignoreThresh {params.min_depth}          \
                                  --verbose > $base_dir/{log} 2>&1
    """