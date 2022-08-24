
rule TKGWV2_downsample_bam:
    input:
        pairA = multiext("{preprocess_dir}/06-mapdamage/{gen}_{pairA}/{gen}_{pairA}.srt.rmdup.rescaled", ".bam", ".bam.bai")
        pairB = multiext("{preprocess_dir}/06-mapdamage/{gen}_{pairB}/{gen}_{pairB}.srt.rmdup.rescaled", ".bam", ".bam.bai")
    output:
        pairA = "out/03-kinship/TKGWV2/{gen}/{pairA}_{pairB}/{gen}_{pairA}.srt.rmdup.rescaled_subsampled.bam"
        pairB = "out/03-kinship/TKGWV2/{gen}/{pairA}_{pairB}/{gen}_{pairB}.srt.rmdup.rescaled_subsampled.bam"
    params:
        workdir = lambda: wildcards, output: basename(output.pairA)
    conda: "../envs/TKGWV2.yml"
    shadow: "minimal"
    shell: """
        ln -sfrt {params.workdir} {input.pairA} {input.pairB} # temporary symlink
        pushd {params.workdir}                                # go to output workdir
        Rscript downsampleBAM.R                               # run TKGWV2 helper-script for downsampling
        find . -maxdepth 1 -type l -delete                    # delete symlinks.
    """


rule TKGWV2_bam2plink:
    input:
        bams           = rules.TKGWV2_downsample_bam.output,
        reference     = config["refgen"]
        bed_targets   = lambda wildcards: splitext(config["kinship"]["targets"])[0] + ".ucscbed"
        plink_targets = lambda wildcards: splitext(config["kinship"]["targets"])[0]
        frequencies   = lambda wildcards: splitext(config["kinship"]["targets"])[0] + ".frq"
        
    output:
    conda: "../envs/TKGWV2.yml"
    params:
        min_MQ    = config["kinship"]["TKGWV2"]["min-MQ"],
        min_BQ    = config["kinship"]["TKGWV2"]["min-BQ"],
        min_depth = config["kinship"]["TKGWV2"]["min-depth"]
    shell: """
        TKGWV2.py bam2plink --referenceGenome {input.reference}    \
                            --gwvList {input.targets}              \
                            --gwvPlink {input.plink_targets}       \
                            --minMQ {params.min_MQ}                \
                            --minBQ {params.min_BQ}                \
                  plink2tkrelated --freqFile {input.frequencies}   \
                                  --ignoreTresh {params.min_depth} \
                                  --verbose
    """