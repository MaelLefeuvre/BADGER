configfile: "config/config.yml"

localrules: create_archive_metadata

import os
import re

def define_archived_bams(wildcards):
    archived_bams = define_dedup_input_bam(wildcards)
    samples       = get_samples_ids(wildcards)
    return expand(expand(archived_bams, sample=samples), generation="ped{gen}")


def define_run_number():
    base_dir = config['archive']['archive-dir']
    os.makedirs(base_dir, exist_ok=True)

    regex          = re.compile("run-([0-9])+")
    found_archives = list(filter(regex.match, os.listdir(base_dir)))
    archive_runs   = sorted( [int(archive.removeprefix("run-")) for archive in found_archives], reverse=True)
    next_archive   = 0 if len(archive_runs) == 0 else archive_runs[0] + 1
    return f"{next_archive:03d}"
    #return next_archive


rule create_archive_metadata:
    output: f"{config['archive']['archive-dir']}/archive-metadata.yml"
    shell: """
        touch {output}
    """


rule archive_config_file:
    input:
        config_file = "config/config.yml"
    output:
        archive = protected("{archive_dir}/run-{{run}}/{file}".format(
            archive_dir = config['archive']['archive-dir'],
            file       = "config/config.yml"
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{file}.md5sum".format(
            archive_dir = config['archive']['archive-dir'],
            file       = "config/config.yml"
        )),
    log: "logs/07-archive/run{run}/archive_pipeline_metadata.log"
    shell: """
        cp {input.config_file} {output.archive} > {log} 2>&1

        # store checksums in metadata file yaml-like format.
        md5sum {input.config_file} > {output.checksum} 2>> {log}
    """


rule archive_pipeline_metadata:
    input:
        metadata = "results/meta/pipeline-metadata.yml"
    output:
        archive = protected("{archive_dir}/run-{{run}}/{file}".format(
            archive_dir = config['archive']['archive-dir'],
            file       = "results/meta/pipeline-metadata.yml"
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{file}.md5sum".format(
            archive_dir = config['archive']['archive-dir'],
            file       = "results/meta/pipeline-metadata.yml"
        )),
    log: "logs/07-archive/run{run}/archive_pipeline_metadata.log"
    shell: """
        cp {input.metadata} {output.archive} > {log} 2>&1

        # store checksums in metadata file yaml-like format.
        md5sum {input.metadata} > {output.checksum} 2>> {log}
    """

rule archive_ped_sim_vcf:
    input:
        merged_vcf = rules.dopplegang_twins.output.merged_vcf.format(POP=config['ped-sim']['params']['pop'])
    output:
        archive = protected("{archive_dir}/run-{{run}}/{ped_sim_dir}.tar.xz".format(
            archive_dir = config['archive']['archive-dir'],
            ped_sim_dir = os.path.dirname(rules.dopplegang_twins.output.merged_vcf)
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{ped_sim_dir}.md5sums".format(
            archive_dir = config['archive']['archive-dir'],
            ped_sim_dir = os.path.dirname(rules.dopplegang_twins.output.merged_vcf)
        ))
    params:
        target_directory = os.path.dirname(rules.dopplegang_twins.output.merged_vcf),
        compress_level=config['archive']['compress-level']
    log: "logs/07-archive/run{run}/archive_ped_sim_vcf.log"
    threads: 16
    shell: """
        XZ_OPT="-{params.compress_level} -T{threads}" tar -cJvf {output.archive} {params.target_directory} > {log} 2>&1 

        # store checksums in metadata file yaml-like format.
        md5sum $(dirname {input.merged_vcf})/* > {output.checksum} 2>> {log}
    """


rule archive_pedigree_bams:
    """
    Archives an entire pedigree replicate into the cram format.
    NOTE: To retrive the archive, "samtools split -f "%!".bam --output-fmt BAM ped1-merged.cram"
    """
    input: 
        bams      = define_archived_bams,
        reference = config["reference"],
        metadata  = rules.create_archive_metadata.output,
    output:
        cram = protected("{archive_dir}/run-{{run}}/{results_dir}/ped{{gen}}-merged.cram".format(
            archive_dir = config['archive']['archive-dir'],
            results_dir = "/".join(os.path.dirname(define_dedup_input_bam([])).strip("/").split('/')[:-1])
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{results_dir}/ped{{gen}}-merged.md5sums".format(
            archive_dir = config['archive']['archive-dir'],
            results_dir = "/".join(os.path.dirname(define_dedup_input_bam([])).strip("/").split('/')[:-1])
        ))
    params:
        opts = "seqs_per_slice=100000,level={level},use_lzma,use_fqz,use_arith".format(
            level=config['archive']['compress-level']
        )
    log: "logs/07-archive/run{run}/ped{gen}/archive_pedigree_bams.log"
    conda: "../envs/samtools-1.15.yml"
    threads: 16
    shell: """
        THREADS=`echo {threads}/2 | bc`
        samtools merge -@ ${{THREADS}} -o - {input.bams} | samtools view -@ ${{THREADS}} -OCRAM -T {input.reference} --output-fmt-option {params.opts} > {output.cram} 2> {log}
        
        # store checksums in metadata file yaml-like format.
        md5sum {input.bams} > {output.checksum} 2>> {log}
    """


rule archive_variant_callset:
    input:
        panel = rules.get_target_panel_intersect.output.targets

    output:
        archive = protected("{archive_dir}/run-{{run}}/{panel}.xz".format(
            archive_dir = config['archive']['archive-dir'],
            panel       = rules.get_target_panel_intersect.output.targets
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{panel}.md5sum".format(
            archive_dir = config['archive']['archive-dir'],
            panel       = rules.get_target_panel_intersect.output.targets
        )),
    params:
        compress_level=config['archive']['compress-level']
    log: "logs/07-archive/run{run}/archive_variant_callset.log"
    threads: 4
    shell: """
        xz -zc -{params.compress_level} -e -T {threads} {input.panel} > {output.archive} 2> {log} 2>&1 

        # store checksums in metadata file yaml-like format.
        md5sum {input.panel} > {output.checksum} 2>> {log}
    """

rule archive_contaminants:
    input:
        contaminants = rules.get_contamination_table.output.cont_table
    output:
        archive = protected("{archive_dir}/run-{{run}}/{file}".format(
            archive_dir = config['archive']['archive-dir'],
            file       = rules.get_contamination_table.output.cont_table
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{file}.md5sum".format(
            archive_dir = config['archive']['archive-dir'],
            file       = rules.get_contamination_table.output.cont_table
        )),
    log: "logs/07-archive/run{run}/archive_contaminants.log"
    shell: """
        cp {input.contaminants} {output.archive} > {log} 2>&1

        # store checksums in metadata file yaml-like format.
        md5sum {input.contaminants} > {output.checksum} 2>> {log}
    """


def make_snakemake_happy_with_READ_checkpoint(wildcards):
    with checkpoints.get_samples.get().output[0].open() as f:
        return [
            rules.run_READ.output.results,
            rules.run_READ.output.means,
            rules.run_READ.output.raw,
            rules.run_READ.output.plot
        ]

rule archive_READ_results:
    input:
        results = make_snakemake_happy_with_READ_checkpoint,
        #results  = [
        #    "results/04-kinship/READ/{generation}/READ_results",
        #    "results/04-kinship/READ/{generation}/meansP0_AncientDNA_normalized",
        #    "results/04-kinship/READ/{generation}/READ_output_ordered",
        #    "results/04-kinship/READ/{generation}/READ_results_plot.pdf"
        #],
    output:
        archive = protected("{archive_dir}/run-{{run}}/{READ_dir}.tar.xz".format(
            archive_dir = config['archive']['archive-dir'],
            READ_dir = "results/04-kinship/READ/{generation}"
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{READ_dir}.md5sums".format(
            archive_dir = config['archive']['archive-dir'],
            READ_dir = "results/04-kinship/READ/{generation}"
        ))
    params:
        target_directory = "results/04-kinship/READ/{generation}",
        compress_level=config['archive']['compress-level']
    log: "logs/07-archive/run{run}/{generation}/archive_READ_results.log"
    threads: 4
    shell: """
        XZ_OPT="-{params.compress_level} -T{threads}" tar -cJvf {output.archive} {params.target_directory} > {log} 2>&1 

        # store checksums in metadata file yaml-like format.
        md5sum {params.target_directory}/* > {output.checksum} 2>> {log}
    """


def make_snakemake_happy_with_GRUPS_checkpoint(wildcards):
    with checkpoints.get_samples.get().output[0].open() as f:
        return rules.run_GRUPS.output.output_dir

rule archive_GRUPS_results:
    input:
        results_dir = make_snakemake_happy_with_GRUPS_checkpoint
        #results_dir = "results/04-kinship/GRUPS/{generation}/{generation}.results",
    output:
        archive = protected("{archive_dir}/run-{{run}}/{GRUPS_dir}.tar.xz".format(
            archive_dir = config['archive']['archive-dir'],
            GRUPS_dir = "results/04-kinship/GRUPS/{generation}"
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{GRUPS_dir}.md5sums".format(
            archive_dir = config['archive']['archive-dir'],
            GRUPS_dir = "results/04-kinship/GRUPS/{generation}"
        ))
    params:
        target_directory = "results/04-kinship/GRUPS/{generation}",
        compress_level   = config['archive']['compress-level']
    log: "logs/07-archive/run{run}/{generation}/archive_GRUPS_results.log"
    threads: 4
    shell: """
        XZ_OPT="-{params.compress_level} -T{threads}" tar -cJvf {output.archive} {params.target_directory} > {log} 2>&1 

        # store checksums in metadata file yaml-like format.
        find {params.target_directory} -type f -exec md5sum {{}} \; > {output.checksum} 2>> {log}
    """


rule archive_TKGWV2_results:
    input:
        results = rules.merge_TKGWV2_results.output.result
    output:
        archive = protected("{archive_dir}/run-{{run}}/{file}".format(
            archive_dir = config['archive']['archive-dir'],
            file       = rules.merge_TKGWV2_results.output.result
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{file}.md5sum".format(
            archive_dir = config['archive']['archive-dir'],
            file       = rules.merge_TKGWV2_results.output.result
        )),
    log: "logs/07-archive/run{run}/ped{gen}/archive_TKGWV2_results.log"
    shell: """
        cp {input.results} {output.archive} > {log} 2>&1

        # store checksums in metadata file yaml-like format.
        md5sum {input.results} > {output.checksum} 2>> {log}
    """

rule check_duplicate_archives:
    """
    Update our global checksums file, with those of the current run. This rule has two objectives:
    1. Keep an automatic record of the MD5sums of each archived file.
    2. Automatically check if there are duplicated files within our archive folders.
    """
    input:
        metadata  = rules.create_archive_metadata.output,

        # ---- Output bam files.
        bam_checksums = expand(rules.archive_pedigree_bams.output.checksum,
            gen    = range(1, config['ped-sim']['replicates'] + 1),
            run = "{run}"
        ),

        # ---- Pipeline config file
        pipeline_config = rules.archive_config_file.output.checksum,

        # ---- Pipeline metadata
        pipeline_meta_checksum = rules.archive_pipeline_metadata.output.checksum,

        # ---- Pedigree simulations
        ped_sim_checksum = rules.archive_ped_sim_vcf.output.checksum,

        # ---- Gargammel contaminant table
        cont_checksum = rules.archive_contaminants.output.checksum,

        # ---- Variant calling callset.
        callset_checksum = rules.archive_variant_callset.output.checksum,

        # ---- READ results callset.
        READ_checksums = expand(rules.archive_READ_results.output.checksum,
            generation    = ["ped" + str(rep) for rep in range(1, config['ped-sim']['replicates'] + 1)],
            run = "{run}"
        ),

        # ---- GRUPS results callset.
        GRUPS_checksums = expand(rules.archive_GRUPS_results.output.checksum,
            generation    = ["ped" + str(rep) for rep in range(1, config['ped-sim']['replicates'] + 1)],
            run = "{run}"
        ),

        # ---- TKGWV2 results callset.
        TKGWV2_checksums = expand(rules.archive_TKGWV2_results.output.checksum,
            gen    = range(1, config['ped-sim']['replicates'] + 1),
            run = "{run}"
        ),
    output:
        metadata = temp(touch(f"{config['archive']['archive-dir']}/run-{{run}}-check_duplicate_archives.done"))
    run:
        import sys
        import yaml
        import hashlib

        def extend_checksums(global_checksums, checksum_file, run_key, subkey):
            with open(checksum_file, "r") as f:
                run_checksums = {bam: md5 for md5, bam in [line.rstrip().split() for line in f]}
            global_checksums[run_key][subkey].update(run_checksums)

        # ---- Load the global checksum yaml
        metadata_file = input.metadata[0]
        with open(metadata_file, "r") as yml:
            global_checksums = yaml.load(yml, Loader=yaml.FullLoader)
        
        # ---- Add this run's id as a key. If the key is already present, this is an error.
        run_key = f"run-{wildcards.run}"
        try:
            if run_key in global_checksums:
                raise RuntimeError(f"{run_key} is already present in {metadata_file}")
        except TypeError as e:
            # NoneType is not iterable -> Probably means our file is empty.
            print(f"WARN: {metadata_file} appears to be empty: Got '{e}'", file=sys.stderr)
            global_checksums = dict()
        finally:
            global_checksums[run_key] = {
                "config-file"     : dict(),
                "metadata"        : dict(),
                "bams"            : dict(),
                "ped-sim"         : dict(),
                "variant-callset" : dict(),
                "contaminants"    : dict(),
                "READ"            : dict(),
                "GRUPS"           : dict(),
                "TKGWV2"          : dict(),
            }
        
        # Extend our global checksum file w/ the bam checksums of this run.
        for checksum in input.bam_checksums:
            extend_checksums(global_checksums, checksum, run_key, "bams")

        # Warn if some archive / runs are duplicates...
        for run in global_checksums.keys():
            for subkey in global_checksums[run]:
                if global_checksums[run][subkey] == global_checksums[run_key][subkey] and run != run_key:
                    print(f"WARN: Archive {subkey} of {run} appears to be the same as the one currently being worked on: {run_key}", file=sys.stderr)

        # Extend our global checksum file w/ the ped-sim checksums of this run.
        extend_checksums(global_checksums, input.pipeline_config, run_key, "config-file")

        # Extend our global checksum file w/ the ped-sim checksums of this run.
        extend_checksums(global_checksums, input.pipeline_meta_checksum, run_key, "metadata")

        # Extend our global checksum file w/ the ped-sim checksums of this run.
        extend_checksums(global_checksums, input.ped_sim_checksum, run_key, "ped-sim")

        # Extend our global checksum file w/ the contamination table of this run.
        extend_checksums(global_checksums, input.cont_checksum, run_key, "contaminants")

        # Extend our global checksum file w/ the variant callset of this run.
        extend_checksums(global_checksums, input.callset_checksum, run_key, "variant-callset")

        # Extend our global checksum file w/ each pedigree's READ result.
        for checksum in input.READ_checksums:
            extend_checksums(global_checksums, checksum, run_key, "READ")

        # Extend our global checksum file w/ each pedigree's GRUPS result.
        for checksum in input.READ_checksums:
            extend_checksums(global_checksums, checksum, run_key, "GRUPS")

        # Extend our global checksum file w/ each pedigree's GRUPS result.
        for checksum in input.TKGWV2_checksums:
            extend_checksums(global_checksums, checksum, run_key, "TKGWV2")

        with open(input.metadata[0], "w") as yml:
            documents = yaml.dump(global_checksums, yml)


rule archive:
    input: 
        # run number must be defined ONCE, at the final rule, or the function will be executed multiple times.
        # Caching does not seem to resolve the trick.
        # @ TODO: READ & GRUPS Checkpoints is causing rule re-triggering...
        checksums =  expand(rules.check_duplicate_archives.output, run=define_run_number())