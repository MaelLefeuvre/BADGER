configfile: "config/config.yml"

localrules: create_archive_metadata

import os, re

# ---- Utility functions.
def define_archived_bams(wildcards):
    archived_bams = define_dedup_input_bam(wildcards)
    samples       = get_samples_ids(wildcards)
    return expand(expand(archived_bams, sample=samples), generation="ped{gen}")


def define_run_number(wildcards, run_prefix="run-"):
    base_dir = config['archive']['archive-dir']
    os.makedirs(base_dir, exist_ok=True)
    regex          = re.compile("run-([0-9])+")
    found_archives = list(filter(regex.match, os.listdir(base_dir)))
    archive_runs   = sorted([int(archive.removeprefix(run_prefix)) for archive in found_archives], reverse=True)
    (last_archive, next_archive) = (archive_runs[0], archive_runs[0] + 1) if len(archive_runs) > 0 else (None, 0)

    # ---- Load the global checksum yaml and check if the last run was completed.
    # Prompt the user if he wishes to overwrite the previous run.
    yml = checkpoints.create_archive_metadata.get().output[0]
    global_checksums = yaml.load(open(yml, "r"), Loader=yaml.FullLoader)
    global_checksums = global_checksums if global_checksums is not None else []
    format_run = lambda x: f"{run_prefix}{x:03d}"
    if last_archive is not None and format_run(last_archive) not in global_checksums:
        while True:
            prompt = input((
                "Last run appears to be incomplete "
                f"({format_run(last_archive)}). "
                "Do you wish to overwrite it? [Yes/no]: "
            ))
            match prompt.lower():
                case "yes":
                    next_archive = last_archive
                    break
                case "no" | "n":
                    next_archive = next_archive
                    break
                case _:
                    continue
    return f"{next_archive:03d}"

# ------------------------------------------------------------------------------------------------ #

checkpoint create_archive_metadata:
    output: f"{config['archive']['archive-dir']}/archive-metadata.yml"
    shell: """
        touch {output}
    """


rule archive_config_file:
    input:
        config_file = "config/config.yml"
    output:
        archive     = protected("{archive_dir}/run-{{run}}/{file}".format(
            archive_dir=config['archive']['archive-dir'],
            file       ="config/config.yml"
        )),
        checksum    = protected("{archive_dir}/run-{{run}}/{file}.md5sum".format(
            archive_dir=config['archive']['archive-dir'],
            file       ="config/config.yml"
        )),
    resources:
        cores       = lambda w, threads: threads
    log:       "logs/07-archive/run{run}/archive_pipeline_metadata.log"
    threads:   1
    shell: """
        cp {input.config_file} {output.archive} > {log} 2>&1

        # store checksums in metadata file yaml-like format.
        md5sum {input.config_file} > {output.checksum} 2>> {log}
    """


rule archive_pipeline_metadata:
    """
    Archive the global pipeline metadata file (containing backup seeds + commit hashes.)
    """
    input:
        metadata = "results/meta/pipeline-metadata.yml"
    output:
        archive  = protected("{archive_dir}/run-{{run}}/{file}".format(
            archive_dir=config['archive']['archive-dir'],
            file       ="results/meta/pipeline-metadata.yml"
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{file}.md5sum".format(
            archive_dir=config['archive']['archive-dir'],
            file       ="results/meta/pipeline-metadata.yml"
        )),
    resources:
        cores    = lambda w, threads: threads
    log:       "logs/07-archive/run{run}/archive_pipeline_metadata.log"
    benchmark: "benchmarks/07-archive/run{run}/archive_pipeline_metadata.tsv"
    threads:   1
    shell: """
        cp {input.metadata} {output.archive} > {log} 2>&1

        # store checksums in metadata file yaml-like format.
        md5sum {input.metadata} > {output.checksum} 2>> {log}
    """


rule archive_ped_sim_vcf:
    """
    Archive ped-sim's pedigree simulations vcf for the entire run. Compressed in lzma1
    """
    input:
        merged_vcf = rules.merge_ped_sim.output.vcf.format(POP=config['ped-sim']['params']['pop'])
    output:
        archive    = protected("{archive_dir}/run-{{run}}/{ped_sim_dir}.tar.xz".format(
            archive_dir=config['archive']['archive-dir'],
            ped_sim_dir=os.path.dirname(rules.merge_ped_sim.output.vcf)
        )),
        checksum   = protected("{archive_dir}/run-{{run}}/{ped_sim_dir}.md5sums".format(
            archive_dir=config['archive']['archive-dir'],
            ped_sim_dir=os.path.dirname(rules.merge_ped_sim.output.vcf)
        ))
    params:
        target_directory = os.path.dirname(rules.merge_ped_sim.output.vcf),
        compress_level=config['archive']['compress-level']
    resources:
        cores     = lambda w, threads: threads
    log:       "logs/07-archive/run{run}/archive_ped_sim_vcf.log"
    benchmark: "benchmarks/07-archive/run{run}/archive_ped_sim_vcf.tsv"
    threads:   16
    shell: """
        XZ_OPT="-{params.compress_level} -T{threads}" tar -cJvf {output.archive} $(find {params.target_directory} -maxdepth 1 -type f) > {log} 2>&1 

        # store checksums in metadata file yaml-like format.
        # md5sum $(dirname {input.merged_vcf})/* > {output.checksum} 2>> {log}
        md5sum $(find $(dirname {input.merged_vcf}) -maxdepth 1 -type f) > {output.checksum} 2>> {log}
    """


rule archive_pedigree_bams:
    """
    Archives an entire pedigree replicate into the cram format.
    NOTE: To retrive the archive, "samtools split -f "%!".bam --output-fmt BAM ped1-merged.cram"
    """
    input: 
        bams      = define_archived_bams,
        reference = ReferenceGenome.get_path(),
        metadata  = rules.create_archive_metadata.output,
    output:
        cram      = protected("{archive_dir}/run-{{run}}/{results_dir}/ped{{gen}}-merged.cram".format(
            archive_dir=config['archive']['archive-dir'],
            results_dir="/".join(os.path.dirname(define_dedup_input_bam([])).strip("/").split('/')[:-1])
        )),
        checksum  = protected("{archive_dir}/run-{{run}}/{results_dir}/ped{{gen}}-merged.md5sums".format(
            archive_dir=config['archive']['archive-dir'],
            results_dir="/".join(os.path.dirname(define_dedup_input_bam([])).strip("/").split('/')[:-1])
        ))
    params:
        opts      = "seqs_per_slice=100000,level={level},use_lzma,use_fqz,use_arith".format(
            level=config['archive']['compress-level']
        )
    resources:
        cores = lambda w, threads: threads
    log:       "logs/07-archive/run{run}/ped{gen}/archive_pedigree_bams.log"
    benchmark: "benchmarks/07-archive/run{run}/ped{gen}/archive_pedigree_bams.tsv"
    conda:     "../envs/samtools-1.15.yml"
    threads:   16
    shell: """
        THREADS=`echo {threads}/2 | bc`
        samtools merge -@ ${{THREADS}} -o - {input.bams} | samtools view -@ ${{THREADS}} -OCRAM -T {input.reference} --output-fmt-option {params.opts} > {output.cram} 2> {log}
        
        # store checksums in metadata file yaml-like format.
        md5sum {input.bams} > {output.checksum} 2>> {log}
    """


rule archive_mapdamage_misincorporation:
    input:
        misincorporation = rules.run_mapdamage.output.misincorporation
    output:
        archive = protected("{archive_dir}/run-{{run}}/{file}".format(
            archive_dir = config['archive']['archive-dir'],
            file        = rules.run_mapdamage.output.misincorporation
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{file}.md5sum".format(
            archive_dir = config['archive']['archive-dir'],
            file        = rules.run_mapdamage.output.misincorporation
        ))
    resources:
        cores = lambda w, threads: threads
    log:     "logs/07-archive/run{run}/{sample}-archive_contaminants.log"
    threads: 1
    shell: """
        cp {input.misincorporation} {output.archive} > {log} 2>&1
        
        # store checksums in metadata file yaml-like format
        md5sum {input.misincorporation} > {output.checksum} 2>> {log}
    """


rule archive_variant_callset:
    """
    Archive this run's interset panel.
    """
    input:
        panel    = rules.samtools_pileup.input.targets

    output:
        archive  = protected("{archive_dir}/run-{{run}}/{panel}.xz".format(
            archive_dir=config['archive']['archive-dir'],
            panel      =rules.samtools_pileup.input.targets
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{panel}.md5sum".format(
            archive_dir=config['archive']['archive-dir'],
            panel      =rules.samtools_pileup.input.targets
        )),
    params:
        compress_level = config['archive']['compress-level']
    resources:
        cores         = lambda w, threads: threads
    log:       "logs/07-archive/run{run}/archive_variant_callset.log"
    benchmark: "benchmarks/07-archive/run{run}/archive_variant_callset.tsv"
    threads:   4
    shell: """
        xz -zc -{params.compress_level} -e -T {threads} {input.panel} > {output.archive} 2> {log} 2>&1 

        # store checksums in metadata file yaml-like format.
        md5sum {input.panel} > {output.checksum} 2>> {log}
    """


rule archive_contaminants:
    """
    Archive this run's contamination table (assigning 1000g sample-ids for each pedigree.)
    """
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
    resources:
        cores = lambda w, threads: threads
    log:       "logs/07-archive/run{run}/archive_contaminants.log"
    threads:   1
    shell: """
        cp {input.contaminants} {output.archive} > {log} 2>&1

        # store checksums in metadata file yaml-like format.
        md5sum {input.contaminants} > {output.checksum} 2>> {log}
    """


def check_READ_checkpoint(wildcards):
    """
    Run the checkpoint for READ
    This additional redirection is required to force checkpoint evaluation...
    """
    with checkpoints.get_samples.get().output[0].open() as f:
        return [
            rules.run_READ.output.results,
            rules.run_READ.output.means,
            rules.run_READ.output.raw,
            rules.run_READ.output.plot
        ]

rule archive_READ_results:
    """
    Create a generation-wise archive of READ's result. Results are compressed in lzma format.
    """
    input:
        results = check_READ_checkpoint,
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
    resources:
        cores = lambda w, threads: threads
    log:       "logs/07-archive/run{run}/{generation}/archive_READ_results.log"
    benchmark: "benchmarks/07-archive/run{run}/{generation}/archive_READ_results.tsv"
    threads:   4
    shell: """
        XZ_OPT="-{params.compress_level} -T{threads}" tar -cJvf {output.archive} {params.target_directory} > {log} 2>&1 

        # store checksums in metadata file yaml-like format.
        md5sum {params.target_directory}/* > {output.checksum} 2>> {log}
    """

def check_READv2_checkpoint(wildcards):
    """
    Run the checkpoint of READv2
    This additional redirection is required to force checkpoint evaluation...
    """
    with checkpoints.get_samples.get().output[0].open() as f:
        return [
            rules.run_READv2.output.results,
            rules.run_READv2.output.meansP0
        ]

rule archive_READv2_results:
    """
    Create a generation-wise archive of READ's result. Results are compressed in lzma format.
    """
    input:
        results = check_READv2_checkpoint,
    output:
        archive = protected("{archive_dir}/run-{{run}}/{READv2_dir}.tar.xz".format(
            archive_dir = config['archive']['archive-dir'],
            READv2_dir = "results/04-kinship/READv2/{generation}"
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{READv2_dir}.md5sums".format(
            archive_dir = config['archive']['archive-dir'],
            READv2_dir = "results/04-kinship/READv2/{generation}"
        ))
    params:
        target_directory = "results/04-kinship/READv2/{generation}",
        compress_level=config['archive']['compress-level']
    resources:
        cores = lambda w, threads: threads
    log:       "logs/07-archive/run{run}/{generation}/archive_READ_results.log"
    benchmark: "benchmarks/07-archive/run{run}/{generation}/archive_READ_results.tsv"
    threads:   4
    shell: """
        XZ_OPT="-{params.compress_level} -T{threads}" tar -cJvf {output.archive} {params.target_directory} > {log} 2>&1 

        # store checksums in metadata file yaml-like format.
        md5sum {params.target_directory}/* > {output.checksum} 2>> {log}
    """


def check_GRUPS_checkpoint(wildcards):
    """
    Run the checkpoint for GRUPS-rs
    This additional redirection is required to force checkpoint evaluation...
    """
    with checkpoints.get_samples.get().output[0].open() as f:
        return rules.run_GRUPS.output

rule archive_GRUPS_results:
    """
    Create a generation-wise archive of GRUPS-rs' result. Results are compressed in lzma format.
    """
    input:
        results_dir = check_GRUPS_checkpoint,
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
    resources:
        cores = lambda w, threads: threads
    log:       "logs/07-archive/run{run}/{generation}/archive_GRUPS_results.log"
    benchmark: "benchmarks/07-archive/run{run}/{generation}/archive_GRUPS_results.tsv"
    threads:   4
    shell: r""" # use of raw string to prevent python escape character SyntaxWarning
        XZ_OPT="-{params.compress_level} -T{threads}" tar -cJvf {output.archive} {params.target_directory} > {log} 2>&1 

        # Store checksums in metadata file yaml-like format. Skip hidden files
        find {params.target_directory} -type f -not -path '*/.*' -exec md5sum {{}} \; > {output.checksum} 2>> {log}
    """


rule archive_TKGWV2_results:
    """
    Create a generation-wise archive of TKGWV2' merged result.
    """
    input:
        results = rules.merge_TKGWV2_results.output.result
    output:
        archive = protected("{archive_dir}/run-{{run}}/{file}".format(
            archive_dir = config['archive']['archive-dir'],
            file        = rules.merge_TKGWV2_results.output.result
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{file}.md5sum".format(
            archive_dir = config['archive']['archive-dir'],
            file        = rules.merge_TKGWV2_results.output.result
        )),
    resources:
        cores = lambda w, threads: threads
    log:     "logs/07-archive/run{run}/{generation}/archive_TKGWV2_results.log"
    threads: 1
    shell: """
        cp {input.results} {output.archive} > {log} 2>&1

        # store checksums in metadata file yaml-like format.
        md5sum {input.results} > {output.checksum} 2>> {log}
    """

rule archive_KIN_results:
    """
    Create a generation-wise archive of KIN's result. Results are compressed in lzma format.
    Note that MD5 checksums are only computed for files located 1-level deep in the results directory.
    """
    input:
        results = rules.run_KIN.output.kin_results,
        #overlap = rules.run_KINgaroo.output.overlap
        #overlap = []
    output:
        archive = protected("{archive_dir}/run-{{run}}/{KIN_dir}.tar.xz".format(
            archive_dir = config['archive']['archive-dir'],
            KIN_dir     = rules.run_KIN.output.kin_results
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{KIN_dir}.md5sums".format(
            archive_dir = config['archive']['archive-dir'],
            KIN_dir     = rules.run_KIN.output.kin_results
        ))
    params:
        compress_level  = config['archive']['compress-level']
    resources:
        cores = lambda w, threads: threads
    log:       "logs/07-archive/run{run}/{generation}/archive_KIN_results.log"
    benchmark: "benchmarks/07-archive/run{run}/{generation}/archive_KIN_results.tsv"
    threads:   4
    shell: r""" # use of raw string to prevent python escape character SyntaxWarning
        XZ_OPT="-{params.compress_level} -T{threads}" tar -cJvf {output.archive} {input.results} > {log} 2>&1 

        # Store checksums in metadata file yaml-like format. Skip hidden files & subdirectories (so many files...)
        find {input.results} -type f -not -path '*/.*' -maxdepth 1 -exec md5sum {{}} \; > {output.checksum} 2>> {log}
    """

rule archive_correctKin_results:
    """
    Create a generation-wise archive of correctKin's result. Results are compressed in lzma format
    """
    input:
        filter_relates = expand(rules.correctKin_filterRelates.output, refpop=config['kinship']['correctKin']['reference-pop'], generation="{generation}"),
        overlap        = expand(rules.correctKin_markerOverlap.output, refpop=config['kinship']['correctKin']['reference-pop'], generation="{generation}")
    output:
        archive = protected("{archive_dir}/run-{{run}}/{correctKin_dir}.tar.xz".format(
            archive_dir    = config['archive']['archive-dir'],
            correctKin_dir = "results/04-kinship/correctKin/{generation}"
        )),
        checksum = protected("{archive_dir}/run-{{run}}/{correctKin_dir}.md5sums".format(
            archive_dir    = config['archive']['archive-dir'],
            correctKin_dir = "results/04-kinship/correctKin/{generation}"
        ))
    params:
        compress_level     = config['archive']['compress-level']
    resources:
        cores = lambda w, threads: threads
    log:       "logs/07-archive/run{run}/{generation}/archive_correctKin_results.log"
    benchmark: "benchmarks/07-archive/run{run}/{generation}/archive_correctKin_results.tsv"
    threads:   4
    shell: """ # use of raw string to prevent python escape character SyntaxWarning
        XZ_OPT="-{params.compress_level} -T{threads}" tar -cJvf {output.archive} {input.filter_relates} {input.overlap} > {log} 2>&1 

        # Store checksums in metadata file yaml-like format. Skip hidden files & subdirectories (so many files...)
        md5sum {input.filter_relates} {input.overlap} > {output.checksum} 2>> {log}
    """


def expand_checksums(
    checksum,
    genlist = get_generations(),
    runlist = "{run}"
):
    """
    Thin wrapper to the expand() function, allows for a bit more concise declaration in rule
    check_duplicate_archives.
    """
    return expand(checksum, generation = genlist, run = runlist)


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
        pipeline_config           = rules.archive_config_file.output.checksum,
        # ---- Pipeline metadata
        pipeline_meta_checksum    = rules.archive_pipeline_metadata.output.checksum,
        # ---- Pedigree simulations
        ped_sim_checksum          = rules.archive_ped_sim_vcf.output.checksum,
        # ---- Gargammel contaminant table
        cont_checksum             = rules.archive_contaminants.output.checksum,
        # ---- MapDamage misincorporation.txt file
        misincorporation_checksum = lambda w: expand(rules.archive_mapdamage_misincorporation.output.checksum,
             sample = get_all_samples_ids(w), run="{run}"
        ),
        # ---- Variant calling callset.
        callset_checksum          = rules.archive_variant_callset.output.checksum,
        # ---- KIN results callset.
        KIN_checksums             = expand_checksums(rules.archive_KIN_results.output.checksum),
        # ---- READ results callset.
        READ_checksums            = expand_checksums(rules.archive_READ_results.output.checksum),
        # ---- READv2 results callset.
        READv2_checksums          = expand_checksums(rules.archive_READv2_results.output.checksum),
        # ---- GRUPS results callset.
        GRUPS_checksums           = expand_checksums(rules.archive_GRUPS_results.output.checksum),
        # ---- TKGWV2 results callset.
        TKGWV2_checksums          = expand_checksums(rules.archive_TKGWV2_results.output.checksum),
        # ---- correctKin results callset.
        correctKin_checksums      = expand_checksums(rules.archive_correctKin_results.output.checksum),
    output:
        metadata = temp(touch(f"{config['archive']['archive-dir']}/run-{{run}}-check_duplicate_archives.done"))
    resources:
        cores = lambda w, threads: threads
    log:     "logs/07-archive/run{run}/check_duplicate_archives.log"
    threads: 1 
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
                "misincorporation": dict(),
                "contaminants"    : dict(),
                "KIN"             : dict(),
                "READ"            : dict(),
                "READv2"          : dict(),
                "GRUPS"           : dict(),
                "TKGWV2"          : dict(),
                "correctKin"      : dict(),
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

        # Extend our global checksum file w/ the misincorporation file of this run
        for checksum in input.misincorporation_checksum:
            extend_checksums(global_checksums, checksum, run_key, "misincorporation")

        # Extend our global checksum file w/ the variant callset of this run.
        extend_checksums(global_checksums, input.callset_checksum, run_key, "variant-callset")

        # Extend our global checksum file w/ each pedigree's KIN result.
        for checksum in input.KIN_checksums:
            extend_checksums(global_checksums, checksum, run_key, "KIN")

        # Extend our global checksum file w/ each pedigree's READ result.
        for checksum in input.READ_checksums:
            extend_checksums(global_checksums, checksum, run_key, "READ")

        # Extend our global checksum file w/ each pedigree's READ result.
        for checksum in input.READv2_checksums:
            extend_checksums(global_checksums, checksum, run_key, "READv2")

        # Extend our global checksum file w/ each pedigree's GRUPS result.
        for checksum in input.GRUPS_checksums:
            extend_checksums(global_checksums, checksum, run_key, "GRUPS")

        # Extend our global checksum file w/ each pedigree's GRUPS result.
        for checksum in input.TKGWV2_checksums:
            extend_checksums(global_checksums, checksum, run_key, "TKGWV2")

        # Extend our global checksum file w/ each pedigree's correctKin result.
        for checksum in input.correctKin_checksums:
            extend_checksums(global_checksums, checksum, run_key, "correctKin")

        with open(input.metadata[0], "w") as yml:
            documents = yaml.dump(global_checksums, yml)


rule archive:
    input: 
        # run number must be defined ONCE, at the final rule, or the function will be executed
        # multiple times.
        # Caching does not seem to resolve the trick.
        # @ TODO: READ & GRUPS Checkpoints is causing rule re-triggering...
        checksums =  lambda w: expand(rules.check_duplicate_archives.output,
            run=define_run_number(w)
        )