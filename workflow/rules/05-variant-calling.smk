localrules: generate_bam_list

# ------------------------------------------------------------------------------------------------------------------- #
# ---- Generic / Utility rules and functions.

def get_samples_ids(wildcards):
    """
    Get each pedigree sample's id based on each unique pedigree sample id + the
    number of pedigree replicates. This list is sorted to match UNIX's 
    """
    # Run through the initial samples files and extract pedigree ids 
    with checkpoints.get_samples.get().output[0].open() as f:
        samples = str.split(f.readline().replace('\n', ''), '\t')
        ids     = set([sample.split('_')[1] for sample in samples])
        return sorted(expand("{generation}_{ids}", ids=ids, generation="{generation}"), key=str.casefold)


rule eigenstrat_to_UCSC_BED:
    """
    Convert an eigensoft's .snp file to a generic bed file 
    (termed 'ucscbed', to avoid collisions with plink's file format).
    """
    input:
        snp   = "{directory}/{eigenstrat}.snp"
    output:
        bed   = "{directory}/{eigenstrat}.ucscbed"
    resources:
        cores = lambda w, threads: threads
    log:       "logs/generics/{directory}/eigenstrat_to_UCSC_BED-{eigenstrat}.log"
    benchmark: "benchmarks/generics/{directory}/eigenstrat_to_UCSC_BED-{eigenstrat}.tsv"
    conda:     "../envs/coreutils-9.1.yml"
    threads:   1
    shell: """
        awk 'BEGIN{{OFS="\t"}}{{print $2, $4-1, $4, $5, $6}}' {input.snp} > {output.bed}
    """


rule plink_bfile_to_tped:
    """
    Convert a binarized PLINK fileset into a human readable transposed set.
    """
    input:
        bfile    = multiext("{directory}/{file}", ".bed", ".bim", ".fam"),
    output:
        tplink   = multiext("{directory}/{file}", ".tped", ".tfam") 
    params:
        basename = "{directory}/{file}"
    resources:
        cores    = lambda w, threads: threads
    log:       "logs/generics/{directory}/plink_bfile_to_tped-{file}.log"
    benchmark: "benchmarks/generics/{directory}/plink_bfile_to_tped-{file}.log"
    conda:     "../envs/plink-1.9.yml"
    threads:   1
    shell: """
        plink \
        --bfile {params.basename} \
        --out {params.basename} \
        --recode transpose tab \
        --allow-no-sex \
        --keep-allele-order > {log} 2>&1
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 00. Download out variant callset.

#module netrules:
#    snakefile: "00-netrules.smk"
#    config: config
#
#use rule download_reich_1240K from netrules

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 01. Generate a list of input BAM list for cohort variant calling (one list per pedigree).

def get_pileup_input_bams(wildcards):
    """
    Define the appropriate input bam for the variant caller, based on which 
    PMD-rescaling method was requested by the user.
    """
    # Run through the initial samples files and extract pedigree ids 
    samples_ids = get_samples_ids(wildcards)

    # Return a list of input bam files for pileup
    rescaler = config['preprocess']['pmd-rescaling']['rescaler']
    match rescaler:
        case "mapdamage":
            print("NOTE: Using MapDamage rescaled bams for variant calling.", file=sys.stderr)
            return expand(
                "results/02-preprocess/06-mapdamage/{sample}/{sample}.srt.rmdup.rescaled.bam",
                sample = samples_ids,
            )
        case "pmdtools":
            print("NOTE: Using PMDTools rescaled bams for variant calling.", file=sys.stderr)
            return expand(
                "results/02-preprocess/06-pmdtools/{sample}/{sample}.srt.rmdup.filtercontam.bam",
                sample = samples_ids,

            )
        case None:
            print("WARNING: Skipping PMD Rescaling for variant calling!", file=sys.stderr)
            return expand(
                define_dedup_input_bam(wildcards),
                sample = samples_ids
            )
    raise RuntimeError(f'Invalid rescaler value "{rescaler}')


rule generate_bam_list:
    """
    Create an ordered list of bam files using the original samples file
    """
    input:
        bamlist = get_pileup_input_bams
    output:
        bamlist = "results/03-variant-calling/01-pileup/{generation}/{generation}.bam.list"
    resources:
        cores   = lambda w, threads: threads
    log:      "logs/03-variant-calling/generate_bam_list/{generation}.log"
    conda:    "../envs/coreutils-9.1.yml"
    threads:  1
    priority: 15
    shell: """
        ls {input.bamlist} > {output.bamlist} 2> {log}
    """


# ------------------------------------------------------------------------------------------------------------------- #
# ---- 02. Generate pileup files using samtools. One pileup per pedigree.

def get_samtools_optargs(wildcards):
    """
    Inject the appropriate flag to disable samtools' Base alignment quality,
    if the user requested it.
    """
    optargs = ""
    if config['variant-calling']['pileup']['disable-BAQ']:
        optargs += "-B "
    return optargs


rule get_target_panel_intersect:
    """
    Get the intersect between our pedigree simulations variant set (A) and our target variant set (B)
    - This is required to prevent bias arising from targeting SNPs that were not simulated. In other terms,
      targeting the B set when n(A∩B) < n(B) would cause a drastic decrease in observed average heterozygocity.

    - This rule is made to ensure n(B) == n(A∩B), at the cost of loosing target coverage.

    MAF support: wildcards 'maf' and 'superpop' respectively represent the requested minor allele frequency and
    the **super**-population on which this frequency is computed.
    - Uses the 1000g {POP}_AF annotations, which are kept in ped-sim's output vcf.
    - Mainly intended for lcMLkin (?).

    # Benchmarks: 
    | depth | max_rss |
    | ----- | ------- |
    | 0.01X |         |
    | 0.05X | 51.91   |
    | 0.10X | 37.91   |
    """
    input:
        ped_vcf = multiext(rules.dopplegang_twins.output.merged_vcf.format(POP=config['ped-sim']['params']['pop']), "", ".tbi"),
        targets = os.path.splitext(config["kinship"]["targets"])[0] + ".ucscbed",
    output:
        targets = "results/03-variant-calling/00-panel/variants-intersect-{superpop}_maf{maf}.ucscbed"
    params:
        exclude = lambda w: f"{w.superpop}_AF<{w.maf} || {w.superpop}_AF>{1-float(w.maf)}"
    resources:
        runtime = 60,
        mem_mb  = 128,
        cores   = lambda w, threads: threads
    log:       "logs/03-variant-calling/00-panel/get_target_panel_intersect-{superpop}-maf{maf}.log"
    benchmark: "benchmarks/03-variant-calling/00-panel/get_target_panel_intersect-{superpop}-maf{maf}.tsv"
    conda:     "../envs/bcftools-1.15.yml"
    threads:   4
    shell: """
        bcftools view --threads {threads} -H -v snps -R {input.targets} -e '{params.exclude}' {input.ped_vcf} \
        | awk 'BEGIN{{OFS="\t"}}{{print $1, $2-1, $2, $4, $5}}' \
        > {output.targets} 2> {log}
    """

rule samtools_pileup:
    """
    Generate a legacy pileup file using samtools mpileup.

    # Benchmarks: 
    | depth | max h:m:s | max_rss |
    | ----- | --------- | ------- |
    | 0.01X |           |         |
    | 0.05X | 0:00:42   | 404     |
    | 0.10X |           | 398     |

    """
    input:
        bamlist   = rules.generate_bam_list.output.bamlist,
        targets   = expand(rules.get_target_panel_intersect.output.targets, 
            maf      = config['variant-calling']['maf'],
            superpop = config['variant-calling']['maf-superpop'],
        ),
        reference = config["reference"],
        fai       = config["reference"] + ".fai",
    output:
        pileup    = "results/03-variant-calling/01-pileup/{generation}/{generation}.pileup"
    params:
        min_MQ    = config['variant-calling']['pileup']['min-MQ'],
        min_BQ    = config['variant-calling']['pileup']['min-BQ'],
        optargs   = get_samtools_optargs
    resources:
        runtime = 10,
        mem_mb  = 512,
        cores   = lambda w, threads: threads
    log:       "logs/03-variant-calling/samtools_pileup/{generation}.log"
    benchmark: "benchmarks/03-variant-calling/samtools_pileup/{generation}.tsv"
    conda:     "../envs/samtools-1.15.yml"
    threads:   1
    shell: """
        (samtools mpileup \
        -R {params.optargs} \
        -q {params.min_MQ} \
        -Q {params.min_BQ} \
        -l {input.targets} \
        -f {input.reference} \
        -b {input.bamlist} \
        | LC_ALL=C sort -k1,1n -k2,2n \
        )> {output.pileup} 2> {log}
        """


# ------------------------------------------------------------------------------------------------------------------- #
# ---- 03-A. Perform pseudo-haploid random variant calling with SequenceTools' PileupCaller

def parse_pileup_caller_flags(wildcards):
    args=""
    if config['variant-calling']['pileupCaller']["skip-transitions"]:
        args+="--skipTransitions "

    seed = config['variant-calling']['pileupCaller']['seed']
    if seed is None:
        with open(rules.meta.output.metadata) as f:
            metadata = yaml.load(f, Loader=yaml.loader.SafeLoader)
            seed     = metadata['seed']

    args += f"--seed {seed} "

    match config['variant-calling']['pileupCaller']["mode"]:
        case "randomHaploid":
            args += "--randomHaploid "
        case "majorityCall":
            args += "--majorityCall "
        case _:
            raise RuntimeError("Incorrect pileupCaller mode selected.")
    return args


rule pileup_caller:
    """
    Perform random sampling variant calling uing SequenceTools' pileupCaller

    # Benchmarks: 
    | depth | max h:m:s | max_rss |
    | ----- | --------- | ------- |
    | 0.01X |           |         |
    | 0.05X | 0:00:39   | 17.64   |
    | 0.10X |           | 30.15   |
    """
    input:
        samples_def       = rules.get_samples.output, 
        pileup            = rules.samtools_pileup.output.pileup,
        bamlist           = rules.generate_bam_list.output.bamlist,
        targets           = os.path.splitext(config["kinship"]["targets"])[0] + ".snp",
        metadata          = "results/meta/pipeline-metadata.yml"
    output:
        plink             = multiext("results/03-variant-calling/02-pileupCaller/{generation}/{generation}", ".bed", ".bim", ".fam"),
        sample_names_file = "results/03-variant-calling/02-pileupCaller/{generation}/{generation}-names.txt"
    params:
        basename          = lambda wildcards, output: os.path.splitext(output.plink[0])[0],
        optargs           = parse_pileup_caller_flags,
        min_depth         = config['variant-calling']["pileupCaller"]["min-depth"],
        seed              = config['variant-calling']['pileupCaller']["seed"],
        sample_names      = lambda wildcards: expand(get_samples_ids(wildcards), generation = wildcards.generation)
    resources:
        runtime = 10,
        mem_mb  = 128,
        cores   = lambda w, threads: threads
    log:       "logs/03-variant-calling/pileup_caller/{generation}.log"
    benchmark: "benchmarks/03-variant-calling/pileup_caller/{generation}.tsv"
    conda:     "../envs/sequencetools-1.5.2.yml"
    threads:   1
    shell: """
        echo {params.sample_names} | tr ' ' '\n' > {output.sample_names_file}

        pileupCaller \
        {params.optargs} \
        --snpFile <( awk '($2<=22)' {input.targets}) \
        --sampleNameFile {output.sample_names_file} \
        --plinkOut {params.basename} \
        --minDepth {params.min_depth} \
        --samplePopName {wildcards.generation} < {input.pileup} > {log} 2>&1 
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 03-B. Perform pseudo-haploid random variant calling with ANGSD

rule ANGSD_random_haploid:
    """
    Perform random pseudo-haploid variant calling using ANGSD
    """
    input:
        pileup    = rules.samtools_pileup.output.pileup,
        bamlist   = rules.generate_bam_list.output.bamlist,
        reference = config['reference'],
        fai       = config["reference"] + ".fai",
    output:
        haplos    = "results/03-variant-calling/02-ANGSD/{generation}/{generation}.haplo.gz"
    params:
        out       = lambda wildcards, output: os.path.splitext(output.haplos)[0]
    resources:
        cores     = lambda w, threads: threads
    log:       "logs/03-variant-calling/ANGSD_random_haploid/{generation}.log"
    benchmark: "benchmarks/03-variant-calling/ANGSD_random_haploid/{generation}.tsv"
    threads:   4
    priority:  15
    conda: "../envs/angsd-0.939.yml"
    shell: """
        angsd -pileup {input.pileup} \
        -fai {input.fai} \
        -nInd $(cat {input.bamlist} | wc -l) \
        -rf <(seq 1 22) \
        -dohaplocall 1 \
        -doCounts 1 \
        -nthreads {threads} \
        -out {params.out} 2> {log}
    """


rule ANGSD_haplo_to_plink:
    """
    Convert a .haplo.gz ANGSD file to a set of PLINK .tped / .tfam files 
    """
    input:
        haplos  = rules.ANGSD_random_haploid.output.haplos,
        bamlist = rules.generate_bam_list.output.bamlist 
    output:
        tped    = "results/03-variant-calling/02-ANGSD/{generation}/{generation}.tped",
        tfam    = "results/03-variant-calling/02-ANGSD/{generation}/{generation}.tfam",
    params:
        out     = lambda wildcards, output: os.path.splitext(output.tped)[0]
    resources:
        cores   = lambda w, threads: threads
    log:       "logs/03-variant-calling/ANGSD_haplo_to_plink/{generation}.log"
    benchmark: "benchmarks/03-variant-calling/ANGSD_haplo_to_plink/{generation}.tsv"
    conda:     "../envs/angsd-0.939.yml"
    threads:   1
    priority:  15
    shell: """
        haploToPlink {input.haplos} {params.out} 2>  {log}
        sed -i 's/N/0/g' {output.tped}           2>> {log}
        cat {input.bamlist} \
        | xargs basename -a \
        | grep -oP '^[^.]+(?=(\.[^.]+)*(\.bam$))' \
        | awk 'BEGIN{{OFS="\t"}}{{print "{wildcards.generation}", $1, 0, 0, 0, 0}}' \
        > {output.tfam} 2>> {log}
    """
