# ---- Utility functions
def correctKin_output(wildcards):
    """
    Returns the list of required generation_wise comparisons from correctKin
    @ TODO: This checkpoint is not needed anymore
    """
    # Get the number of expected generations.
    refpop   = config['kinship']['correctKin']['reference-pop']
    template = "results/04-kinship/correctKin/{generation}/{generation}-merge-{refpop}-subset.rels.tsv"
    return expand(template, generation = get_generations(), refpop = refpop)

# ------------------------------------------------------------------------------------------------ #

rule format_convertf_parfile:
    input:
        template  = "resources/convertf/par.EIGENSTRAT.PED"
    output:
        parfile = "{directory}/{eigenstrat}.par.EIGENSTRAT.PED"
    template_engine: "jinja2"

         
rule eigenstrat_to_binary_plink:
    """
    Convert an eigensoft' [.snp|.geno|.ind] fileset to a generic set of [.bed|.bim|.fam] fileset
    """
    input:
        eigensoft = multiext("{directory}/{eigenstrat}", ".geno", ".snp", ".ind"),
        parfile   = rules.format_convertf_parfile.output.parfile
    output:
        plink     = multiext("{directory}/{eigenstrat}", ".bed", ".bim", ".fam")
    resources:
        cores = lambda w, threads: threads
    log:       "logs/generics/{directory}/eigenstrat_to_binary_plink-{eigenstrat}.log"
    benchmark: "benchmarks/generics/{directory}/eigenstrat_to_binary_plink-{eigenstrat}.tsv"
    conda:     "../envs/eigensoft-8.0.0.yml"
    threads: 1
    shell: """
        convertf -p {input.parfile} > {log} 2>&1
    """

rule subset_reich_dataset:
    input:
        bfile = multiext("data/Reich-dataset/1240K/v52.2_1240K_public", ".bed", ".bim", ".fam"),
        panel = "data/vcf/1000g-phase3/samples-list/integrated_call_samples_v3.20130502.ALL.panel"
    output:
        subset = multiext("data/Reich-dataset/1240K/subsets/{refpop}/v52.2_1240K_public", ".bed", ".bim", ".fam"),
    params:
        bfile_basename  = lambda wildcards, input:  splitext(input.bfile[0])[0],
        subset_basename = lambda wildcards, output: splitext(output[0])[0]
    log:       "logs/04-kinship/correctKin/subset_reich_dataset-{refpop}.log"
    benchmark: "benchmarks/04-kinship/correctKin/subset_reich_dataset-{refpop}.tsv"
    conda:     "../envs/plink-1.9.yml"
    threads:   1
    shell: """
        plink \
        --bfile {params.bfile_basename} \
        --keep <(grep -f <(grep "{wildcards.refpop}" {input.panel} | cut -f1) {params.bfile_basename}.fam | awk '{{print $1, $2}}') \
        --allow-no-sex \
        --keep-allele-order \
        --make-bed \
        --threads {threads} \
        --out {params.subset_basename} > {log} 2>&1
    """

rule plink_coverage:
    input:
        plink = multiext("{directory}/{basename}", ".bed", ".bim", ".fam")
    output:
        imiss     = "{directory}/{basename}.imiss",
        lmiss     = "{directory}/{basename}.lmiss",
        genocount = "{directory}/{basename}.genocount"
    log:       "logs/04-kinship/correctKin/plink_coverage/{directory}/{basename}.log"
    benchmark: "benchmarks/04-kinship/correctKin/plink_coverage/{directory}/{basename}.tsv"
    threads:    1
    conda:     "../envs/plink-1.9.yml"
    shell: """
        plink \
        --bfile {wildcards.directory}/{wildcards.basename} \
        --missing \
        --threads {threads} \
        --out {wildcards.directory}/{wildcards.basename} > {log} 2>&1
        awk '(NR>1){{print $5 - $4}}' {output.imiss} > {output.genocount} 2>> {log}
    """


rule correctKin_deplete_inds:
    input:
        subset   = rules.subset_reich_dataset.output.subset,
        coverage = expand(rules.plink_coverage.output.genocount,
            directory = dirname(rules.pileup_caller.output.plink[0]),
            basename  = splitext(basename(rules.pileup_caller.output.plink[0]))[0]
        )
    output:
        depleted = multiext("results/04-kinship/correctKin/{generation}/{refpop}-v52.2_1240K_public-depleted", ".bed", ".bim", ".fam")
    params:
        subset_basename   = lambda wildcards, input: splitext(input.subset[0])[0],
        depleted_basename = lambda wildcards, output: splitext(output.depleted[0])[0]
    log:       "logs/04-kinship/correctKin/correctKin_deplete_inds/{generation}-{refpop}.log"
    benchmark: "benchmarks/04-kinship/correctKin/correctKin_deplete_inds/{generation}-{refpop}.tsv"
    conda:     "../envs/correctKin.yml"
    threads:   1
    shell: """
        depleteIndiv --out {params.depleted_basename} {input.subset[0]} $(sort -n {input.coverage} | head -n1) $(sort -n {input.coverage} | tail -n1) > {log} 2>&1
        for ext in bim fam; do
            cp {params.subset_basename}.$ext {params.depleted_basename}.$ext
        done >> {log} 2>&1
    """


def parse_merge_plink_reich_input(wildcards):
    should_deplete = config['kinship']['correctKin']['deplete-indivs']
    return rules.correctKin_deplete_inds.output.depleted if should_deplete else rules.subset_reich_dataset.output.subset

rule merge_plink_reich_pop:
    input:
        bfile  = rules.pileup_caller.output.plink,
        subset = parse_merge_plink_reich_input,
        ids    = expand(rules.merge_ped_sim.output.ids, POP=config['ped-sim']['params']['pop'])
        #ids    = expand(rules.run_ped_sim.output.ids, POP=config['ped-sim']['params']['pop'])
        #ids    = f"results/00-ped-sim/{config['ped-sim']['params']['pop']}-pedigrees.ids"
    output:
        merged = multiext("results/04-kinship/correctKin/{generation}/{generation}-merge-{refpop}-subset", ".bed", ".bim", ".fam")
    params:
        bfile_basename  = lambda wildcards, input: splitext(input.bfile[0])[0],
        subset_basename = lambda wildcards, input: splitext(input.subset[0])[0],
        merged_basename = lambda wildcards, output: splitext(output.merged[0])[0],
        merge_mode      = 1
    log:       "logs/04-kinship/correctKin/merge_plink_reich_pop/{generation}-{refpop}.log"
    benchmark: "benchmarks/04-kinship/correctKin/merge_plink_reich_pop/{generation}-{refpop}.tsv"
    conda:     "../envs/plink-1.9.yml"
    threads:   1
    shell: """
        plink \
        --bfile {params.bfile_basename} \
        --bmerge {params.subset_basename} \
        --merge-mode {params.merge_mode} \
        --remove <(grep -f <(grep "^{wildcards.generation}" {input.ids} | cut -f2) {params.subset_basename}.fam) \
        --allow-no-sex \
        --keep-allele-order \
        --make-bed \
        --threads {threads} \
        --out {params.merged_basename} > {log} 2>&1
    """

rule run_pcangsd:
    input:
        bfile   = rules.merge_plink_reich_pop.output.merged
    output:
        kinship = multiext(splitext(rules.merge_plink_reich_pop.output.merged[0])[0], ".cov", ".kinship.npy", ".inbreed.npy")
        #kinship = multiext("results/04-kinship/correctKin/{generation}/{generation}-merge-{refpop}-subset", ".cov", ".kinship.npy", ".inbreed.npy")
    params:
        input_basename  = lambda wildcards, input: splitext(input.bfile[0])[0],
        output_basename = lambda wildcards, output: splitext(output.kinship[0])[0] # //Input basename is the same as the output
    log:       "logs/04-kinship/correctKin/run_pcangsd/{generation}-{refpop}.log"
    benchmark: "benchmarks/04-kinship/correctKin/run_pcangsd/{generation}-{refpop}.tsv"
    conda:     "../envs/pcangsd-0.99.yml"
    priority:  75
    threads:   8
    shell: """
        pcangsd -plink {params.input_basename} -o {params.output_basename} -inbreed 1 -kinship -threads {threads} > {log} 2>&1
    """

rule correctKin_markerOverlap:
    input:
        bed     = rules.merge_plink_reich_pop.output.merged[0]
    output:
        overlap = "results/04-kinship/correctKin/{generation}/{generation}-merge-{refpop}-subset.overlap"
    log:       "logs/04-kinship/correctKin/correctKin_markerOverlap/{generation}-{refpop}.log"
    benchmark: "benchmarks/04-kinship/correctKin/correctKin_markerOverlap/{generation}-{refpop}.tsv"
    conda:     "../envs/correctKin.yml"
    threads:   1
    shell: """
        markerOverlap {input.bed} > {log} 2>&1
    """

rule correctKin_filterRelates:
    input:
        overlap = rules.correctKin_markerOverlap.output.overlap,
        kinship = rules.run_pcangsd.output.kinship[1]
    output:
        rels  = "results/04-kinship/correctKin/{generation}/{generation}-merge-{refpop}-subset.rels.tsv",
        corr  = "results/04-kinship/correctKin/{generation}/{generation}-merge-{refpop}-subset.corr.tsv",
        stats = "results/04-kinship/correctKin/{generation}/{generation}-merge-{refpop}-subset.stats.tsv"
    log:       "logs/04-kinship/correctKin/correctKin_filterRelates/{generation}-{refpop}.log"
    benchmark: "benchmarks/04-kinship/correctKin/correctKin_filterRelates/{generation}-{refpop}.tsv"
    conda:     "../envs/correctKin.yml"
    threads:   1
    shell: """
        filterRelates {input.overlap} {input.kinship} > {output.rels} 2> {log}
    """
