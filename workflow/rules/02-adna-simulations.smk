import random
import sys
from os.path import dirname

configfile: "./config/config.yml"

localrules: get_contamination_table

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 00. Simulate modern human contamination for Gargammel.

def get_contaminants(wildcards):
    """
    Returns a list of random ID(s) from a list of potential individuals to use as contamination.
    This will give out a single contaminating individual for each pedigree replicate.
    @ TODO: Maybe adding this info to the global metadata.yml file would be a good idea ?
    """

    # Prevent unneccessary rule re-run trigger events: If the contamination table already exists,
    # simply print out its contents and leave...
    if Path(rules.get_contamination_table.output.cont_table).exists():
        with open(rules.get_contamination_table.output.cont_table, "r") as contamination_table:
            samples = [line.strip("\n").split()[1] for line in contamination_table.readlines()]
            return samples

    # else, open the pedigree generation file and 1000g panel definition file and return a list
    # of random samples.
    # @ TODO: This might be where the FTP OS error is located ?
    #with open(rules.run_ped_sim.input.definition) as f, open(rules.fetch_samples_panel.output.panel) as samples:

    gen_no              = config['ped-sim']['replicates']
    contam_samples_file = rules.get_target_pop_samples.output.target_list.format(
        POP=config["gargammel"]["params"]["contam-pop"]
    )
    with open(contam_samples_file) as samples:
        # get the number of replicates / generations.
        #gen_no = int(list(dropwhile(lambda x: x.startswith('#'), [line for line in f]))[0].split(' ')[2])

        # Fetch sample-IDs matching the user-provided contaminating population tag ('EUR', 'AFR', 'YRI', etc...)
        cont_pop_tag = config["gargammel"]["params"]["contam-pop"]
        #contaminants = [sample.strip("\n").split("\t")[0] for sample in samples.readlines() if cont_pop_tag in sample]
        contaminants = [sample.strip("\n").split("\t")[0] for sample in samples.readlines()]
        # Return a random list of size gen_no (one contaminating individual per replicate).
        return random.sample(contaminants, gen_no)


checkpoint get_contamination_table:
    """
    Call `get_contaminants` and output the results into a tsv file.
    @ FIXME: This might be where my FTP Error is stemmed. If this the case, please take the time to submit an issue.
    """
    input:
        samples_panel = expand(
            rules.get_target_pop_samples.output.target_list,
            POP=config['gargammel']['params']['contam-pop']
        )
    output:
        cont_table    = "results/01-gargammel/contaminants/contaminants.tsv",
    params:
        cont          = get_contaminants
    resources:
        cores = lambda w, threads: threads
    log:      "logs/01-gargammel/get_contamination_table.log"
    conda:    "../envs/coreutils-9.1.yml"
    priority: 99
    threads:  1 
    shell: """
        echo {params.cont} | awk 'BEGIN{{RS=\" \"}}{{print \"ped\"NR, $1}}' > {output.cont_table}
    """


rule create_human_contamination:
    """
    Generate phased haplotypes for contminating individuals using a preprocessed merged 1000g VCF and a provided 
    reference fasta file.

    Sample IDs are generated, stored and fetched once per run from the contamination table. (contaminants.tsv)

    (Copy Number Variants, i.e. 'CN[0-9]+' are unsupported by bcftools, thus they must be excluded.)
    """
    input:
        vcf = lambda wildcards: expand(rules.concat_1000_genomes.output.merged_vcf, 
            POP=config["gargammel"]["params"]["contam-pop"]
        ),
        tbi = lambda wildcards: expand(f"{rules.concat_1000_genomes.output.merged_vcf}.tbi",
            POP=config["gargammel"]["params"]["contam-pop"]
        ),
        chr_ref = lambda wildcards: dirname(ReferenceGenome.get_path()) + "/splitted/{chr}.fasta"
    output:
        hap1 = temp("results/01-gargammel/contaminants/{cont}/{chr}/{cont}_chr{chr}_haplo1.fasta"),
        hap2 = temp("results/01-gargammel/contaminants/{cont}/{chr}/{cont}_chr{chr}_haplo2.fasta")
    params:
        exclude = 'ALT~"<CN[0-9].*>"||ALT~"<INS:.*>" || ALT~"<INV>"'
    resources:
        tmpdir  = config["tempdir"],
        runtime = 10,
        mem_mb  = 128,
        cores   = lambda w, threads: threads
    log: 
        hap1 = "logs/01-gargammel/create_human_contamination/{cont}_chr{chr}_haplo1.log",
        hap2 = "logs/01-gargammel/create_human_contamination/{cont}_chr{chr}_haplo2.log"
    benchmark: "benchmarks/01-gargammel/create_human_contamination/{cont}_chr{chr}.tsv"
    conda:     "../envs/bcftools-1.15.yml"
    priority:  99
    threads:   2
    shell: """
    bcftools consensus -e '{params.exclude}' -H 1 -f {input.chr_ref} --sample {wildcards.cont} {input.vcf} -o {output.hap1} 2> {log.hap1} \
    & \
    bcftools consensus -e '{params.exclude}' -H 2 -f {input.chr_ref} --sample {wildcards.cont} {input.vcf} -o {output.hap2} 2> {log.hap2} \
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 01. Fetch the bacterial contamination database for Gargammel.

rule fetch_bacterial_contamination:
    """
    Placeholder rule. This is the step where we'll download contaminating bacterial fastas
    As of now, stolen from gargammel's Makefile 'makebacterialex' target
    """
    output:
        out_dir        = directory("data/gargammel/bact"),
        fasta_dir      = directory("data/gargammel/bact/bact"),
        abundance_list = "data/gargammel/bact/bact/list",
        taxa_list      = "data/gargammel/bact/all_taxa.tsv"
    params:
        url = "https://www.dropbox.com/s/obmr48d72ahjvhp/clovis.tar.gz?dl=1"
        #url = "https://www.dropbox.com/s/1pdbqbguw0jfzib/k14.tar.gz?dl=1"
    resources:
        cores = lambda w, threads: threads
    log:     "logs/01-gargammel/fetch_bacterial_contamination.log"
    conda:   "../envs/coreutils-9.1.yml"
    threads: 1
    shell: """
        root_dir=$(pwd)
        cd {output.out_dir} > {log} 2>&1 

        # Download bacterial data
        wget --quiet -O- {params.url} | tar xvzf - --strip-components=1 >> $root_dir/{log} 2>&1 

        # Rename "fasta" directory to "bact"
        mv fasta/* bact/ >> $root_dir/{log} 2>&1 
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 02. Generate ancient DNA consensus sequences for our simulated samples.

rule get_consensus:
    """
    Generate phased haplotypes for pedigree individuals using ped-sim's output VCF and a provided reference fasta file.
    (Copy Number Variants, i.e. 'CN[0-9]+' are unsupported by bcftools, thus they must be excluded.)
    """
    input:
        vcf     = expand(rules.merge_ped_sim.output.vcf, POP=config["ped-sim"]["params"]["pop"]),
        tbi     = expand(rules.merge_ped_sim.output.vcf + ".tbi", POP=config["ped-sim"]["params"]["pop"]),
        chr_ref = lambda wildcards: dirname(ReferenceGenome.get_path()) + "/splitted/{chr}.fasta"
    output:
        hap1 = temp("results/01-gargammel/{sample}/{chr}/endo/{sample}_chr{chr}_haplo1.fasta"),
        hap2 = temp("results/01-gargammel/{sample}/{chr}/endo/{sample}_chr{chr}_haplo2.fasta"),
    params:
        exclude = 'ALT~"<CN[0-9].*>"||ALT~"<INS:.*>" || ALT~"<INV>"'
    resources:
        tmpdir  = config["tempdir"],
        runtime = 10,
        mem_mb  = 128,
        cores   = lambda w, threads: threads
    log: 
        hap1 = "logs/01-gargammel/get_consensus/{sample}-chr{chr}_haplo1.log",
        hap2 = "logs/01-gargammel/get_consensus/{sample}-chr{chr}_haplo2.log"
    benchmark: "benchmarks/01-gargammel/get_consensus/{sample}-chr{chr}_haplo1.tsv",
    conda:     "../envs/bcftools-1.15.yml"
    group: "scatter"
    threads: 2
    shell: """
        bcftools consensus -e '{params.exclude}' -H 1 -f {input.chr_ref} --sample {wildcards.sample} {input.vcf} -o {output.hap1} 2> {log.hap1} \
        & \
        bcftools consensus -e '{params.exclude}' -H 2 -f {input.chr_ref} --sample {wildcards.sample} {input.vcf} -o {output.hap2} 2> {log.hap2} \
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 03. Choose between DWGSIM or Gargammel to simulate reads.

def choose_simulator(wildcards, orientation):
    """
    Read the config file and resolve the requested FASTQ simulator (gargammel, dwgsim).
    """
    allowed_orientations = ['forward', 'reverse']
    if orientation not in allowed_orientations:
        raise RuntimeError("Invalid orientation: {orientation}. Valid values are {allowed_orientations}")

    match config["simulator"]:
        case "gargammel":
            fastq = rules.run_gargammel.output.forwd if orientation == "forward" else rules.run_gargammel.output.revrs,
        case "dwgsim":
            fastq = rules.run_dwgsim.output.forwd if orientation == "forward" else rules.run_dwgsim.output.revrs,

    return expand(fastq, chr=range(1,23), sample="{sample}")

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 04-A Run Gargammel to simulate raw *ancient DNA* FASTQ files from our simulated samples.

def find_contaminant(wildcards):
    """
    Assign the correct contaminating individual to a given pedigree sample, according to its generation #.
    """
    # Get the current pedigree sampleID, chromosome and generation.
    sample = wildcards.sample
    chromo = wildcards.chr
    gen    = sample.split("_")[0] 
    
    # ---- Do not trigger the checkpoint if the corresponding file already exists...
    if Path(rules.get_contamination_table.output.cont_table).exists():
        contamination_table = rules.get_contamination_table.output.cont_table
    else:
        contamination_table = checkpoints.get_contamination_table.get().output.cont_table
    marked = False
    with open(contamination_table) as f:
        for line in f.readlines():
            if line.startswith(f"{gen} "):
                marked = True
                contaminant = line.strip("\n").split(" ")[1]
                break
    
    # TODO: Refactor this.
    if not marked:
        raise RuntimeError("Failed to find {gen} pattern in {contamination_table}")

    out = "results/01-gargammel/contaminants/{cont}/{{chr}}/{cont}_chr{{chr}}_haplo{haplo}.fasta"
    if float(config['gargammel']['comp-cont']) > 0.0:
        return expand(out, cont=contaminant, chr=chromo, haplo=[1,2])
    else:
        return []


def get_pmd_model(wildcards):
    """
    Choose between the Briggs model, or a raw MapDamage misincorporation file 
    to simulate PMD.
    """
    pmd_model = config['gargammel']['pmd-model']
    match pmd_model:
        case "misincorporation":
            mis_file = config['gargammel']['misincorporation']['file']
            protocol = config['gargammel']['misincorporation']['protocol']
            if protocol not in ["double", "single"]:
                raise RuntimeError(f"Invalid Gargammel misincorporation protocol: '{protocol}'")
            return f"-mapdamage {mis_file} double"

        case "briggs":
            v = config['gargammel']['briggs']['nick-frequency']
            l = config['gargammel']['briggs']['overhang-length']
            d = config['gargammel']['briggs']['ds-deaminations']
            s = config['gargammel']['briggs']['ss-deaminations']
            return f"-damage {v},{l},{d},{s}"
        case other:
           raise RuntimeError(f"Invalid Gargammel pmd-model value: '{pmd_model}'")


def parse_gargammel_comp(wildcards):
    cont = config['gargammel']['comp-cont']
    endo = config['gargammel']['comp-endo']
    bact = config['gargammel']['comp-bact']


    if config['gargammel']['params']['contaminate-samples'] is None:
        pass
    else: 
        try:
            with open(config['gargammel']['params']['contaminate-samples']) as f:
                candidates = [line.strip('\n') for line in f.readlines()]
            if wildcards.sample.split("_")[1] not in candidates:
                (endo, cont) = (1.0, 0.0)
        except Exception as e:
            raise e
    return f"{bact},{cont},{endo}"


rule run_gargammel:
    """
    Perform raw ancient DNA sequencing data using Gargammel.
    """
    input:
        hap1                    = rules.get_consensus.output.hap1,
        hap2                    = rules.get_consensus.output.hap2,
        cont_table              = rules.get_contamination_table.output.cont_table,
        human_contamination     = find_contaminant,
        bacterial_contamination = rules.fetch_bacterial_contamination.output.fasta_dir,
    output:
        forwd = "results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}_s1.fq.gz",
        revrs = "results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}_s2.fq.gz",
        a     = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}_a.fa.gz"),
        b     = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}.b.fa.gz"),
        c     = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}.c.fa.gz"),
        d     = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}_d.fa.gz"),
        e     = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}.e.fa.gz")
    params:
        comp               = parse_gargammel_comp,
        coverage           = config['gargammel']['coverage'],
        size_freq          = config['gargammel']['sizefreq'],
        qshift             = config['gargammel']['qshift'],
        misincorporation   = get_pmd_model,
        output_base_name   = "results/01-gargammel/{sample}/{chr}",
        input_directory    = directory("results/01-gargammel/{sample}/{chr}"),
    resources:
        tmpdir  = config["tempdir"],
        runtime = 10,
        mem_mb  = 128, 
        cores   = lambda w, threads: threads
    log:       "logs/01-gargammel/run_gargammel/{sample}_chr{chr}.log"
    benchmark: "benchmarks/01-gargammel/run_gargammel/{sample}_chr{chr}.tsv"
    conda:     "../envs/gargammel-1.1.4.yml"
    group:     "scatter"
    priority:  2
    threads:   1
    shell: """
        mkdir -p {params.output_base_name}/cont                             >  {log} 2>&1
        ln -sfrt {params.output_base_name} {input.bacterial_contamination}  >> {log} 2>&1
        if [ "{input.human_contamination}" ]; then
            ln -sfrt {params.output_base_name}/cont {input.human_contamination} >> {log} 2>&1
        fi

        gargammel \
        --comp {params.comp} \
        {params.misincorporation} \
        -c {params.coverage} \
        -f {params.size_freq} \
        -o {params.output_base_name}/{wildcards.sample}_chr{wildcards.chr} \
        -qs {params.qshift} -qs2 {params.qshift} \
        {params.input_directory} >> {log} 2>&1
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 04-B Run DWGSIM to simulate raw *modern DNA* FASTQ files from our simulated samples.

def set_dwgsim_seed(wildcards):
    seed = config['dwgsim']['seed']
    if seed is None:
        with open(rules.meta.output.metadata) as f:
            metadata = yaml.load(f, Loader=yaml.loader.SafeLoader)
            seed     = metadata['seed']
    return seed

rule run_dwgsim:
    """
    simulate raw, modern sequencing data, using dwgsim.
    @TODO: Incomplete implementation - downstream BAM preprocessing is absolutely broken when using this software (picard Rmdup especially)
    """
    input:
        metadata                = "results/meta/pipeline-metadata.yml",
        hap1                    = rules.get_consensus.output.hap1,
        hap2                    = rules.get_consensus.output.hap2,
    output:
        forwd  = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}_s1.dwgsim.fq.gz"),
        revrs  = temp("results/01-gargammel/{sample}/{chr}/{sample}_chr{chr}_s2.dwgsim.fq.gz")
    params:
        coverage   = config['gargammel']['coverage'],
        output_dir = "results/01-gargammel/{sample}/{chr}/",
        seed       = set_dwgsim_seed
    resources:
        cores = lambda w, threads: threads
    log:       "logs/01-gargammel/run_dwgsim/{sample}_chr{chr}.log"
    benchmark: "benchmarks/01-gargammel/run_dwgsim/{sample}_chr{chr}.log"
    conda: "../envs/dwgsim-1.1.13.yml"
    threads: 2
    shell: """
        BASENAME="{params.output_dir}/{wildcards.sample}_chr{wildcards.chr}"
        dwgsim -H -C {params.coverage} -z {params.seed} -r0 -R0 -X0 -F 0 -y 0 -n0 -c0 -S0 -q~ {input.hap1} ${{BASENAME}}_haplo1 >  {log} 2>&1
        dwgsim -H -C {params.coverage} -z {params.seed} -r0 -R0 -X0 -F 0 -y 0 -n0 -c0 -S0 -q~ {input.hap2} ${{BASENAME}}_haplo2 >> {log} 2>&1
        cat ${{BASENAME}}_haplo1.bwa.read1.fastq.gz ${{BASENAME}}_haplo2.bwa.read1.fastq.gz > ${{BASENAME}}_s1.dwgsim.fq.gz 2>> {log}
        cat ${{BASENAME}}_haplo1.bwa.read2.fastq.gz ${{BASENAME}}_haplo2.bwa.read2.fastq.gz > ${{BASENAME}}_s2.dwgsim.fq.gz 2>> {log}
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 05. Merge chromosomes for each samples

rule merge_chromosomes:
    """
    Gather all simulated the chromosomes-splitted FASTQ simulation data of a given sample, and merge them.
    """
    input:
        forwd = lambda wildcards: choose_simulator(wildcards, "forward"),
        revrs = lambda wildcards: choose_simulator(wildcards, "reverse"),
    output:
        forwd = temp("results/02-preprocess/00-raw/{sample}_s1.fq.gz"),
        revrs = temp("results/02-preprocess/00-raw/{sample}_s2.fq.gz")
    log: 
        forwd = "logs/01-gargammel/merge_chromosomes/{sample}_s1.log",
        revrs = "logs/01-gargammel/merge_chromosomes/{sample}_s2.log"
    resources:
        tmpdir  = config["tempdir"],
        runtime = 10,
        mem_mb  = 128, 
        cores   = lambda w, threads: threads
    benchmark: "benchmarks/01-gargammel/merge_chromosomes/{sample}.tsv"
    conda:     "../envs/coreutils-9.1.yml"
    group:     "scatter"
    priority:  3
    threads:   2
    shell: """
        zcat {input.forwd} | gzip > {output.forwd} 2> {log.forwd}
        zcat {input.revrs} | gzip > {output.revrs} 2> {log.revrs}
    """

ruleorder: multiqc_adna_simulations > multiqc
use rule multiqc as multiqc_adna_simulations with:
    input:
        required_files = lambda w: temp(
            expand(rules.run_fastqc_fq.output.data,
                directory="results/02-preprocess/00-raw",
                file=[f"{sample}_{s}" for sample in get_all_samples_ids(w) for s in ("s1", "s2")]
            )
        ),
    output: 
        html = report("results/02-preprocess/00-raw/multiqc-report.html",
            caption     = "../report/02-adna-simulations/multiqc-raw.rst",
            category    = "01. Ancient DNA simulations",
        )
    params:
        extra_args = rules.multiqc.params.extra_args,
        extra_dirs = directory(expand("results/02-preprocess/{subdir}", subdir=["00-raw"]))
