
#from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
#from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
#FTP  = FTPRemoteProvider(retry=config['FTP']['retries']) # Anonymous 
#HTTP = HTTPRemoteProvider() # Anonymous 

from os.path import dirname

configfile: "config/config.yml"

# ---- storage definition
storage g1k:
    provider                = "ftp",
    max_requests_per_second = 1,
    username                = None,
    password                = None,
    active_mode             = False

storage ensembl:
    provider                = "ftp",
    max_requests_per_second = 1,
    username                = None,
    password                = None,
    active_mode             = False

storage ncbi:
    provider                = "ftp",
    max_requests_per_second = 1,
    username                = None,
    password                = None,
    active_mode             = False

storage reich:
    provider                = "http",
    max_requests_per_second = 1,
    auth                    = None,
    allow_redirects         = True

storage github:
    provider                = "http",
    max_requests_per_second = 1,
    auth                    = None,
    allow_redirects         = True

# ------------------------------------------------------------------------------------------------ #
# ---- 1000g-phase3 dataset

rule download_1000_genomes:
    """
    Download 1000genomes phase 3 SNPs
    http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
    """
    input:
        #vcf = FTP.remote(expand("{url}/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        #    url=config["FTP"]["1000g"],
        #    chr = "{chr}"
        #))
        storage.g1k(expand("ftp://{url}/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
            url=config["FTP"]["1000g"],
            chr = "{chr}"
        ))
    output:
        vcf = "data/vcf/1000g-phase3/00-original/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    log: "logs/00-netrules/download_1000_genomes/download_1000_genomes-chr{chr}.log"
    resources:
        cores=lambda w, threads: threads
    threads: 1
    shell: """
        mv {input.vcf} {output.vcf} 2> {log}
    """

rule fetch_samples_panel:
    """
    Download samples metadata from the 1000g FTP website
    """
    input:
        #panel = FTP.remote("ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel")
        panel = storage.g1k("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel")
    output:
        panel = "data/vcf/1000g-phase3/samples-list/integrated_call_samples_v3.20130502.ALL.panel"
    log: "logs/00-netrules/fetch_samples_panel.log"
    resources:
        cores=lambda w, threads: threads
    threads: 1
    shell: """
        mv {input.panel} {output.panel} > {log} 2>&1
    """

# ------------------------------------------------------------------------------------------------ #
# ---- Reich Lab's 1240K-Chip variant callset and panel

rule download_reich_1240K:
    """
    Download the 1240K dataset from Reich Lab's website.
    """
    input:
        #tarball    = HTTP.remote(config["FTP"]["1240K"])
        tarball    = storage.reich(config["FTP"]["1240K"])
    output:
        eigenstrat = multiext("data/Reich-dataset/1240K/v52.2_1240K_public", ".snp", ".ind", ".geno")
    params:
        output_dir = lambda wildcards, output: dirname(output.eigenstrat[0])
    resources:
        cores=lambda w, threads: threads
    log: "logs/00-netrules/download_reich_1240K.log"
    threads: 1
    shell:"""
        tar -xvf {input.tarball} -C {params.output_dir} 2> {log}
    """

# ------------------------------------------------------------------------------------------------ #
# ---- Reference genomes

rule download_reference_genome:
    """
    Download a reference genome from a predefined ftp URL
    human_g1k_v37 is malformed... see: https://github.com/hammerlab/biokepi/issues/117
    """
    input:
        #refgen = FTP.remote("ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/{reference}.fa.gz")
        refgen = storage.ensembl("ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/{reference}.fa.gz")
    output:
        refgen = "data/refgen/GRCh37/{reference}.fa.gz"
    resources:
        cores=lambda w, threads: threads
    log: "logs/00-netrules/download_reference_genome/{reference}.log"
    threads: 1
    shell: """
        mv {input.refgen} {output.refgen}
    """


# ------------------------------------------------------------------------------------------------ #
# ---- Genetic maps

rule download_HapMapII_recombination_map:
    """
    Download the 2011 HapMapII recombination map from ncbi.
    """
    input:
        #tarball = HTTP.remote("http://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz")
        tarball = storage.ncbi("ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz")
    output:
        map     = expand("data/recombination-maps/HapMapII_GRCh37/genetic_map_GRCh37_chr{chr}.txt", chr=range(1, 23)),
        exclude = temp(expand("data/recombination-maps/HapMapII_GRCh37/genetic_map_GRCh37_chr{chr}.txt", chr=["X", "X_par1", "X_par2"])),
        readme  = temp("data/recombination-maps/HapMapII_GRCh37/README.txt")
    params:
        output_dir = lambda wildcards, output: dirname(output.map[0])
    resources:
        cores=lambda w, threads: threads
    log: "logs/00-netrules/download_HapMapII_recombination_map.log"
    threads: 1
    shell: """
        tar -xvzf {input.tarball} -C {params.output_dir} 2>  {log}
        rm {input.tarball}                               2>> {log}
    """


rule fetch_sex_specific_recombination_map:
    """
    See: Bhérer, C., Campbell, C. & Auton, A. Refined genetic maps reveal sexual dimorphism in human meiotic
         recombination at multiple scales. Nat Commun 8, 14994 (2017). https://doi.org/10.1038/ncomms14994

     - https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination
    """
    input:
        #gen_map = HTTP.remote(config["ped-sim"]["input"]["refined-genetic-map-url"])
        gen_map = storage.github("https://" + config["ped-sim"]["input"]["refined-genetic-map-url"])
    output:
        map_dir  = directory("data/recombination-maps/Refined_genetic_map_b37"),
        gen_maps = expand("data/recombination-maps/Refined_genetic_map_b37/{sex}_chr{chrom}.txt", chrom=range(1,23), sex=["female", "male", "sexavg"])
    resources:
        cores = lambda w, threads: threads
    log:     "logs/00-netrules/fetch_sex_specific_gen_map.log"
    threads: 1
    shell: """
        tar --strip-components=1 -xvzf {input.gen_map} -C {output.map_dir} 2> {log} && rm -rf {input.gen_map} >> {log} 2>&1
    """

# ------------------------------------------------------------------------------------------------ #
# ---- Interference maps.

rule fetch_interference_map:
    """
    Download a sex-refined genetic interference map. Required for ped-sim
    See: Campbell, C., Furlotte, N., Eriksson, N. et al. Escape from crossover interference increases with maternal age.
         Nat Commun 6, 6260 (2015). https://doi.org/10.1038/ncomms7260
     - https://github.com/williamslab/ped-sim/blob/master/interfere/nu_p_campbell.tsv
    """
    input:
        #intf_map = HTTP.remote(config["ped-sim"]["input"]["interference-map-url"])
        intf_map = storage.github("https://" + config["ped-sim"]["input"]["interference-map-url"])
    output:
        intf_map = config['ped-sim']['data']['interference']
    resources:
        cores = lambda w, threads: threads
    log: "logs/00-netrules/fetch_sex_specific_gen_map.log"
    threads: 1
    shell: """
        mv {input.intf_map} {output.intf_map} > {log} 2>&1
    """

# ------------------------------------------------------------------------------------------------ #
# ---- Miscellaneous

rule download_TKGWV2_support_files:
    """
    Download Daniel Fernandes' 22M SNP panel from his public google drive.
    """
    output:
        support_files = expand("{directory}/{dataset}", 
            directory = config['kinship']['TKGWV2']['support-files-dir'],
            dataset = [
                "1240K/1000GP3_EUR_1240K.frq",
                "genomeWideVariants_hg19/1000GP3_22M_noFixed_noChr.bed",
                "genomeWideVariants_hg19/DummyDataset_EUR_22M_noFixed.bed",
                "genomeWideVariants_hg19/DummyDataset_EUR_22M_noFixed.bim",
                "genomeWideVariants_hg19/DummyDataset_EUR_22M_noFixed.fam"
            ]
        )
    params:
        url        = config['kinship']['TKGWV2']['support-files-url'],
        output_dir = config['kinship']['TKGWV2']['support-files-dir']
    resources:
        cores = lambda w, threads: threads
    log:     "logs/04-kinship/TKGWV2/TKGWV2_download_support_files.log"
    conda:   "../envs/gdown-4.6.0.yml"
    threads: 1
    shell: """
        gdown "{params.url}" -O {params.output_dir} --folder > {log} 2>&1
    """