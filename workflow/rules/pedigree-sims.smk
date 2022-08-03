import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

# ---- Set config variables

configfile: "./config/config.yml"


rule fetch_sex_specific_gen_map:
    """
    See: Bhérer, C., Campbell, C. & Auton, A. Refined genetic maps reveal sexual dimorphism in human meiotic
         recombination at multiple scales. Nat Commun 8, 14994 (2017). https://doi.org/10.1038/ncomms14994

     - https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination
    """
    input:
        gen_map = HTTP.remote(config["ped-sim"]["input"]["refined-genetic-map-url"], keep_local=True)
    output:
        map_dir  = directory("data/ped-sim/Refined_genetic_map_b37"),
        gen_maps = expand("data/ped-sim/Refined_genetic_map_b37/{sex}_chr{chrom}.txt", chrom=range(1,23), sex=["female", "male", "sexavg"])
    shell:
        "tar --strip-components=1 -xvzf {input.gen_map} -C {output.map_dir} && rm -rf {input.gen_map}"


rule format_sex_specific_gen_map:
    """
    Download a sex-specific genetic map. Required by ped-sim
    """
    input:
        gen_map = rules.fetch_sex_specific_gen_map.output.gen_maps
    output:
        sim_map = config["ped-sim"]['data']["map"]
    params:
        map_dir = rules.fetch_sex_specific_gen_map.output.map_dir
    shell:
        """
        printf "#chr\tpos\tmale_cM\tfemale_cM\n" > {output.sim_map};                    \
        for chr in {{1..22}}; do                                                        \
            paste {params.map_dir}/male_chr$chr.txt {params.map_dir}/female_chr$chr.txt \
                | awk -v OFS="\t" 'NR > 1 && $2 == $6 {{print $1,$2,$4,$8}}'            \
                | sed 's/^chr//' >> {output.sim_map};                                   \
        done
        """

rule fetch_interference_map:
    """
    Download a sex-refined genetic interference map. Required for ped-sim
    See: Campbell, C., Furlotte, N., Eriksson, N. et al. Escape from crossover interference increases with maternal age.
         Nat Commun 6, 6260 (2015). https://doi.org/10.1038/ncomms7260
     - https://github.com/williamslab/ped-sim/blob/master/interfere/nu_p_campbell.tsv
    """
    input:
        intf_map = HTTP.remote(config["ped-sim"]["input"]["interference-map-url"], keep_local=True)
    output:
        intf_map = config['ped-sim']['data']['interference']
    shell:
        "mv {input.intf_map} {output.intf_map}"


rule run_ped_sim:
    """
    Run pedigree-simulator, using individuals from a specified population of the 1000g project as founder individuals.
    """
    input:
        vcf          = rules.concat_1000_genomes.output.merged_vcf,
        definition   = config['ped-sim']['data']['definition'],
        map          = config['ped-sim']['data']['map'],
        interference = config['ped-sim']['data']['interference']
    output:
        fam = config['ped-sim']['output']['directory'] +"/{POP}-pedigrees-everyone.fam",
        log = config['ped-sim']['output']['directory'] +"/{POP}-pedigrees.log",
        seg = config['ped-sim']['output']['directory'] +"/{POP}-pedigrees.seg",
        vcf = config['ped-sim']['output']['directory'] +"/{POP}-pedigrees.vcf.gz",
    params:
        output_basename  = config['ped-sim']['output']['directory']+ "/{POP}-pedigrees",
        error_rate       = config['ped-sim']['params']['error_rate'],
        missingness      = config['ped-sim']['params']['missingness']
    conda: "../envs/ped-sim-1.3.yml"
    shell:
        """
        ped-sim -d {input.definition}            \
                -m {input.map}                   \
                -i {input.vcf}                   \
                -o {params.output_basename}      \
                --intf {input.interference}      \
                --fam                            \
                --keep_phase                     \
                --miss_rate {params.missingness} \
                --err_rate {params.error_rate}
        #sed -i 's/0/NA/g' {output.fam}
		"""

rule filter_snps:
    input:
        rules.run_ped_sim.output.vcf
    output:
        vcf = config['ped-sim']['output']['directory'] +"/{POP}-pedigrees-M2-m2-snps.vcf.gz"
    threads: 16
    shell:
        "bcftools view --threads {threads} -M2 -m2 -v snps {input} | bgzip -c > {output} && tabix {output}"



checkpoint get_samples:
    input:
        expand(rules.filter_snps.output.vcf, POP=config["ped-sim"]["params"]["POP"])
    output:
        "results/00-ped-sim/sample_names.tsv"
    priority: 50
    shell: 
        """
        zcat {input} | grep '#CHROM' | cut -f 10- > {output}
        """
