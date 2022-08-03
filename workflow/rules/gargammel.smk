configfile: "./config/config.yml"

wildcard_constraints:
    chr="\d+"

rule get_consensus:
    input:
        vcf     = expand(rules.filter_snps.output.vcf, POP="CEU"),
        chr_ref = "data/refgen/splitted/{chr}.fasta"
    output:
        hap1=temp("results/01-gargammel/chr{chr}/{sample}/endo/chr{chr}_{sample}_haplo1.fasta"),
        hap2=temp("results/01-gargammel/chr{chr}/{sample}/endo/chr{chr}_{sample}_haplo2.fasta"),
    group: "scatter"
    threads: 4
    shell:
        """
        bcftools consensus -H 1 -f {input.chr_ref} --sample {wildcards.sample} {input.vcf} > {output.hap1} \
        & \
        bcftools consensus -H 2 -f {input.chr_ref} --sample {wildcards.sample} {input.vcf} > {output.hap2} \
        """

rule get_fastq:
    input:
        hap1=rules.get_consensus.output.hap1,
        hap2=rules.get_consensus.output.hap2
    output:
        forw =temp("results/01-gargammel/chr{chr}_{sample}_s1.fq.gz"),
        rev = temp("results/01-gargammel/chr{chr}_{sample}_s2.fq.gz")
    group: "scatter"
    conda: "../envs/gargammel-1.1.2.yml"
    threads: 4
    priority: 2
    shell:
        """
        ln -sfrt results/01-gargammel/chr{wildcards.chr}/{wildcards.sample} in/garga/cont \
        && \
        ln -sfrt results/01-gargammel/chr{wildcards.chr}/{wildcards.sample} in/garga/bact \
        && \
        gargammel -c 0.005                                                      \
                  --comp 0,0,1                                                  \
                  -l 40                                                         \
                  -o results/01-gargammel/chr{wildcards.chr}_{wildcards.sample} \
                  results/01-gargammel/chr{wildcards.chr}/{wildcards.sample}    
        """


rule merge_chromosomes:
    input:
        forw=expand(rules.get_fastq.output.forw, chr=range(1,23), sample="{sample}"),
        rev =expand(rules.get_fastq.output.rev,  chr=range(1,23), sample="{sample}")
    output:
        forw="results/01-gargammel/fastqs/{sample}_s1.fq.gz",
        rev ="results/01-gargammel/fastqs/{sample}_s2.fq.gz"
    group: "scatter"
    priority: 3
    threads: 4
    shell:
        """
        zcat {input.forw} | gzip > {output.forw}
        zcat {input.rev} | gzip > {output.rev} 
        #rm {input.forw} {input.rev} 
        """
