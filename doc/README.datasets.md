# Detailled download instructions

Here, a detailled list of URLs as well as instructions to download each dataset required by BADGER may be found.

## Dataset summary:
| Dataset                                         | version      | URL                                                                                                                                                |
| ----------------------------------------------- | ------------ | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1000g-phase3                                    | v5b-20130502 | http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/                                                                                        |
| AADR 1240K                                      | v52.2        | https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V52/V52.2/SHARE/public.dir/                                               |
| GRCh37 reference genome                         | release-113  | http://ftp.ensembl.org/pub/grch37/release-113/fasta/homo_sapiens/dna/                                                                              |
| HapMapII genetic map                            | 2012         | http://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/                                                                              |
| Refined sex-specific genetic Map (Bhérer et al) | N/A          | https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination                                                                               |
| cross-over interference maps (Campbell et al)   | N/A          | https://github.com/williamslab/ped-sim/blob/1b58e53f23eb61d1e9429c4c14804b2b15d5d928/interfere/nu_p_campbell.tsv                                   |
| TKGWV2 support files (Fernandes 2021)           | N/A          | https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination/blob/0456f74ab8a2dce4d5172914bb799ff3d26d0b43/Refined_genetic_map_b37.tar.gz  |

## Data directory structure.

BADGER currently expects that these datasets adhere to the following directory structure, 
inside of a `./data` directory, located at the root of this repository. 

Using a mix of symbolic and/or hard links should, however, be of no problem at all,
in cases where the user already has these files at hand on their workstation.

```diff
data
├── ped-sim
│   ├── interference_maps
│   │   └── nu_p_campbell.tsv
│   └── recombination_maps
│       └── refined_mf.simmap
├── recombination-maps
│   ├── HapMapII_GRCh37
│   │   ├── genetic_map_GRCh37_chr1.txt
│   │   ├── genetic_map_GRCh37_chr2.txt
│   │   ├── genetic_map_GRCh37_chr3.txt
│   │   ├── genetic_map_GRCh37_chr4.txt
│   │   ├── genetic_map_GRCh37_chr5.txt
│   │   ├── genetic_map_GRCh37_chr6.txt
│   │   ├── genetic_map_GRCh37_chr7.txt
│   │   ├── genetic_map_GRCh37_chr8.txt
│   │   ├── genetic_map_GRCh37_chr9.txt
│   │   ├── genetic_map_GRCh37_chr10.txt
│   │   ├── genetic_map_GRCh37_chr11.txt
│   │   ├── genetic_map_GRCh37_chr12.txt
│   │   ├── genetic_map_GRCh37_chr13.txt
│   │   ├── genetic_map_GRCh37_chr14.txt
│   │   ├── genetic_map_GRCh37_chr15.txt
│   │   ├── genetic_map_GRCh37_chr16.txt
│   │   ├── genetic_map_GRCh37_chr17.txt
│   │   ├── genetic_map_GRCh37_chr18.txt
│   │   ├── genetic_map_GRCh37_chr19.txt
│   │   ├── genetic_map_GRCh37_chr20.txt
│   │   ├── genetic_map_GRCh37_chr21.txt
│   │   └──  genetic_map_GRCh37_chr22.txt
│   └── Refined_genetic_map_b37
│       ├── female_chr1.txt
│       ├── female_chr2.txt
│       ├── female_chr3.txt
│       ├── female_chr4.txt
│       ├── female_chr5.txt
│       ├── female_chr6.txt
│       ├── female_chr7.txt
│       ├── female_chr8.txt
│       ├── female_chr9.txt
│       ├── female_chr10.txt
│       ├── female_chr11.txt
│       ├── female_chr12.txt
│       ├── female_chr13.txt
│       ├── female_chr14.txt
│       ├── female_chr15.txt
│       ├── female_chr16.txt
│       ├── female_chr17.txt
│       ├── female_chr18.txt
│       ├── female_chr19.txt
│       ├── female_chr20.txt
│       ├── female_chr21.txt
│       ├── female_chr22.txt
│       ├── female_chrX.txt
│       ├── male_chr1.txt
│       ├── male_chr2.txt
│       ├── male_chr3.txt
│       ├── male_chr4.txt
│       ├── male_chr5.txt
│       ├── male_chr6.txt
│       ├── male_chr7.txt
│       ├── male_chr8.txt
│       ├── male_chr9.txt
│       ├── male_chr10.txt
│       ├── male_chr11.txt
│       ├── male_chr12.txt
│       ├── male_chr13.txt
│       ├── male_chr14.txt
│       ├── male_chr15.txt
│       ├── male_chr16.txt
│       ├── male_chr17.txt
│       ├── male_chr18.txt
│       ├── male_chr19.txt
│       ├── male_chr20.txt
│       ├── male_chr21.txt
│       ├── male_chr22.txt
│       ├── sexavg_chr1.txt
│       ├── sexavg_chr2.txt
│       ├── sexavg_chr3.txt
│       ├── sexavg_chr4.txt
│       ├── sexavg_chr5.txt
│       ├── sexavg_chr6.txt
│       ├── sexavg_chr7.txt
│       ├── sexavg_chr8.txt
│       ├── sexavg_chr9.txt
│       ├── sexavg_chr10.txt
│       ├── sexavg_chr11.txt
│       ├── sexavg_chr12.txt
│       ├── sexavg_chr13.txt
│       ├── sexavg_chr14.txt
│       ├── sexavg_chr15.txt
│       ├── sexavg_chr16.txt
│       ├── sexavg_chr17.txt
│       ├── sexavg_chr18.txt
│       ├── sexavg_chr19.txt
│       ├── sexavg_chr20.txt
│       ├── sexavg_chr21.txt
│       └── sexavg_chr22.txt
├── refgen
│   └── GRCh37
│       └── Homo_sapiens.GRCh37.dna.primary_assembly.fa
├── Reich-dataset
│   └── 1240K
│       ├── v52.2_1240K_public.anno
│       ├── v52.2_1240K_public.geno
│       ├── v52.2_1240K_public.ind
│       └── v52.2_1240K_public.snp
├── TKGWV2
│   ├── 1240K
│   │   └── 1000GP3_EUR_1240K.frq
│   └── genomeWideVariants_hg19
│       ├── 1000GP3_22M_noFixed_noChr.bed
│       ├── 1000GP3_EUR_22M_noFixed.frq
│       ├── DummyDataset_EUR_22M_noFixed.bed
│       ├── DummyDataset_EUR_22M_noFixed.bim
│       └── DummyDataset_EUR_22M_noFixed.fam
└── vcf
    └── 1000g-phase3
        ├── 00-original
        │   ├── ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr5.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   ├── ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        │   └── ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        └── samples-list
            └── integrated_call_samples_v3.20130502.ALL.panel
```
