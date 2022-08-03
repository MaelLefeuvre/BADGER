## Execution
```Bash
snakemake --cores `nproc` --use-conda --conda-frontend mamba --describe channel_priority
```

## @TODO:
- `rule bwa_mem` : use `multiext()` for the input.reference
- `rule get_consensus`: error because the output vcf of `rule run_ped_sim` is not BGZF compressed.

### Dependencies:
#### Pedigree simulator>=1.4
 - Source code: https://github.com/williamslab/ped-sim
 - Citation:
     > Caballero M, Seidman DN, Qiao Y, Sannerud J, Dyer TD, Lehman DM, et al. (2019) Crossover interference and sex-specific genetic maps shape identical by descent sharing in close relatives. PLoS Genet 15(12): e1007979. https://doi.org/10.1371/journal.pgen.1007979

 - sub-dependencies: GSL (recommended) (or boost, alternatively)
   `sudo apt-get install libgsl-dev`
 
 - installation example: `resources/utilities/00-install-ped-sim`
