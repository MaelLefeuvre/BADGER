# aDNA-kinship-simulations

## Installation

0. Clone this repository
  ```Bash
  user@desktop:~$ git clone --recursive git@github.com:MaelLefeuvre/aDNA-kinship-simulations.git
  ```

1. Install [Conda](https://docs.conda.io/en/latest/)
  - Check the documentation of [miniconda3](https://docs.conda.io/en/latest/miniconda.html) and review the detailled [installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) beforehand.
  - On a x86_64 bits Linux architecture:
    ```Bash
    user@desktop:~$ MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    user@desktop:~$ wget $MINICONDA_URL && bash Miniconda3-latest-Linux-x86_64.sh
    ```

2. [Optional] Install [Mamba](https://github.com/mamba-org/mamba)
  - This integrates seamlessly with conda, and will greatly speed-up its environment solver.
  - On a x86_64 bit Linux architecture:
    ```Bash
    user@desktop:~$ conda install -n base -c conda-forge mamba
    ```

3. Install [Snakemake](https://snakemake.github.io/)
  - Check the [documentation](https://snakemake.readthedocs.io/en/stable/) and review the detailled [installation instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
  - A dedicated environment is available within this repository. On a x86_64 bit Linux architecture:
    ```Bash
    (base) user@desktop:~$ mamba env install -f ./workflow/envs/snakemake-7.12.0.yml
    ```

## Execution
1. Activate snakemake
  ```Bash
mamba activate snakemake-7.12.0.yml
  ``` 
2. Download data through FTP
```Bash
snakemake fetch_data --cores `nproc` --use-conda --conda-frontend mamba --printshellcmds --rerun-incomplete
```

3. Run the pipeline  
  ```Bash
snakemake all --cores `nproc` --use-conda --conda-frontend mamba --printshellcmds --rerun-incomplete
  ```

4. [Optional] Archive the results to the archive directory
```Bash
snakemake archive --cores `nproc` --use-conda --conda-frontend mamba --printshellcmds --rerun-incomplete
```

# Development

## Testing
1. Copy the provided dummy dataset into a `data` directory at the root of the project.
  ```Bash
  (snakemake-7.12.0) user@desktop:~$ cp -r resources/dummy-data ./data
  ```
2. Run the pipeline using the provided `test-config.yml` file
  ```Bash
  (snakemake-7.12.0) user@desktop:~$ snakemake all --cores `nproc` --use-conda --conda-frontend mamba --configfile config/test-config.yml
  ```

## CHANGELOG:
 - 2022-12-09: Added support for samtools markdup
 - 2022-10-05: Added ***Limited*** Support for dedup 

## @TODO:
- [ ] (URGENT) Restructure config file and fully isolate snakefiles into their own submodules. Especially `00-netrules` and `07-archive`
- [ ] (URGENT) Fix bug regarding FTP Permission denied...This is likely caused by snakemake running the connection **before** the checkpoint is validated Thus, the connection is times-out. Generate a MR example and submit an issue to github ?
- [ ] Bump to ped-sim-1.4 and source compile.
- [ ] Finish implementing `dedup`. A bit time consuming, since I need to add FixAdapterRemovalPrefix before running bwa. 
- [ ] There's no need to simulate contaminating genome when gargammel's cont value are @ zero...
- [ ] Unit testing not yet implemented.
- [ ] Use jinja2 engine for pedigree definition file templating.
- [x] Add Indel support for simulations
- [x] Add optional MAF filtration support for pedigree simulations 
- [x] Add microbial database to gargammel.
- [x] `bcftools consensus` has trouble handling indels. We need to find a way to add those.
- [x] `rule bwa_mem` : use `multiext()` for the input.reference
- [x] `rule get_consensus`: error because the output vcf of `rule run_ped_sim` is not BGZF compressed.
- [x] `rule extract_twins` is needlessly complicated and generates unnecessary intermediary files.
- [x] `find_contaminants()`: there's a redundant context manager with a seemingly needless "dummy" file. try to refactor this ? 
 

### Dependencies:
#### [OPTIONAL] Pedigree simulator>=1.4
 - Source code: https://github.com/williamslab/ped-sim
 - Citation:
     > Caballero M, Seidman DN, Qiao Y, Sannerud J, Dyer TD, Lehman DM, et al. (2019) Crossover interference and sex-specific genetic maps shape identical by descent sharing in close relatives. PLoS Genet 15(12): e1007979. https://doi.org/10.1371/journal.pgen.1007979

 - sub-dependencies: GSL (recommended) (or boost, alternatively)
   `sudo apt-get install libgsl-dev`
