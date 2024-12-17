# v0.5.1
## Fixes
- Tweak `grups-rs-0.3.2`, `pmd-mask-0.3.2` and `pcangsd-0.99` conda environment definition files
# v0.5.0
## Changes
- **(Breaking)** remove underscores in `config.yml`:
  - In `config.gargammel`: `comp_endo` &rarr; `comp-endo` | `comp_cont` &rarr; `comp-cont` | `comp_bact` &rarr; `comp-bact`)
  - In `config.ped-sim.params`: `error_rate` &rarr; `error-rate`
- **(Breaking)** Rename misleading variables names in `config.yml`:
  - In `config.kinship.TKGWV2`: `min-depth` &rarr; `min-overlap`

# v0.4.2
## Features
- [badger-plots] Implement parsing parameters changes detection procedure
- [badger-plots] Allow user-defined cm-levels during parsing
- Allow plotting of jagged/ragged input data (see: `ragged-input` parameter)
- Allow plotting of individual subcomponents of the plots (see: `plot` parameters)
- [Performance plot]: Allow specification of marker symbols and size (see `markers`, `size` and `symbols` parameters)
- [Accuracy plot]: Allow *visual* transposition of nMBE and nRMSD subplots (see: `flip` parameter)
## Documentation
- [badger-plots] Add documentation entry for `border` parameter
- [badger-plots] Add documentation entry for `split` parameter
- [badger-plots] Add documentation entry for `ragged-input` parameter)
## Fix
- [badger-plots] Fix incorrect CM plotting when row sums equal zero. (which caused division by zero when computing colorscale values)
- [badger][READv2] Fix erroneous setting of `--norm_method` when config value is set to `norm-method: "value"`
- [badger-plots] Fix transposition bug which caused incorrect plotting when setting `transpose=TRUE` on performance plot
- [badger] Fix pip conflicts by explicitly calling the binary found within badger's conda environment
- [badger-plots] Fix incorrect calculation of MBE confidence intervals (normalization value was previously not taken into account when calculating lower and upper bounds)

# CHANGELOG:
 - 2022-12-09: Added support for samtools markdup
 - 2022-10-05: Added ***Limited*** Support for dedup 

