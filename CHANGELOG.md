# v0.4.2
## Features
- Allow plotting of jagged/ragged input data (see: `ragged-input` parameter)
- Allo plotting of individual subcomponents of the plots (see: `plot` parameters)
- [Performance plot]: Allow specification of marker symbols and size (see `markers`, `size` and `symbols` parameters)
- [Accuracy plot]: Allow *visual* transposition of nMBE and nRMSD subplots (see: `flip` parameter)
# Documentation
- [badger-plots] Add documentation entry for `border` parameter
- [badger-plots] Add documentation entry for `split` parameter
- [badger-plots] Add documentation entry for `ragged-input` parameter)
# Fix
- [badger-plots] Fix incorrect CM plotting when row sums equal zero. (which caused division by zero when computing colorscale values)
- [badger][READv2] Fix erroneous setting of `--norm_method` when config value is set to `norm-method: "value"`
- [badger-plots] Fix transposition bug which caused incorrect plotting when setting `transpose=TRUE` on performance plot
- [badger] Fix pip conflicts by explicitly calling the binary found within badger's conda environment
- [badger-plots] Fix incorrect calculation of MBE confidence intervals (normalization value was previously not taken into account when calculating lower and upper bounds)

# CHANGELOG:
 - 2022-12-09: Added support for samtools markdup
 - 2022-10-05: Added ***Limited*** Support for dedup 

