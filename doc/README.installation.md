# Manual or custom installation of BADGER

These are detailled installation instructions for users wishing to manually install BADGER within a custom environment.

#### 1. Create a base conda environment for `badger-0.4.1`.
```bash
conda env create -f envs/badger-0.4.1.yml
```

#### 2. Update the base environment with the proper dependencies of `badger-plots`
```bash
conda env update --name badger-0.4.1 --file workflow/scripts/badger-plots/envs/badger-plots.yml
```

#### 3. Install the badger command line program using pip
```bash
conda activate badger-0.4.1
pip install ./badger
```

#### 4. Install the badger.plots library and command line program
```bash
conda activate badger-0.4.1
R --slave -e 'devtools::install("badger/src/badger-plots/badger.plots")'
cp badger/src/badger-plots/badger-plots.R $CONDA_PREFIX/bin/badger-plots
```

#### 5. Configure r-reticulate dependency

```bash
conda activate badger-0.4.1
R_ENVIRON="$(R --slave -e 'cat(R.home(component="home"))')/etc/Renviron"
echo BADGER_PLOTS_CONDA_ENV="badger-0.4.1" >> $R_ENVIRON
echo BADGER_PLOTS_CONDA_EXE=$CONDA_EXE >> $R_ENVIRON
```
