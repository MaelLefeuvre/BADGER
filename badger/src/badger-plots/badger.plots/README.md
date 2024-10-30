# badger.plots

`badger.plots` is a pure `R` companion command line interface and library designed to generate plots summarizing the benchmarking results of [BADGER](https://github.com/MaelLefeuvre/BADGER).

## Dependencies

`badger.plots` is currently designed with the intent of being embedded within a self contained environment, and thus will require the use use of an environment management system such as [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/en/latest/#)

Quick command line installation instructions of the miniconda may be found through the link here: [miniconda-quick-install](https://docs.anaconda.com/miniconda/#quick-command-line-install)

## Installation 

The base conda environment of `badger.plots` may be created using the provided `.yml` definition file:
```bash
conda env create -f envs/badger-plots.yml
```

This should create a default conda environment named `badger-plots-0.4.6`

Once successfully created, configuring the environment and embedding the `badger.plots` library within is achieved by simply executing the provided post-deployment bash script
```bash
./envs/badger-plots.post-deploy.sh
```

Once executed, one may ensure the command line interface works by executing the following commands
```bash
conda activate badger-plots-0.4.6
badger-plots help
```

**Expected output**:
```
USAGE: badger-plots <MODULE> [options]

MODULES:
- make-input: Create a .yml specifying inputs for the 'plot' module
- template  : Create a .yml specifying params for the 'plot' module
- plot      : Generate accuracy and performance comparison plots
- help      : Show this help message and exit
```

## Usage

Generating plots using the `badger-plots` command line interface will typically require three steps, which all have their dedicated module.

1. Create a `yaml` file specifying the paths to the archived `BADGER` archive data files. This is done through the `make-input` module.
2. Create and configure a `yaml` file, specifying the parameters for the plots. The `template` module is here to easily provide with a template `yaml` containing all of the parameters, with sensible defaults.
3. Generate the summary plots, using the previously generated `input.yml` and `params.yml` files. This is achieved through the `plot` module.

The following step-by-step procedure assumes you already have a set of archived BADGER results, located within a single archive directory (`path/to/badger/archive-dir`) and with the following structure:

<pre class="ansi2html-content">
  tree -L 1 path/to/badger/archive-dir
  <span style="color: blue">path/to/archive/dir</span>
  ├── <span style="color: blue">Ust-Ishim-GBR-0.02X</span>
  ├── <span style="color: blue">Ust-Ishim-GBR-0.04X</span>
  ├── <span style="color: blue">Ust-Ishim-GBR-0.06X</span>
  ├── <span style="color: blue">Ust-Ishim-GBR-0.08X</span>
  ├── <span style="color: blue">Ust-Ishim-GBR-0.1X</span>
  └── <span style="color: blue">Ust-Ishim-GBR-0.2X</span>
</pre>

In other terms, this example archive directory contains 6 distinct BADGER runs, using a constant Ust-Ishim PMD-damage profile, and randomly sampled `GBR` individuals as pedigree founders. Here the biological impact of genome coverage is most likely to be studied, since every run is parametrized with a different simulated sequencing depth.

---

### 1. Generating an input yaml

Provided these results were archived using `BADGER`'s `archive` target rule, creating an `input.yml` for the `badger-plots` command line interface is simply a matter of specifying the path to the desired archive directory, along with a list of targetted subdirectories.

```bash
badger-plot make-input --archive-dir path/to/badger/archive-dir --subdirs Ust-Ishim-GBR-0.02X,Ust-Ishim-GBR-0.04X,Ust-Ishim-GBR-0.06X > input.yml
```

This will generate an `input.yml`, targetting all of the archived BADGER results found in the three subdirectories `Ust-Ishim-GBR-0.02X`, `Ust-Ishim-GBR-0.04X` and `Ust-Ishim-GBR-0.06X`, ie.:

As the argument `--subdir` leverages regular expressions to match directory suffixes, the previous command may be simplified as
```bash
badger-plot make-input --archive-dir path/to/badger/archive-dir --subdirs 0.02X,0.04X,0.06X,0.08X,0.1X,0.2X > input.yml
```

One may also specify a subset list of kinship estimation methods to target, using the `--methods` argument (e.g. `--methods "KIN,READ,READv2"`)

The generated `input.yml` is expected to be structured according to `subdirs` > `methods` > `[archive-paths]`. e.g.:
```
Ust-Ishim-GBR-0.02X:
  correctKin:
  - path/to/archive/dir/Chan-Meso-CEU-002X/run-000/results/04-kinship/correctKin/ped1.tar.xz
  - path/to/archive/dir/Chan-Meso-CEU-002X/run-000/results/04-kinship/correctKin/ped2.tar.xz
  - ...
  GRUPS:
  - path/to/archive/dir/Chan-Meso-CEU-002X/run-000/results/04-kinship/GRUPS/ped1.tar.xz
  - path/to/archive/dir/Chan-Meso-CEU-002X/run-000/results/04-kinship/GRUPS/ped2.tar.xz
  - ...
  KIN:
  - ...
  READ:
  - ...
  READv2:
  - ...
  TKGWV2:
  - ...
Ust-Ishim-GBR-0.04X:
  ...  
```

### 2. Generate and configure a template plotting parameters yaml file

Plots generated through `badger-plots` can be customized in multiple way, and thus requires the use of a simple `params.yml` file to specify parameters. As most of these parameters contain multiple sensible defaults, a good start is to simply generate a template file, using the dedicated `template` module

```bash
badger-plots template > params.yml
```

We recommend that users directly specify the path to the previously created `input.yml` and `pedigree-codes.txt` files, as these are the only two parameters that are strictly required by `badger.plots`

```bash
badger-plots template -i input.yml -p pedigree-codes.txt > params.yml
```

### 3. Plot the results.

```bash
badger-plots plot --yaml params.yml --threads 32
```

Note that running the `plot` module on the first time may take several minutes, as the program is first tasked with decompressing, parsing and formatting all of the archived BADGER results specified within the `input.yml` file. Preempting several worker threads using the `--threads` is a good way to decrease the runtime.

However, subsequent runs of `badger.plots plot` on the same `params.yml` file (i.e., when tweaking the plotting parameters and design), is expected to be much faster, as the program will store the parsed results within a `data-backup.Rdata`. Thus, leveraging multiple cores using `--threads` is unnecessary when reapplying the `plot` module on the same `params.yml` file.

## Output

At this state, `badger.plots` generates two plots, both provided in the form of an interactive `html` file and a static `svg`:

1. A classification performance plot (`<timestamp>-OCI-performance-plot.[html|svg]`), which summarises confusion matrices for every tested method and biological condition, and highlights their general classification performance through the use of the Uniform Ordinal Classification Index [(Silva W. et al, 2018)](https://doi.org/10.1109/IJCNN.2018.8489327)

2. A set of bartraces, summarizing the normalized Root Mean Square Error and Mean Bias Error of every method, biological condition and degree of relatedness, when attempting r-coefficient estimation.


## Creating a `pedigree-codes.tsv` file

`pedigree-codes.tsv` is a simple unheaded tab-separated file used within badger to specify which relationships should be tested within the simulated pedigrees of the run. In almost every occasion, `badger.plots` will simply require the same `pedigree-codes.tsv` file that was used during the BADGER run the user wishes to plot. Users are thus not expected to extensively modify this file.

Multiple example files may be found and reused as templates within the `resources` directory of BADGER here : [ped-definition](https://github.com/MaelLefeuvre/BADGER/tree/main/resources/ped-sim/ped-definition) 

This file must contains the following four fields.
- `<LABEL>`  User-defined label of the relationship (e.g. "Siblings", "Unrelated", "Half-Sibling").
- `<Ind1>`   ped-sim identifier of the first individual composing the pair (e.g. 'g2-b1-i1')
- `<Ind2>`   ped-sim identifier of the second individual composing the pair (e.g. 'g2-b2-i1')
- `<R>`      the expected theoretical relatedness coefficient for this pair of individuals (e.g. '0.5' for siblings).

**Example**:

```text
#LABEL          Ind1      Ind2      R
Unrelated       g1-b1-i1  g1-b3-i1  0
Parental        g1-b1-i1  g2-b1-i1  0.5
Siblings        g2-b1-i1  g2-b2-i1  0.5
Avuncular       g3-b1-i1  g2-b2-i1  0.25
Half-siblings   g2-b1-i1  g2-b3-i1  0.25
First-cousins   g3-b1-i1  g3-b2-i1  0.125
Half-Avuncular  g2-b3-i1  g3-b2-i1  0.125
```