# BADGER `config.yml` parameters reference

This page lists and describes all of the `badger` parameters that are available and configurable through its [`config/config.yml`](/config/config.yml) configuration file.

## <a name="reference-genome"></a>reference-genome
> <ins>**Description**</ins>: Name of the reference genome used throughout the pipeline.

> <ins>**Allowed values**</ins>: (String)
    - `GRCh37`: Download and use Ensembl's GRCh37 release-113 reference genome.
    - `hs37d5`: Download and use IGSR's hs37d5 reference genome, from the 1000g-phase2 project.
    - `human_g1k_v37`: Download and use IGSR's human_g1k_v37 reference genome from the 1000g-phase3 project.

## <a name="tempdir"></a>tempdir
> <ins>**Description**</ins>: Set a custom temporary directory for methods and softwares requiring the creation of temporary file.

> <ins>**Allowed values**</ins>: (String) Any valid path to a directory.

## <a name="archive"></a>archive
Categorizes arguments related to the `badger archive` module.
- ### <a name="archive-dir"></a>archive-dir
  > <ins>**Description**</ins>: Set a custom directory where BADGER will archive its results, using the `badger archive` module.

  > <ins>**Allowed values**</ins>: (String) Any valid path to a directory.

- ### <a name="compress-level"></a>compress-level
  > <ins>**Description**</ins>: Specify a custom compression level for the [`xz`](https://linux.die.net/man/1/xz) compression algorithm. Higher values will require a higher compression rate implying a longer runtime, but a lighter disk-usage footprint.

  > <ins>**Allowed values**</ins>: (Integer) An integer in the range $[0, 9]$
## <a name="ped-sim"></a>ped-sim
Categorizes arguments related to pedigree simulations, and the `ped-sim` software
- ### <a name="replicates"></a>replicates
  > <ins>**Description**</ins>: Number of pedigree replicates for this given badger run. Increasing this value will optimize the overall runtime of BADGER, but will also linearly increase the number of required snakemake jobs. This behavior will be of little consequence on HPC and computer-cluster environments, but may become unmanageable on single-server environments. Values in the range $[1, 10]$, combined with the use of `loop-pipeline` are thus highly recommended if you plan to run BADGER on a desktop computer.

  > <ins>**Allowed values**</ins>: (Integer) An integer in the range $[1, +\infty[$

- ### <a name="filter-indels"></a>filter-indels
  > <ins>**Description**</ins>: Apply indel filtration on the input `vcf` dataset that is given to `ped-sim` as a source of founder individuals during pedigree simulations.

  > <ins>**Allowed values**</ins>: (boolean) `True` or `False`

- ### <a name="filter-maf"></a>filter-maf
  > <ins>**Description**</ins>: Apply minor allele filtration on the input `vcf` dataset that is given to `ped-sim` as a source of founder individuals during pedigree simulations.

  > <ins>**Allowed values**</ins>: (boolean) `True` or `False`

- ### <a name="filter-maf-threshold"></a>filter-maf-threshold
  > <ins>**Description**</ins>: Minor allele frequency filtration threshold. This value is ignored if [filter-maf](#filter-maf) is set to `False`

  > <ins>**Allowed values**</ins>: (float) any floating point value in the range $[0, 1]$, representing an allele frequency.

- ### <a name="data"></a>data
  Categorizes keywords and arguments providing `ped-sim` with its required input datasets.
  - #### <a name="codes"></a>codes
    > <ins>**Description**</ins>: Path to a user-created pedigree-codes file. Examples of these files can be found in the [`resources/ped-sim/ped-definition`](resources/ped-sim/ped-definition) directory. Detailled specifications on the format and purpose of these files can be found here: [pedigree-codes-files](README.ped-sim-config.md#pedigree-codes-files)

    > <ins>**Allowed values**</ins>:  (String) Any valid path pointing to a pedigree-codes file.

  - #### <a name="definition"></a>definition
    > <ins>**Description**</ins>: Path to a user-created pedigree definition file (`.def`). Examples of these files can be found in the [`resources/ped-sim/ped-definition`](resources/ped-sim/ped-definition) directory.

    > <ins>**Allowed values**</ins>: (String) Any valid path pointing a *"`.def`"* pedigree definition file.

  - #### <a name="map"></a>map
    > <ins>**Description**</ins>: Path to a user-created ped-sim genetic map file. If set to `None`, BADGER will attempt to download and create a default file from the HapMapII recombination map. See the corresponding entry in `ped-sim`'s documentation here: [Map File](https://github.com/williamslab/ped-sim/tree/v1.4.2?tab=readme-ov-file#map-file-)

    > <ins>**Allowed values**</ins>: (String) (Optional) Any valid path pointing to a valid *"`.map`"* genetic map file, or *`None`* (`~` character).

  - #### <a name="interference"></a>interference
    > <ins>**Description**</ins>: Path to a user-created ped-sim *`.intf`* file. If set to `None`, BADGER will attempt to download a default file from `ped-sim`'s code repository. See the corresponding entry in `ped-sim`'s documentation here: [intf file](https://github.com/williamslab/ped-sim/tree/v1.4.2?tab=readme-ov-file#crossover-interference-model---intf-file)

    > <ins>**Allowed values**</ins>: (String) (Optional) Any valid path pointing to a valid *`.intf`* interference file, or *`None`* (`~` character)

- ### <a name="params"></a>params
  Categorizes optional parameters and flags for the `ped-sim` program. These optional parameters are described in a dedicated section within ped-sim's documentation: [Other optional arguments](https://github.com/williamslab/ped-sim/tree/v1.4.2?tab=readme-ov-file#other-optional-arguments). Note that this section does not currently map every optional arguments of `ped-sim` exhaustively.

  - #### <a name="error_rate"></a>error_rate
    > <ins>**Description**</ins>: Sets `ped-sim`'s internal genotyping error rate (`--err_rate` parameter)

    > <ins>**Allowed values**</ins>: (float) Any valid floating point value found in the range $[0, 1]$

  - #### <a name="missingness"></a>missingness
    > <ins>**Description**</ins>: Sets `ped-sim`'s internal rate of missing genotype calls (`--missingness` parameter)

    > <ins>**Allowed values**</ins>: (float) Any valid floating point value found in the range $[0, 1]$
  - #### <a name="retain-extra"></a>retain-extra
    > <ins>**Description**</ins>: Include randomly selected reference individuals as additional samples, which are not related to the main template pedigree file. (`--retain_extra`) This parameter may be useful to evaluate the impact of the number of available Unrelated individuals on the sensitivity of cohort-normalization kinship estimation methods, such as `READ` or `KIN`. 

    > <ins>**Allowed values**</ins>: (Integer) a valid integer number, specifying the number of additional samples to include. Negative values (e.g.: `-1`) will have the impact of printing all of the available input samples.
  - #### <a name="pop"></a>pop
    > <ins>**Description**</ins>: Select which 1000g-project population should be used as a source for founder individuals within `ped-sim`'s pedigree simulations. A detailled list of the available population codes can be seen in IGSR's FAQ section here: [What do the population codes mean ?](https://www.internationalgenome.org/faq/what-do-the-population-codes-mean/)

    > <ins>**Allowed values**</ins>: (String) Any valid 1000 genome population, or super-population label.

  - #### <a name="seed"></a>seed
    > <ins>**Description**</ins>: Provide `ped-sim` with a fixed RNG seed. If set to None (`~`), BADGER will pick a random seed for you, and store it within `results/meta/pipeline-metadata.yml`

    > <ins>**Allowed values**</ins>: (integer) Any valid 64bit integer.

## <a name="simulator"></a>simulator
> <ins>**Description**</ins>: Select which DNA fragment to use. `gargammel` is highly recommended for almost all intents and purposes, `dwgsim` can be useful for debugging purposes, or to simulate modern sequencing data.

> <ins>**Allowed values**</ins>: (String) Either *`"gargammel"`* or *`"dwgsim"`*

## <a name="gargammel"></a>gargammel
This sub-section categorizes all parameters related to the `gargammel` aDNA fragment simulations software. Note that these parameters are ignored unless [`simulator`](#simulator) is set to *`"gargammel"`*
- ### <a name="coverage"></a>coverage
  > <ins>**Description**</ins>: Average simulated sequencing depth (X).

  > <ins>**Allowed values**</ins>: (float) Any non-negative floating point value.

- ### <a name="comp_endo"></a>comp_endo
  > <ins>**Description**</ins>: Fractional composition of endogenous sequences.

  > <ins>**Allowed values**</ins>: (float) Any valid floating point value in the range $[0, 1]$

- ### <a name="comp_cont"></a>comp_cont
  > <ins>**Description**</ins>: Fractional composition of modern human contaminating sequences. Note that this value should be equal to (1 - *[`comp_endo`](#com)*)

  > <ins>**Allowed values**</ins>: (float) Any valid floating point value in the range $[0, 1]$

- ### <a name="comp_bact"></a>comp_bact
  > <ins>**Description**</ins>: Fractional composition of contaminating bacterial sequences. Note that setting this value to anything other than 0 will trigger BADGER to automatically download the example bacterial database proposed in gargammel's source code.

  > <ins>**Allowed values**</ins>: (float) Any valid floating point value in the range $[0, 1]$
- ### <a name="pmd-model"></a>pmd-model
  > <ins>**Description**</ins>: Select which post-mortem damage simulation model to use with `gargammel`.
  > - *`"misincorporation"`: Provide the software with a table of misincorporation occurrences. Usually obtained by applying [`MapDamage`](https://ginolhac.github.io/mapDamage/) on a template sample (`misincorporation.txt` output file).
  > - *`"briggs"`: Simulate post-mortem deamination using the Briggs model *[(Briggs et al 2007)](https://doi.org/10.1073/pnas.0704665104)*.

  > <ins>**Allowed values**</ins>: (String) *`"briggs"`* or *`"misincorporation"`*
- ### <a name="misincorporation"></a>misincorporation
  Categorizes arguments and inputs required to use the `"misincorporation"` post-mortem damage model of gargammel. These arguments are ignored unless [`pmd-model`](#pmd-model) is set to *`"misincorporation"`*
  - #### <a name="file"></a>file
    > <ins>**Description**</ins>: Path to a user-generated *`misincorporation.txt`* file. Note that these files can be generated by applying [`MapDamage`](https://ginolhac.github.io/mapDamage/) on a template sample.

    > <ins>**Allowed values**</ins>: (String) Any valid path pointing to a valid *"`misincorporation.txt`"* file. Usable examples of such files can be found in the [`resources/gargammel/misincorporations/`](resources/gargammel/misincorporations/) directory of this repository.

  - #### <a name="protocol"></a>protocol
    > <ins>**Description**</ins>: Specify whether the misincorporation.txt file was generated from a double-strand or single-strand library.

    > <ins>**Allowed values**</ins>: (String) *`"single"`* or *`"double"`*

- ### <a name="briggs"></a>briggs
  Categorizes arguments and inputs required to use the `"briggs"` post-mortem damage model of gargammel. These arguments are ignored unless [`pmd-model`](#pmd-model) is set to *`"briggs"`*
  - #### <a name="nick-frequency"></a>nick-frequency
    > <ins>**Description**</ins>: Per-base nick-frequency rate ($\nu$ parameter.)

    > <ins>**Allowed values**</ins>:  (float) Any valid floating point value in the range $[0, 1]$

  - #### <a name="overhang-length"></a>overhang-length
    > <ins>**Description**</ins>: length of single-stranded overhangs ($\lambda$ parameter)

    > <ins>**Allowed values**</ins>: (float) Any valid floating point value in the range $[0, 1]$
  - #### <a name="ds-deaminations"></a>ds-deaminations
    > <ins>**Description**</ins>: Rate of double-stranded DNA deaminations ($\delta$ parameter)

    > <ins>**Allowed values**</ins>: (float) Any valid floating point value in the range $[0, 1]$
  - #### <a name="ss-deaminations"></a>ss-deaminations
    > <ins>**Description**</ins>: Rate of single-stranded DNA deaminations ($\delta_{SS}$ parameter)

    > <ins>**Allowed values**</ins>: (float) Any valid floating point value in the range $[0, 1]$
- ### <a name="sizefreq"></a>sizefreq
  > <ins>**Description**</ins>: Path to a user-generated fragment size frequency table. Usable exemples of such files may be found in the [`resources/gargammel/sizefreqs/`](resources/gargammel/sizefreqs/) directory of this repository. This file should be unheaded, and contain two columns
  >  - `<n>`: length of the fragment
  >  - `<prob>`: probability density of occurence, for the given fragment length.

  > <ins>**Allowed values**</ins>: (String) Any valid path pointing to an existing size-frequency distribution table. 

- ### <a name="qshift"></a>qshift
  > <ins>**Description**</ins>: Shift the error rate of reads by a factor of ${1}/{(10^{{qshift}\over{10}})}$. Notice that higher values will carry the effect of ***decreasing*** the error rate on both forward and reverse strands.

  > <ins>**Allowed values**</ins>: (Integer) Any valid integer value in the range $[0, 93]$

- ### <a name="params-1"></a>params
  Categorises parameters related to modern human contamination.
  - #### <a name="contam-pop"></a>contam-pop
    > <ins>**Description**</ins>: Select which 1000g-project' population should be used as a source of modern contaminating individuals within `gargammel`'s aDNA fragment simulations. A detailled list of the available population codes can be seen in IGSR's FAQ section here: [What do the population codes mean ?](https://www.internationalgenome.org/faq/what-do-the-population-codes-mean/). Note that using the same population label as the one found within [`pedsim:params:pop`](#pop) parameter is acceptable, as badger will ensure that any sample previously selected as a founder individual by `ped-sim` is filtered out from the candidate list of `gargammel` when selecting contaminating individual(s).

    > <ins>**Allowed values**</ins>: (String) Any valid 1000 genome population, or super-population label.

  - #### <a name="contaminate-samples"></a>contaminate-samples
    > <ins>**Description**</ins>: If you only wish to contaminate a subset of individuals, provide this parameter with a simple line-separated `.txt` file containing the sample ids of the individuals you wish to contaminate. An example of this file can be seen here: [`resources/gargammel/contaminate-samples-list.txt`](resources/gargammel/contaminate-samples-list.txt)

    > <ins>**Allowed values**</ins>: (String) Any valid path pointing to an existing `.txt` file.

## <a name="dwgsim"></a>dwgsim
Categorises parameters related to the use of the [`dwgsim`](https://github.com/nh13/DWGSIM.git) reads simulator. These parameters are ignored unless [`simulator`](#simulator) is set to *`"dwgsim"`*
- ### <a name="seed-1"></a>seed
  > <ins>**Description**</ins>: Provide `dwgsim` with a fixed RNG seed. If set to None (`~`), BADGER will pick a random seed for you, and store it within `results/meta/pipeline-metadata.yml`

  > <ins>**Allowed values**</ins>: (integer) Any valid 64bit integer.

## <a name="preprocess"></a>preprocess
Categorises parameters related to the alignment of reads and pre-processing of alignment files.
- ### <a name="trimming"></a>trimming
  Categorises parameters related to adapter trimming of reads, prior to alignment.
  - #### <a name="min-overlap"></a>min-overlap
    > <ins>**Description**</ins>: Only trim adapters having an overlap of $n$ nucleotides with its corresponding read (not counting ambiguous nucleotides). See the corresponding section in the documentation of `AdapterRemoval`: [`--minadapteroverlap`](https://adapterremoval.readthedocs.io/en/2.3.x/manpage.html#cmdoption-AdapterRemoval-minadapteroverlap)

    > <ins>**Allowed values**</ins>: (integer) any positive non-zero integer representing a fragment length.
  - #### <a name="min-length"></a>min-length
    > <ins>**Description**</ins>: Discard any read shorter than the provided length following adapter trimming. See the corresponding section in the documentation of `AdapterRemoval`: [`--minlength`](https://adapterremoval.readthedocs.io/en/2.3.x/manpage.html#cmdoption-AdapterRemoval-minlength)

    > <ins>**Allowed values**</ins>: (integer) any positive non-zero integer, representing a fragment length.
  - #### <a name="min-quality"></a>min-quality
    > <ins>**Description**</ins>: Trim low-quality bases that are lower than the provided threshold. See the corresponding section in the documentation of `AdapterRemoval`: [`--minquality`](https://adapterremoval.readthedocs.io/en/2.3.x/manpage.html#cmdoption-AdapterRemoval-minquality)

    > <ins>**Allowed values**</ins>: (integer) any positive non-zero integer, representing a PHRED scale base quality score [(Ewing and Green, 1998)](https://doi.org/10.1101/gr.8.3.186).

  - #### <a name="qualitymax"></a>qualitymax
    > <ins>**Description**</ins>: Specifies the maximum expected Phred score found within the input fastq files. Modifying this parameter should only be useful for debugging purposes and/or when modifying [gargammel:qshift](#qshift) See the corresponding section in the documentation of `AdapterRemoval`: [`--qualitymax`](https://adapterremoval.readthedocs.io/en/2.3.x/manpage.html#cmdoption-AdapterRemoval-qualitymax)

    > <ins>**Allowed values**</ins>: (integer) any positive non-zero integer, representing a PHRED scale base quality score [(Ewing and Green, 1998)](https://doi.org/10.1101/gr.8.3.186).

  - #### <a name="seed-2"></a>seed
    > <ins>**Description**</ins>: Provide `AdapterRemoval` with a fixed RNG seed. If set to None (`~`), BADGER will pick a random seed for you, and store it within `results/meta/pipeline-metadata.yml`. See the corresponding section in the documentation of `AdapterRemoval`: [`--seed`](https://adapterremoval.readthedocs.io/en/2.3.x/manpage.html#cmdoption-AdapterRemoval-seed)

    > <ins>**Allowed values**</ins>: (integer) Any valid 64bit integer.

- ### <a name="bwa"></a>bwa
  Categorises parameters for the burrows-wheeler alignment software. See the documentation of `bwa` here: [`Burrows-Wheeler Aligner`](https://bio-bwa.sourceforge.net/)
  - #### <a name="aligner"></a>aligner
    > <ins>**Description**</ins>: Select the main alignment algorithm for `bwa`:
    > - *`"aln"`*: Run the `bwa aln` gapped/ungapped alignment algorithm, followed by `samse`|`sampe` to generate alignment. This algorithm is probably the sounder option out of the two, and is recommended for most cases, but is much more demanding in terms of average runtime.
    > - *`"mem"`*: Run the `bwa mem` algorithm. Note that this method is much faster than `bwa aln`, but is expected to return less accurate results *[(Oliva et al. 2021)](https://doi.org/10.1002/ece3.8297)*.

    > <ins>**Allowed values**</ins>: (string) *`"aln"`* or *`"mem"`*

  - #### <a name="collapsed-only"></a>collapsed-only
    > <ins>**Description**</ins>: Ignore and filter out uncollapsed reads before alignment when simulating paired-end data.

    > <ins>**Allowed values**</ins>: `True`|`False`

  - #### <a name="bwa-aln"></a>bwa-aln
    Categorises parameters for the `aln` module of `bwa`. These parameters are ignored unless the [aligner](#aligner) parameter value is set to *`"aln"`*. For a good starting set of recommendations regarding these parameters, see *[(Oliva et al. 2021)](https://doi.org/10.1093/bib/bbab076)*

  - ##### <a name="seed-length"></a>seed-length
    > <ins>**Description**</ins>: Minimum seed length (this parameter is mapped to the `-l` argument of `bwa aln`). The provided integer value will require `bwa` to use subsequence of the corresponding length as seed. However, if the provided value is ***larger*** than the query sequence, then seeding will be disabled altogether. This latter behavior is usually recommended in most use-cases when aligning aDNA whole-genome shotgun sequencing data.

    > <ins>**Allowed values**</ins>: (integer) Any positive non-zero integer, representing a seed length (in bases).

  - ##### <a name="max-open-gap"></a>max-open-gap
    > <ins>**Description**</ins>: Maximum number of gap openings when aligning a fragment (i.e. a continuous indel sequence found along the fragment). This parameter is mapped to the `-o` argument of `bwa aln`. 

    > <ins>**Allowed values**</ins>: (integer) Any positive non-zero integer, representing a maximum allowed number of occurrences.

  - ##### <a name="max-seed-diff"></a>max-seed-diff
    > <ins>**Description**</ins>: Maximum within-seed edit distance. This parameter is mapped to the `-k` argument of `bwa aln`.

    > <ins>**Allowed values**</ins>: Any positive non-zero integer, representing a maximum number of allowed substitutions

  - ##### <a name="max-miss-prob"></a>max-miss-prob
    > <ins>**Description**</ins>: This parameter is mapped to the `-n` argument of `bwa aln`, which can either control the maximum within-fragment edit distance, if the provided value is an integer, or the maximum fraction of missing alignments if the provided value is a floating point number. In the latter case, the maximum allowed edit distance is chosen dynamically, according to the length of the fragment.

    > <ins>**Allowed values**</ins>: Either a floating point value in the range $[0, 1]$, representing a fraction of missing alignments, or an integer value in the range $[0, +\infty]$, representing the maximum allowed Levenstein distance. 

- ### <a name="filter"></a>filter
  Categorises post-alignment quality filtration parameters.
  - #### <a name="min-mq"></a>min-MQ
    > <ins>**Description**</ins>: Filter out any fragment carrying a mapping quality that is lower than the provided threshold

    > <ins>**Allowed values**</ins>: (integer) Any positive non-zero integer in the range $[0, 33]$, representing a PHRED33 mapping quality score

  - #### <a name="min-length-1"></a>min-length
    > <ins>**Description**</ins>: Filter out any fragment with a length lower than the provided threshold
  
    > <ins>**Allowed values**</ins>: (integer) Any positive non-zero integer in the range $[0, +\infty]$, representing a fragment length.
- ### <a name="dedup"></a>dedup
  Categorises arguments related to PCR and/or optical duplicates removal softwares.
  - #### <a name="method"></a>method
    > <ins>**Description**</ins>:  Select the main duplicate removal program:
    > - *`"samtools"`*: Mark and remove duplicates using `samtools markdup` *[(Danecek et al 2021)](https://doi.org/10.1093/gigascience/giab008)*. See the corresponding documentation here: [markdup](http://www.htslib.org/doc/1.15/samtools-markdup.html)
    > - *`"picard"`*: Mark and remove duplicates using `picard MarkDuplicates` *[(Broad Institute)](http://broadinstitute.github.io/picard/)*. See the corresponding documentation here: [MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) 
    > - *`"dedup"`*: Mark and remove duplicates using `DeDup` *[(Peltzer et al. 2016)](https://doi.org/10.1186/s13059-016-0918-z)*. See the corresponding documentation here: [Dedup](https://dedup.readthedocs.io/en/latest/)

    > <ins>**Allowed values**</ins>: (string) one of *`"samtools"`*|*`"picard"`*|*`"dedup"`*
- ### <a name="pmd-rescaling"></a>pmd-rescaling
  Categorises arguments related to post-mortem damage estimation and rescaling softwares
  - #### <a name="rescaler"></a>rescaler
    > <ins>**Description**</ins>:  Select the main pmd-rescaling software to use. This pre-processing step is optional and may be replaced by a null value (*`~`*), to skip pmd-rescaling altogether.
    > - *`"mapdamage"`*: Estimate and rescale post-mortem deaminations using `mapDamagev2` *[(Jо́nsson et al. 2013)](https://doi.org/10.1093/bioinformatics/btt193)*. See the corresponding documentation here: [mapDamagev2](https://ginolhac.github.io/mapDamage/).
    > - *`"pmdtools"`*: Estimate and rescale post-mortem deaminations using `PMDtools` *[(Skoglund et al. 2014)](https://doi.org/10.1073/pnas.1318934111)* . See the corresponding documentation here: [PMDtools](https://github.com/pontussk/PMDtools/blob/master/README.md).
    > <ins>**Allowed values**</ins>: (string)(optional) either *`"mapdamage"`*, *`"pmdtools"`*, or *`None`* (*`~`* character).

  - #### <a name="apply-masking"></a>apply-masking
    > <ins>**Description**</ins>: Boolean specifying if putatively deaminated nucleotides should be masked using the software `pmd-mask`. Note that masking terminal deaminations may be applied on either non-rescaled or previously rescaled bam files, depending on whether or not `rescaler` was set to None`~`. In both cases, setting `apply-masking` to *`True`* will trigger the execution of `MapDamagev2` regardless, as `pmd-mask` currently requires the misincorporation frequencies estimates of this software as an input.

    > <ins>**Allowed values**</ins>: (boolean) either *`True`* or *`False`*

  - #### <a name="map-damage"></a>map-damage
    Categorises arguments related to the `mapDamage-v2` software. These parameters are ignored unless [rescaler](#rescaler) is set to *`"mapdamage"`* or [apply-masking](#apply-masking) is set to *`True`*.
    - ##### <a name="rescale"></a>rescale
      > <ins>**Description**</ins>: Specify whether or not pmd-rescaling should be applied. Setting this parameter to `False` may be useful to run mapdamage in order to only obtain PMD-estimates as a form of QC-validation.

      > <ins>**Allowed values**</ins>: (boolean) either *`True`* or *`False`*.
    - ##### <a name="downsample"></a>downsample
      > <ins>**Description**</ins>: Estimate deamination frequencies on a randomly selected fraction of the file (if the provided value is a floating point), or a fixed number of reads (if the provided value is an integer). This can speed up computation, but may return inaccurate estimates.

      > <ins>**Allowed values**</ins>: (float) A floating point value in the range $[0, 1]$, representing a fraction of the number of reads within a given file, or a non-negative integer value in the range $[0, +\infty]$, representing a set number of reads.

    - ##### <a name="downsample-seed"></a>downsample-seed
      > <ins>**Description**</ins>: Provide `MapDamagev2` with a fixed RNG seed for downsampling. If set to None (`~`), BADGER will pick a random seed for you, and store it within `results/meta/pipeline-metadata.yml`.

      > <ins>**Allowed values**</ins>: (integer) Any non-negative 32bit integer value.

  - #### <a name="pmdtools"></a>pmdtools
    Categorises arguments related to the `PMDtools` software. These parameters are ignored unless [rescaler](#rescaler) is set to *`"pmdtools"`*
    - ##### <a name="threshold"></a>threshold
      > <ins>**Description**</ins>: Filter out reads carrying a PMD-score that is lower than the provided threshold. This parameter may be useful to filter out reads originating from modern human contamination. Note that higher PMD-score values indicate a higher likelihood of a fragment having been subject do post-mortem damage. Setting this parameter to a value that is less than or equal to *`-19999`* will disable PMD-score filtering altogether.

      > <ins>**Allowed values**</ins>: (integer) Any integer value in the range $[-\infty,  +\infty]$, representing a Log-likelihood for the presence of PMD contamination.

    - ##### <a name="mask-terminal-deams"></a>mask-terminal-deams
      > <ins>**Description**</ins>: Mask a set number of putatively deaminated nucleotides along fragments.

      > <ins>**Allowed values**</ins>: (integer) Any positive integer in the range $[0, +\infty]$, representing a number of nucleotides.

  - #### <a name="pmd-mask"></a>pmd-mask
    Categorises arguments related to the `pmd-mask` software. These parameters are ignored unless [apply-masking](#apply-masking) is set to *`True`*
    - ##### <a name="threshold-1"></a>threshold
      > <ins>**Description**</ins>: Mask all `5'C` and `3'G` along the reads, until the estimated nucleotide misincorporation rate is lower than the provided threshold

      > <ins>**Allowed values**</ins>: (float) Any floating point value in the range $[0, 1]$, representing a rate of misincorporation.

## <a name="variant-calling"></a>variant-calling
Categorises parameters related to random pseudo-haploid variant calling.
- ### <a name="caller"></a>caller
  > **Description:** Select the main random pseudo-haploid variant caller to use.
  > - *`"pileupCaller"`*: Apply random pseudo-haploid variant calling using `sequenceTools pileupCaller` program *[(Schiffels, 2015)](https://github.com/stschiff/sequenceTools)*. 
  > - *`"ANGSD"`*: Apply random pseudo-haploid variant calling using the `ANGSD` program *[(Korneliussen et al. 2014)](https://doi.org/10.1186/s12859-014-0356-4)*. See the documentation of the corresponding haploid calling module here: [`ANGSD -doHaploCall`](https://www.popgen.dk/angsd/index.php/Haploid_calling)

  > **Allowed values:** (String) Either *`"pileupCaller"`* or *`"ANGSD"`*

- ### <a name="maf"></a>maf
  > **Description:** Pre-filter the provided SNP callset by removing variants that are below the provided minor allele-frequency threshold, when looking at the 1000g-phase3 dataset.

  > **Allowed values:** Any non-negative floating point value in the range $[0, 1]$, representing a minor allele frequency threshold.

- ### <a name="maf-superpop"></a>maf-superpop
  > **Description:** Specify which super-population is used when evaluating allele frequencies, during pre-filtration step. (See [variant-calling:maf](#maf-1))

  > **Allowed values:** (String) A three-letter 1000genomes super-population tag, i.e: *`"AFR"`*, *`"AMR"`*,*`"EAS"`*,*`"EUR"`*, or *`"SAS"`*.

- ### <a name="pileup"></a>pileup
  Categorises arguments related to the `samtools mpileup` module. The extensive documentation of this specific module of the `samtools` program can be found here: [mpileup](http://www.htslib.org/doc/1.15/samtools-mpileup.html)
  - #### <a name="disable-baq"></a>disable-BAQ
    > **Description:** Specify whether per-Base Alignment Quality should be recalculated. Setting this value to `True` is recommended in most use-cases, when simulating ancient DNA fragments.

    > **Allowed values:** (boolean) *`True`*, or *`False`*

  - #### <a name="min-bq"></a>min-BQ
    > **Description:**  Skip nucleotides carrying a base quality that is lower than the provided threshold

    > **Allowed values:** (integer) a non-negative integer in the range $[0, 40]$, representing a PHRED quality score.

  - #### <a name="min-mq-1"></a>min-MQ
    > **Description:** Skip alignments carrying a mapping quality that is lower than the provided threshold

    > **Allowed values:** (integer) a non-negative integer in the range $[0, 40]$, representing a PHRED quality score.

- ### <a name="pileupcaller"></a>pileupCaller
  Categorises arguments related to the `pileupCaller` random pseudo-haploid variant calling software. These parameters are ignored unless [caller](#caller) is set to *`"pileupCaller"`*
  - #### <a name="skip-transitions"></a>skip-transitions
    > **Description:** Filter out transitions from the SNP variant callset.

    > **Allowed values:** (boolean) *`True`* or *`False`*

  - #### <a name="mode"></a>mode
    > **Description:** Select the main variant calling mode of `pileupCaller`:
    > - *`"randomHaploid"`*: apply random-haploid variant calling. This mode is recommended for most use-cases.
    > - *`"majorityCall"`*: Always pick the allele that has the highest number of observations. Ties are resolved through random selection.

    > **Allowed values:** *`"randomHaploid"`* or *`"majorityCall"`*

  - #### <a name="min-depth"></a>min-depth
    > **Description:** Filter out positions where the local sequencing depth is lower than the provided threshold.

    > **Allowed values:** (integer) Any non-negative number in the range $[1, +\infty]$, representing a sequencing depth.

  - #### <a name="seed-3"></a>seed
    > **Description:**  Provide `pileupCaller` with a fixed RNG seed when applying random sampling. If set to None (`~`), BADGER will pick a random seed for you, and store it within `results/meta/pipeline-metadata.yml`.

    > **Allowed values:** (integer) Any non-negative 32bit integer value.

## <a name="kinship"></a>kinship
Categorises arguments related to all of the benchmarked kinship estimation methods
- ### <a name="targets"></a>targets
  > **Description:** Path pointing to a user-provided list of bi-allelic SNP targets file (preferably in EIGENSTRAT `.snp` format). If set to `None`, BADGER will attempt to download and preprocess the *AADR-1240K-v52.2* SNP dataset as a default. More information regarding this default dataset may be found here: [AADR-v52.2](https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/index_v52.2.html).

  > **Allowed values:** (String) (optional) Any valid path pointing to an existing and valid *"`.snp`"* SNP variant callset file, or *`None`* (`~` character).

- ### <a name="exclude-samples"></a>exclude-samples
  > **Description:** Specify a list of ped-sim samples to exclude from the benchmark during kinship estimation. This parameter may be useful to apply kinship estimation on a subset of a larger pedigree, or to estimate the impact of estimating kinship when key individuals are unsampled.

  > **Allowed values:** (List<Strings>) A list of ped-sim sample ids, (e.g.: `["ped1_g1-b1-i1", "ped1_g1-b2-i1", "ped1_g1-b3-i1"]`)

- ### <a name="correctkin"></a>correctKin
  Categorises arguments related to the `correctKin` kinship estimation software *[(Nyerki et al. 2024)](https://doi.org/10.1186/s13059-023-02882-4)*
  - #### <a name="deplete-indivs"></a>deplete-indivs
    > **Description:** Apply random depletion of genotyped markers across all samples and reference genomes, to ensure equal average coverage.

    > **Allowed values:** (boolean) either *`True`* or *`False`*

  - #### <a name="reference-pop"></a>reference-pop
    > **Description:** Specify which 1000g-phase3 population should be used as a reference population

    > **Allowed values:** (String) Any valid three-lettered 1000g-phase3 population, or super-population label. A detailled list of the available population codes can be seen in IGSR's FAQ section here: [What do the population codes mean ?](https://www.internationalgenome.org/faq/what-do-the-population-codes-mean/) 

- ### <a name="grups"></a>GRUPS
  Categorises arguments related to the `GRUPS-rs` kinship estimation software *[(Lefeuvre et al. 2024)](https://doi.org/10.47248/hpgg2404010001)*
  - #### <a name="pedigree"></a>pedigree
    > **Description:** Path pointing to a `GRUPS-rs` pedigree definition file. Useable examples of such files may be found in this repository's resource directory ([resources/grups-rs/pedigrees/](/resources/grups-rs/pedigrees/)). For users wishing to provide GRUPS-rs with custom pedigrees, detailled specifications of the format may be found in GRUPS-rs' documentation here : [Defining custom pedigrees](https://github.com/MaelLefeuvre/grups-rs/tree/v0.3.2#Defining-custom-pedigrees)

    > **Allowed values:** (Path) (Optional) Any valid path pointing to an existing valid GRUPS-rs pedigree definition file. 

  - #### <a name="pedigree-pop"></a>pedigree-pop
    > **Description:** Specify which 1000g-phase3 population should be used as a reference population, within GRUPS-rs' internal pedigree simulations.

    > **Allowed values:** (String) Any valid three-lettered 1000g-phase3 population, or super-population label. A detailled list of the available population codes can be seen in IGSR's FAQ section here: [What do the population codes mean ?](https://www.internationalgenome.org/faq/what-do-the-population-codes-mean/) 

  - #### <a name="contam-pop-1"></a>contam-pop
    > **Description:** Specify which 1000g-phase3 population should be used as a source for modern human contamination, within GRUPS-rs' internal pedigree simulations.

    > **Allowed values:** (String) Any valid three-lettered 1000g-phase3 population, or super-population label. A detailled list of the available population codes can be seen in IGSR's FAQ section here: [What do the population codes mean ?](https://www.internationalgenome.org/faq/what-do-the-population-codes-mean/) 

  - #### <a name="min-depth-1"></a>min-depth
    > **Description:** Specify the minimum allowed sequencing depth, when comparing two overlapping positions on a pair of individuals.
    
    > **Allowed values:** (integer) any integer value in the range $[1, +\infty]$. 
  
  - #### <a name="mode-1"></a>mode
    > **Description:** Specify the input data format, and memory model of GRUPS-rs
    > - *`"fst"`*: Encode the 1000g-phase3 dataset as a set of FSA-encoded files; load these files in memory at runtime. This mode is recommended for most use cases, as it provides good runtime performance for most use-cases. 
    > - *`"fst-mmap"`*: Encode the 1000g-phase3 dataset as a set of FSA-encoded files; access these files at runtime using memory-mapped files. This mode is highly recommended when the data is stored on an SSD hard-drive, as this provides with the best runtime performance, while mitigating memory usage.
    > - *`"vcf"`*: Directly use the 1000g-phase3 dataset in `.vcf` format. Not recommended unless for debugging purposes.

    > **Allowed values:** (String) Either *"`fst`"*, *"`fst-mmap`"*, or *"`vcf`"*

  - #### <a name="reps"></a>reps
    > **Description:** Specify the number of simulation replicates for GRUPS-rs.

    > **Allowed values:** (integer) Any integer value in the range $[1, +\infty]$.

  - #### <a name="maf-1"></a>maf
    > **Description:** Apply minor allele filtration of overlapping positions, allele frequencies from the 1000g phase3 project. Setting this value to 0 will disable maf filtration altogether.

    > **Allowed values:** (float) Any floating point value in the range $[0, 1]$

  - #### <a name="min-qual"></a>min-qual
    > **Description:** Ignore positions and nucleotides whose base quality is found lower than the provided threshold, when comparing two individuals.

    > **Allowed values:** (integer) Any non-negative value in the range $[0, 40]$, representing a PHRED base-quality score

  - #### <a name="seq-error-rate"></a>seq-error-rate
    > **Description:** Specify a rate of simulated sequencing error rate. When set to `None` (`~`), `GRUPS-rs` will instead simulate sequencing errors using the average base quality scores of the sample, at every position.

    > **Allowed values:** (Float) (Optional) Any non-negative floating point in the range $[0, 1]$, representing a rate of sequencing errors, or `None` (`~`).

  - #### <a name="seed-4"></a>seed
    > **Description:** Provide `GRUPS-rs` with a fixed RNG seed during pedigree simulations. If set to None (`~`), BADGER will pick a random seed for you, and store it within `results/meta/pipeline-metadata.yml`.

    > **Allowed values:** (integer) Any non-negative 32bit integer value. 

- ### <a name="kin"></a>KIN
  Categorises arguments related to the `KIN` kinship estimation software *[(Popli et al. 2023)](https://doi.org/10.1186/s13059-023-02847-7)*
  - #### <a name="interval"></a>interval
    > **Description:** Specify the length of the non overlapping sliding windows when inferring patterns of IBD-sharing between two individuals. This argument is passed on to the `KINgaroo` module of the software. (Default recommended value is `10 000 000`).

    > **Allowed values:** (integer) Any integer in the range $[1, +\infty]$, representing a genomic length.

  - #### <a name="p0-threshold"></a>p0-threshold
    > **Description:** Minimum number of distinct  non-zero genomic windows for a sample to be included when estimating $P_{0}$ (i.e. the median number of pairwise differences when comparing all samples). This argument is passed on to the `KINgaroo` module of the software. (Default recommended value is `10`)

    > **Allowed values:** (integer) Any value in the range $[1, +\infty]$, representing a number of non-zero genomic windows.

  - #### <a name="roh-threshold"></a>roh-threshold
    > **Description:** Minimum required number of genotyped positions within a given genomic window. This argument is passed on to `KIN`. (Default recommended value is `10`)

    > **Allowed values:** (integer) Any value in the range $[1, +\infty]$.

  - #### <a name="diversity-parameter"></a>diversity-parameter
    > **Description:** Provide the software with a user-defined estimate of $P_{0}$. Setting this value to `None` (`~`) will order KIN to estimate this value from the data (default.)

    > **Allowed values:** (float) (optional) Any floating point value in the range $[0,1]$ representing a median genetic distance between two individuals.

  - #### <a name="noisy-windows"></a>noisy-windows
    > **Description:** Optionally provide `KINgaroo` with a list of 0-based genomic windows to filter them out.

    > **Allowed values:** (Path) (optional) Any path pointing to an existing and valid file containing 0-based window indexes, or `None` (`~`).

- ### <a name="read"></a>READ
  Categorises arguments related to the `READ` kinship estimation software *[(Kuhn et al. 2018)](https://doi.org/10.1371/journal.pone.0195491)*
  - #### <a name="window-size"></a>window-size
    > **Description:** Specify the length of the non-overlapping sliding windows when estimating $P_{0}$. The default recommended value is `1 000 000`

    > **Allowed values:** (integer) Any value in the range $[1, +\infty]$, representing a genomic length.

  - #### <a name="norm-method"></a>norm-method
    > **Description:** Specify which summarry statistic should be used to compute a normalized $\overline{P}_{0}$:
    > - *`"median"`*: Use the median of $P_{0}$ across all pairwise comparisons.
    > - *`"mean"`*: Use the mean of $P_{0}$ across all pairwise comparisons.
    > - *`"max"`* Use the maximum computed $P_{0}$ across all pairwise comparisons.
    > - *`"value"`* Specify a user-defined value for normalization. (See [norm-value](#norm-value))

    > **Allowed values:** (String) One of either `"median"`, `"mean"`, `"max"`, or `"value"`.

  - #### <a name="norm-value"></a>norm-value
    > **Description:** Specify a set $P_{0}$ value to compute normalized $\overline{P}_{0}$ estimates. This value is ignored unless [norm-method](#norm-method) is set to *`"value"`*

    > **Allowed values:** (float) (optional) Any floating point value in the range $[0, 1]$

- ### <a name="readv2"></a>READv2
  Categorises arguments related to the `READv2` kinship estimation software *[(Alaçamlı et al. 2024)](https://doi.org/10.1186/s13059-024-03350-3)*
  - #### <a name="window-est"></a>window-est
    > **Description:** Require a window-based estimation of $P_{0}$ values, as opposed to a genome-wide estimate. Default recommended setting is `False`. Setting this value to `True` essentially reverts back to the behaviour of `READv1`.

    > **Allowed values:** (boolean) either `True` or `False`

  - #### <a name="window-size-1"></a>window-size
    > **Description:** Specify the length of non-overlapping sliding windows when estimating $P_{0}$. The default recommended value is `5 000 000`

    > **Allowed values:** (integer) Any integer value in the range $[1, +\infty]$, representing a genomic length.

  - #### <a name="norm-method-1"></a>norm-method
    > **Description:** Specify which summarry statistic should be used to compute a normalized $\overline{P}_{0}$:
    > - *`"median"`*: Use the median of $P_{0}$ across all pairwise comparisons.
    > - *`"mean"`*: Use the mean of $P_{0}$ across all pairwise comparisons.
    > - *`"max"`* Use the maximum computed $P_{0}$ across all pairwise comparisons.
    > - *`"value"`* Specify a user-defined value for normalization. (See [norm-value](#norm-value))

    > **Allowed values:** (String) One of either `"median"`, `"mean"`, `"max"`, or `"value"`.
  - #### <a name="norm-value-1"></a>norm-value
    > **Description:** Specify a set $P_{0}$ value to compute normalized $\overline{P}_{0}$ estimates. This value is ignored unless [norm-method](#norm-method) is set to *`"value"`*

    > **Allowed values:** (float) (optional) Any floating point value in the range $[0, 1]$

  - #### <a name="2pow"></a>2pow
    > **Description:** Request the use of alternate kinship classification thresholds, using powers of two.
    > - When set to `False` (this is the recommended default), the program will use the classification thresholds described in the paper, i.e.:
    >    | Degree  | Threshold |
    >    | ------- | --------- |
    >    | `3°`    | `0.96875` |
    >    | `2°`    | `0.90625` |
    >    | `1°`    | `0.8125`  |
    >    | `Twins` | `0.625`   | 
    > - When set to `True`, the program will instead use the following thresholds:
    >   | Degree  | Threshold      |
    >   | ------- | -------------- |
    >   | `3°`    | `1-1/(2**4.5)` |
    >   | `2°`    | `1-1/(2**3.5)` | 
    >   | `1°`    | `1-1/(2**2.5)` |
    >   | `Twins` | `1-1/(2**1.5)` |

    > **Allowed values:** (boolean) either `True` or `False`
- ### <a name="tkgwv2"></a>TKGWV2
  Categorises arguments related to the `TKGWV2` kinship estimation software *[(Fernandes et al. 2021)](https://doi.org/10.1038/s41598-021-00581-3)*
  - #### <a name="downsample-1"></a>downsample
    > **Description:** Apply read downsampling of the input BAM files prior to using `TKGWV2` (using the `downsampleBam.R` helper script provided by the authors)

    > **Allowed values:** (boolean) eiher `True` or `False`

  - #### <a name="downsample-n"></a>downsample-N
    > **Description:** Specify the number of requested remaining aligned fragments after downsampling. This parameter is ignored unless [downsample](#downsample-1) is set to `True`.

    > **Allowed values:** (integer) (optional) Any integer value in the range $[1, +\infty]$, representing a number of aligned reads.

  - #### <a name="downsample-seed-1"></a>downsample-seed
    > **Description:** Provide `TKGWV2` with a fixed RNG seed during downsampling. If set to None (`~`), BADGER will pick a random seed for you, and store it within `results/meta/pipeline-metadata.yml`. 

    > **Allowed values:** (integer) Any non-negative 32bit integer value. 

  - #### <a name="target-frequencies"></a>target-frequencies
    > <ins>**Description**</ins>: Path pointing to a user-provided file containing allele frequencies, in PLINK `.freq` format (See the specifications of this file format here: [.frq file format](https://www.cog-genomics.org/plink/1.9/formats#frq)). If set to `None`, BADGER will attempt to download and use the example dataset provided by the authors of `TKGWV2`, which contains allele frequencies of the 1000g-phase3 `EUR` population, along bi-allelic positions of the *AADR-1240K* dataset.. More information regarding this default dataset may be found here: [Support files](https://github.com/danimfernandes/tkgwv2/blob/175b0d6674b7b19c5833b66b8a801ea9a2716b51/README.md#description-and-generation-of-support-files).

    > **Allowed values:** (String) (optional) Any valid path pointing to an existing and valid PLINK *"`.freq`"* file, or *`None`* (`~` character).

  - #### <a name="min-bq-1"></a>min-BQ
    > **Description:**  Skip nucleotides carrying a base quality that is lower than the provided threshold. Note that `TKGWV2` internally calls samtools and passes this argument directly to it. The default recommended value of the software is 30

    > **Allowed values:** (integer) a non-negative integer in the range $[0, 40]$, representing a PHRED quality score. 
  - #### <a name="min-mq-2"></a>min-MQ
    > **Description:**  Skip aligned fragments carrying a mapping quality that is lower than the provided threshold. Note that `TKGWV2` internally calls samtools and passes this argument directly to it. The default recommended value of the software is 30

    > **Allowed values:** (integer) a non-negative integer in the range $[0, 40]$, representing a PHRED quality score. 
  - #### <a name="min-depth-2"></a>min-depth
    > **Description:** Specify the minimum number of SNPs allowed to estimate relatedness. Default value is 1.

    > **Allowed values:** (integer) Any integer value in the range $[1, +\infty]$ representing a number of overlapping snps.
