# `badger-plot` plotting parameters reference

This page lists and describes all the graphical parameters that can be configured via the `badger-plots template` module, and the [YAML](https://yaml.org/spec/) file it produces.

Note that Required parameters are flagged by a '$\mathsf{[\color{red}Required\color{reset}]}$' tag

--- 

- ## $\mathsf{[\color{red}Required\color{reset}]}$ input
  > **Description**:  Path to a user-generated `badger-plots` input definition file. This file can be generated using the `badger-plots input` module. See the dedicated section in the main README here : [Generating an input definition file for `badger-plots`](/README.md#01-generating-an-input-definition-file-for-badger-plots).

  > **Allowed Values**: (String) Any path linking to an existing and valid `badger-plots` input definition file (in [YAML format](https://yaml.org/spec/)).

- ## $[\mathsf{\color{red}Required}]$ pedigree-codes
  > **Description**: Path pointing to a `badger` `pedigree-codes` input definition file. The  provided file should be the the same as the one used during ***all*** of the archived simulations results pointed through [input](#input). i.e. The provided file should be identical to the one specified in your `config.yml` file when running `badger` (See the corresponding parameter: [codes](/doc/README.badger-config.md#codes)).

  > **Allowed Values**: (String) Any path linking to an existing and valid `badger`. `pedigree-codes` definition file. Detailled specifications on the format and purpose of pedigree-codes definition files can be found here: [pedigree-codes-files](README.ped-sim-config.md#pedigree-codes-files).

- ## output-dir
  > **Description**: Path pointing to a directory, where output plots and data should be stored. If this directory is non existent, `badger-plots` will attempt to create it for you, provided it has the correct filesystem permissions. 

  > **Allowed Values**: (String) Any path pointing to a valid directory.

- ## reorder
  > **Description**: Change the plotting order of the benchmarked kinship estimation methods. The provided list of method names ***must*** match those found within the provided [input` definition file](#input). Providing this argument with an empty list (`[]`) will preserve the order originally found within the input definition file.

  > **Allowed Values**: (List[Strings]) A comma separated list of kinship estimation method names.

  > **Default**: `[]`

  > **Examples**:  
  > `reorder: []`  
  > `reorder: ["TKGWV2", "KIN", "READv2", "READ"]`  

- ## rename
  > **Description**: Change the plotting labels of the benchmarked kinship estimation methods. Note that the renaming of labels is performed *after* the reordering step (see [reorder](#reorder)). Thus, the provided list of method labels ***must*** match the order provided with the [reorder](#reorder) parameter. Providing this argument with an empty list (`[]`) will preserve the labels originally found within the input definition file

  > **Allowed Values**: (List[Strings]) A comma separated list of kinship estimation method names.

  > **Default**: `[]`

  > **Examples**:  
  > `rename: ["tkgwv2", "KIN", "read-v2", "read-v1"]`  
  > `rename: []`

- ## exclude
  > **Description**: Ignore and exclude the specified kinship estimation methods found within [input](#input) from plotting. Note that the exclusion of kinship estimation methods is performed *after* [reorder](#reorder) ***and*** [rename](#rename). Thus, the provided list of method names ***must*** match the labels provided with the [rename](#rename) parameter. Providing this argument with an empty list (`[]`) will plot all of the available results found in [input](#input).

  > **Allowed Values**: (List[Strings]) A comma separated list of kinship estimation method names.

  > **Default**: `[]`

  > **Examples**:  
  > `exclude: ["read-v1", "tkgwv2"]`  
  > `exclude: []`

- ## performance-plot
  Categorizes plotting parameters arguments for the summarized Ordinal classification index performance plot. 
  - ### filename
    > **Description**: Base name of the output plot. The specified filename should be provided *without* any file extension. The corresponding files will be created within the directory specified through [output-dir](#output-dir).

    > **Allowed Values**: (String) Any valid string representing a filename. Use of special characters is not recommended. (Prefer using the following character set: `[a-zA-Z0-9-_.]`)
  
    > **Default**: *`"OCI-performance-plot"`

  - ### transpose
    > **Description**: Transpose the plotting order of benchmarked methods and biological conditions. When the value is set to `no`, The grid of confusion matrices within the plot is ordered by arranging each kinship estimation method as a distinct row, and each biological condition as a distinct column. Setting this value to `yes` will instead arrange biological conditions on a distinct *row*, and methods on a distinct *column*.

    > **Allowed Values**: (boolean) Any value parseable as a boolean, according to the YAML format specifications (See: [YAML boolean](https://yaml.org/spec/1.2.2/#10212-boolean), for more details.)

    > **Default**: `no`

    > **Examples**:  
    > `transpose: yes`  
    > `transpose: true`  
    > `transpose: false`  

  - ### scatterplot_ratio
    > **Description**: Specify the display ratio between the height grid of confusion matrices, and the summarized scatterplot of OCI values. Increasing this value will *increase* the relative heigth of the scatterplot.

    > **Allowed Values**: (float) Any non negative floating point value in the range $[0, 1]$, representing a ratio.

    > **Default**: `0.2`

  - ### horizontal_margin
    > **Description**: Specify the display margin separating each *row* of the plot. 

    > **Allowed Values**: (float) Any non negative floating point value in the range $[0, 1]$, representing a display ratio.

    > **Default**: `0.02`

  - ### vertical_margin
    > **Description**: Specify the display margin separating each *column* of the plot

    > **Allowed Values**: (float) Any non negative floating point value in the range $[0, 1]$, representing a display ratio.
    
    > **Default**: `0.02`

  - ### axis_fontsize
    > **Description**: Specify the font size of every axis label and annotation within the plot (in px.)

    > **Allowed Values**: (integer) Any non negative integer in the range $[1, +\infty]$, representing a pixel size.

    > **Default**: `14`

  - ### legend:
    Categorizes plotting parameters related to the figure legend of the OCI-performance-plot.
    - #### size: 10
      > **Description**: Specify the font size of the legend labels (in px.)

      > **Allowed Values**: (integer) Any non negative integer in the range $[1, +\infty]$, representing a pixel size.

    - #### xpos
      > **Description**: Specify the horizontal position of the figure legend.

      > **Allowed Values**: (float) Any floating point value in the range $[-\infty, +\infty]$, representing a relative ratio, in relation to the left of the plotting surface.

      > **Default**: `-0.1`

    - #### ypos
      > **Description**: Specify the vertical position of the figure legend.

      > **Allowed Values**: (float) Any floating point value in the range $[-\infty, +\infty]$, representing a relative ratio, in relation to the bottom of the plotting surface.

      > **Default**: `-0.1`

  - ### cm:
      Categorizes plotting parameters related to the display of every individual confusion matrix found within the main OCI performance plot.
    - #### condense
      > **Description**: Hide columns containing no observations within the confusion matrix. i.e. hide prediction classes where no predictions were ever made by the method.

      > **Allowed Values**: (boolean) Any value parseable as a boolean, according to the YAML format specifications (See: [YAML boolean](https://yaml.org/spec/1.2.2/#10212-boolean), for more details.)

      > **Default**: `no`

    - #### ratio:
      > **Description**: Specify whether the amount of observations within each class of the confusion matrix should be displayed as a within-class ratio.

      > **Allowed Values**: (boolean) Any value parseable as a boolean, according to the YAML format specifications (See: [YAML boolean](https://yaml.org/spec/1.2.2/#10212-boolean), for more details.)

      > **Default**: `no`

    - #### colorscale
      > **Description**: 

      > **Allowed Values**: (String) Any valid (ColorBrewer)[https://en.wikipedia.org/wiki/ColorBrewer], or [Viridis](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html) color palette label.

      > **Default**: `"Blues"`

    - #### ticklabels
      > **Description**: Specify custom ticklabels for the rows and columns of every confusion matrix. Note that the provided list of labels must match the dimensions of the matrices. Setting this value to None (`~`) will instead use the expected relatedness coefficients found in the provided [`pedigree-codes`](#pedigree-codes) file.

      > **Allowed Values**: (List[String]) A list of strings, representing labels for every row of a confusion matrix. 
      > **Default**: `~`

      > **Examples**:  
      > `ticklabels: ~`  
      > `ticklabels: ["U", "3°", "2°", "1°", "S"]`  
      > `ticklabels: ["U", "3°", "2°", "1°", "S"]`  
      > `ticklabels: ["Unr", "Third", "Second", "First", "Self"]`  

    - #### tickfont
      > **Description**: Specify the fontsize of every confusion matrix' ticklabels. (in px.)

      > **Allowed Values**: (integer) Any non negative value in the range $[1, +\infty]$, representing a pixel size.

      > **Default**: `10`

    - #### tickangle
      > **Description**: Specify the display angle of every confusion matrix' ticklabel (in degrees.)

      > **Allowed Values**: (integer) Any integer value in the range $[0, 360]$, representing an angle. 

      > **Default**: `0`

    - #### fontsize
      > **Description**: Specify the size of every prediction count within the cells of the confusion matrices (in px.).

      > **Allowed Values**: (integer) Any integer value in the range $[1, +\infty]$, representing a pixel size.

      > **Default**: `10`

    - #### show_xaxis: no
      > **Description**: Unilaterally choose to display or not the horizontal axes of the confusion matrices.

      > **Allowed Values**: (boolean) Any value parseable as a boolean, according to the YAML format specifications (See: [YAML boolean](https://yaml.org/spec/1.2.2/#10212-boolean), for more details.)

      > **Default**: `no`

  - ### scatter:
      Categorizes plotting parameters related to the display of the summary scatterplot of OCI performance values.
    - #### dash: solid
      > **Description**: Specify A list of drawing style for every line found within the scatterplot. Note that values of the list will be recycled, if the number of lines to plot is greater than the number of elements within the list. 

      > **Allowed Values**: (List[strings]) A list of plotly-compatible line dash styles. (i.e.: *`"solid"`*, *`"dot"`*, *`"dash"`*, *`"longdash"`*, *`"dashdot"`*, or *`"longdashdot"`*). For an extensive description of allowed values, see the corresponding plotly reference entry: [scatter line dash](https://plotly.com/r/reference/scatter/#scatter-line-dash)

      > **Default**: `"solid"`

      > **Examples**:  
      > `dash: "solid"`  
      > `dash: ["solid", "dash"]`  
      > `dash: ["solit", "dash", "dashdot", "dot"]`

    - #### mode
      > **Description**: Specify the general drawing mode of the scatterplot. 

      > **Allowed Values**: (String) A plotly-compatible string, representing a drawing mode (i.e.: "lines", "markers", "line+markers"). See the corresponding plotly reference entry, for an extensive description of allowed values: [scatter mode](https://plotly.com/r/reference/scatter/#scatter-mode)

      > **Default**: *`"lines+markers"`*

    - #### yaxis:
      Categorizes plotting parameters related to the display of the summary scatterplot's y-axis.

      - ##### range:
        > **Description**: Specify the display range of the y axis.

        > **Allowed Values**: (List[float]) A list of two floating point values, representing the minimum and maximum displayed value of the yaxis, respectively.

        > **Default**: `[0.38, 1.02]`

      - ##### tickfont
        > **Description**: Specify the fontsize of the yaxis tick labels. (in px.)

        > **Allowed Values**: (integer) Any non negative integer value in the range $[1, +\infty]$, representing a pixel size.

        > **Default**: `10`

      - ##### dtick
        > **Description**: Specify the interval between every displayed tick labels on the yaxis.

        > **Allowed Values**: (float) Any floating point value in the range $[0, 1]$, representing a y-axis value.

        > **Default**: `0.1`

    - #### colors
      > **Description**: Specify a custom list of colors for every line of the scatterplot. If set to None (`~`), `badger.plots` will automatically generate a default color palette, using RColorBrewer's `Set2` palette.

      > **Allowed Values**: (List<String>) A list of hexadecimal color codes, one for each plotted line.

      > **Default**: `~`

      > **Examples**:  
      > `colors: ~`  
      > `colors: ["#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"]`  
      > `colors: ["#00204DFF]`  

  - ### width
    > **Description**: Specify the default width of the output plot (in px.)

    > **Allowed Values**: (integer) Any non negative value in the range $[1, +\infty]$, represending a displayed width (in pixels).

    > **Default**: `1920`

  - ### height
    > **Description**: Specify the default height of the output plot (in px.)

    > **Allowed Values**: (integer) Any non negative value in the range $[1, +\infty]$, represending a displayed width (in pixels).

    > **Default**: `1080`

- ## accuracy-plot:
  - ### filename
    > **Description**: Base name of the output plot. The specified filename should be provided *without* any file extension. The corresponding files will be created within the directory specified through [output-dir](#output-dir).

    > **Allowed Values**: (String) Any valid string representing a filename. Use of special characters is not recommended. (Prefer using the following character set: `[a-zA-Z0-9-_.]`)
  
    > **Default**: *`"nRMSD-accuracy-plot"`


  - ### transpose

    > **Description**: Transpose the plotting order of benchmarked methods and biological conditions. When the value is set to `no`, Every tested biological condition is plotted within a separate plot, while individual bars will correspond to a given kinship estimation method. Setting this value to `yes` will instead split the results of every kinship estimation method in a separate subplot, while individual bars will correspond to a given biological condition.

    > **Allowed Values**: (boolean) Any value parseable as a boolean, according to the YAML format specifications (See: [YAML boolean](https://yaml.org/spec/1.2.2/#10212-boolean), for more details.)

    > **Default**: `no`

  - ### fixed_axis
    > **Description**: Specify whether every `nRMSD` and `nMBE` subplot should share the same y-axis range. When set to `no`, the y-axis range of every subplot is tailored to the given values. When set to `yes`, the program will instead search for the maximum and minimum value across every subplot, and use these as a shared plotting range across all subplots.

    > **Allowed Values**: (boolean) Any value parseable as a boolean, according to the YAML format specifications (See: [YAML boolean](https://yaml.org/spec/1.2.2/#10212-boolean), for more details.)

    > **Default**: `yes`

  - ### title
    > **Description**: Specify a main title for the plot.

    > **Allowed Values**: (String) Any string representing a main title for the plot.

    > **Default**: `""`

  - ### rmsd:
    Categorizes plotting parameters shared across every `nRMSD` subplots
    - #### title: nRMSD
      > **Description**: Specify a subtitle for the `nRMSD` subplots

      > **Allowed Values**: (String) Any string representing a main subtitle for the `nRMSD` subplots.

      > **Default**: `nRMSD`

    - #### dtick: 0.2
        > **Description**: Specify the interval between every displayed tick labels on the yaxis.

        > **Allowed Values**: (float) Any floating point value in the range $[0, 1]$, representing a y-axis value.

        > **Default**: `0.2`

    - #### tickangle
      > **Description**: Specify the display angle of every displayed ticklabel (in degrees.)

      > **Allowed Values**: (integer) Any integer value in the range $[0, 360]$, representing an angle. 

      > **Default**: `0`

  - ### mbe:
    Categorizes plotting parameters shared across every `nMBE` subplots

    - #### title
      > **Description**: Specify a subtitle for the `nMBE` subplots

      > **Allowed Values**: (String) Any string representing a main subtitle for the `nMBE` subplots.

      > **Default:** `nMBE`
    
    - #### dtick
        > **Description**: Specify the interval between every displayed tick labels on the yaxis.

        > **Allowed Values**: (float) Any floating point value in the range $[0, 1]$, representing a y-axis value.

        > **Default**: `0.2`

    - #### tickangle
      > **Description**: Specify the display angle of every displayed ticklabel (in degrees.)

      > **Allowed Values**: (integer) Any integer value in the range $[0, 360]$, representing an angle. 

      > **Default**: `45`

    - #### scale_factor
      > **Description**: Multiply the displayed values of `nMBE` by the provided scaling factor. Note that setting this value to anything other than `1` will force the program to display the scaling factor in the corresponding `nMBE` [title](#title-1).

      > **Allowed Values**: (float) Any floating point value in the range $[-\infty, +\infty]$, representing a scaling factor for the y-axis values.

      > **Default**: `1`

  - ### mbe_plot_ratio
    > **Description**: Specify the display ratio between the column of `nMBE` subplots and the column of `nRMSD` subplots. Note that increasing this value will have the effect of increasing the display width of every `nMBE` subplot.

    > **Allowed Values**: (float) Any floating point value in the range $[0, 1]$, reprenting a display width ratio. 

    > **Default**: `0.3`

  - ### ticksize:
    Categorizes plotting parameters related to the display size of tick labels within the accuracy plot
    - #### xaxis
      > **Description**: Specify the fontsize of every label found on the x-axis (in px.).

      > **Allowed Values**: (integer) Any non negative integer value in the range $[1, +\infty]$, reprenting a pixel size.

      > **Default**: `16`

    - #### yaxis
      > **Description**: Specify the fontsize of every label found on the y-axis (in px.).

      > **Allowed Values**: (integer) Any non negative integer value in the range $[1, +\infty]$, reprenting a pixel size.

      > **Default:** `12`

  - ### legend
    Categorizes plotting parameters related to the display of the legend within the main accuracy plot.
    - #### size
      > **Description**: Specify the fontsize of every figure legend key. (in px.)

      > **Allowed Values**: (integer) Any non negative integer value in the range $[1, +\infty]$, reprenting a pixel size.

      > **Default**: `12`

    - #### title
      > **Description**: Specify a title categorizing every item within the legend

      > **Allowed Values**: (string) Any string representing a legend title.

      > **Default**: `"Method"`

  - ### marker
    Categorizes plotting parameters related to the display of every marker displayed within the plot.
    - #### colors
      > **Description**: Specify a custom list of colors for every grouped bar within the subplots. If set to None (`~`), `badger.plots` will automatically generate a default color palette. Depending on the value of [transpose](#transpose-1), the default color palette will either be RColorBrewer's `Set2`, or the `viridis` color palette.

      > **Allowed Values**: (List<String>) A list of hexadecimal color codes, one for each grouped bar.

      > **Default**: `~`

      > **Examples**:  
      > `colors: ~`  
      > `colors: ["#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"]`  
      > `colors: ["#00204DFF]`  

    - #### patterns:
      Categorizes arguments related to the fill pattern of bars within the plot. See the corresponding entry in plotly's reference: [bar marker pattern](https://plotly.com/r/reference/bar/#bar-marker-pattern)
      - ##### shape
        > **Description**: Specify the shape(s) of the hatching-pattern fill of every grouped.  Note that values of the list will be recycled, if the number of grouped bars is greater than the number of elements within the list.

        > **Allowed Values**: (List[strings]) A list of plotly-compatible fill pattern labels (i.e.: *`""`*, *`"/"`*, *`"\"`*, *`"x"`*, *`"-"`*, *`"|"`*, *`"+"`*, *`"."`*), where each symbol corresponds to a different hatch-pattern. For an extensive description of allowed values, see the corresponding plotly reference entry: [pattern shape](https://plotly.com/r/reference/bar/#bar-marker-pattern-shape)

        > **Default**: `""`

        > **Examples**:  
        > `colors: ~`  
        > `colors: ["", "/"]`  
        > `colors: ["", "x", "/" "+"]`  

      - ##### solidity
        > **Description**: Specify the fraction of the area filled by the provided patterns. Increasing this value will *increase* the filled area of the hatching pattern.

        > **Allowed Values**: (float) Any floating point value in the range $[0, 1]$, reprenting a fraction of filled area.

        > **Default**: `0.5`

      - ##### size
        > **Description**: Specify the size of the shapes used to generate the hatching pattern (in px.)

        > **Allowed Values**: (integer) Any integer value in the range $[1, +\infty]$, representing a pixel size. 

        > **Default**: `5`

  - ### width
    > **Description**: Specify the default width of the output plot (in px.)

    > **Allowed Values**: (integer) Any non negative value in the range $[1, +\infty]$, represending a displayed width (in pixels).

    > **Default**: `1120`

  - ### height
    > **Description**: Specify the default height of the output plot (in px.)

    > **Allowed Values**: (integer) Any non negative value in the range $[1, +\infty]$, represending a displayed width (in pixels).

    > **Default**: `1190`
