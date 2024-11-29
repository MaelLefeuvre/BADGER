# `badger-plot` plotting parameters reference

This page lists and describes all the graphical parameters that can be configured via the `badger-plots template` module, and the [YAML](https://yaml.org/spec/) file it produces.

> [!IMPORTANT]
> Required parameters are flagged by a $[\mathsf{\color{red}Required}]$ tag  

---

## <a name="input"></a>input
> <ins>**Description**</ins>: $[\mathsf{\color{red}Required}]$  Path to a user-generated `badger-plots` input definition file. This file can be generated using the `badger-plots input` module. See the dedicated section in the main README here : [Generating an input definition file for `badger-plots`](/README.md#01-generating-an-input-definition-file-for-badger-plots).

> <ins>**Allowed Values**</ins>: (String) Any path linking to an existing and valid `badger-plots` input definition file (in [YAML format](https://yaml.org/spec/)).

## <a name="pedigree-codes"></a>pedigree-codes
> <ins><ins>**Description**</ins></ins>: $[\mathsf{\color{red}Required}]$ Path pointing to a `badger` `pedigree-codes` input definition file. The  provided file should be the same as the one used during ***all*** of the archived simulations results pointed through [input](#input). i.e. the provided file should be identical to the one specified in your `config.yml` file when running `badger` (See the corresponding parameter: [codes](/doc/README.badger-config.md#codes)).

> <ins>**Allowed values**</ins>: (String) Any path linking to an existing and valid `badger`. `pedigree-codes` definition file. Detailled specifications on the format and purpose of pedigree-codes definition files can be found here: [pedigree-codes-files](README.ped-sim-config.md#pedigree-codes-files).

## <a name="output-dir"></a>output-dir
> <ins>**Description**</ins>: Path pointing to a directory, where output plots and data should be stored. If this directory is non-existent, `badger-plots` will attempt to create it for you, provided it has the correct filesystem permissions. 

> <ins>**Allowed values**</ins>: (String) Any path pointing to a valid directory.

## <a name="reorder"></a>reorder
> <ins>**Description**</ins>: Change the plotting order of the benchmarked kinship estimation methods. The provided list of method names ***must*** match those found within the provided [input` definition file](#input). Providing this argument with an empty list (`[]`) will preserve the order originally found within the input definition file.

> <ins>**Allowed values**</ins>: (List[Strings]) A comma separated list of kinship estimation method names.

> <ins>**Default**</ins>: `[]`

> <ins>**Examples**</ins>:  
> `reorder: []`  
> `reorder: ["TKGWV2", "KIN", "READv2", "READ"]`  

## <a name="rename"></a>rename
> <ins>**Description**</ins>: Change the plotting labels of the benchmarked kinship estimation methods. Note that the renaming of labels is performed *after* the reordering step (see [reorder](#reorder)). Thus, the provided list of method labels ***must*** match the order provided with the [reorder](#reorder) parameter. Providing this argument with an empty list (`[]`) will preserve the labels originally found within the input definition file

> <ins>**Allowed values**</ins>: (List[Strings]) A comma separated list of kinship estimation method names.

> <ins>**Default**</ins>: `[]`

> <ins>**Examples**</ins>:  
> `rename: ["tkgwv2", "KIN", "read-v2", "read-v1"]`  
> `rename: []`

## <a name="exclude"></a>exclude
> <ins>**Description**</ins>: Ignore and exclude the specified kinship estimation methods found within [input](#input) from plotting. Note that the exclusion of kinship estimation methods is performed *after* [reorder](#reorder) ***and*** [rename](#rename). Thus, the provided list of method names ***must*** match the labels provided with the [rename](#rename) parameter. Providing this argument with an empty list (`[]`) will plot all of the available results found in [input](#input).

> <ins>**Allowed values**</ins>: (List[Strings]) A comma separated list of kinship estimation method names.

> <ins>**Default**</ins>: `[]`

> <ins>**Examples**</ins>:  
> `exclude: ["read-v1", "tkgwv2"]`  
> `exclude: []`

## <a name="performance-plot"></a>performance-plot
Categorizes plotting parameters arguments for the summarized Ordinal classification index performance plot. 
- ### <a name="filename"></a>filename
  > <ins>**Description**</ins>: Base name of the output plot. The specified filename should be provided *without* any file extension. The corresponding files will be created within the directory specified through [output-dir](#output-dir).

  > <ins>**Allowed values**</ins>: (String) Any valid string representing a filename. Use of special characters is not recommended. (Prefer using the following character set: `[a-zA-Z0-9-_.]`)
  
  > <ins>**Default**</ins>: *`"OCI-performance-plot"`

- ### <a name="transpose"></a>transpose
  > <ins>**Description**</ins>: Transpose the plotting order of benchmarked methods and biological conditions. When the value is set to `no`, the grid of confusion matrices within the plot is ordered by arranging each kinship estimation method as a distinct row, and each biological condition as a distinct column. Setting this value to `yes` will instead arrange biological conditions on a distinct *row*, and methods on a distinct *column*.

  > <ins>**Allowed values**</ins>: (boolean) Any value parseable as a boolean, according to the YAML format specifications (See: [YAML boolean](https://yaml.org/spec/1.2.2/#10212-boolean), for more details.)

  > <ins>**Default**</ins>: `no`

  > <ins>**Examples**</ins>:  
  > `transpose: yes`  
  > `transpose: true`  
  > `transpose: false`  

- ### <a name="scatterplot_ratio"></a>scatterplot_ratio
  > <ins>**Description**</ins>: Specify the display ratio between the height grid of confusion matrices, and the summarized scatterplot of OCI values. Increasing this value will *increase* the relative heigth of the scatterplot.

  > <ins>**Allowed values**</ins>: (float) Any non negative floating point value in the range $[0, 1]$, representing a ratio.

  > <ins>**Default**</ins>: `0.2`

- ### <a name="horizontal_margin"></a>horizontal_margin
  > <ins>**Description**</ins>: Specify the display margin separating each *row* of the plot. 

  > <ins>**Allowed values**</ins>: (float) Any non negative floating point value in the range $[0, 1]$, representing a display ratio.

  > <ins>**Default**</ins>: `0.02`

- ### <a name="vertical_margin"></a>vertical_margin
  > <ins>**Description**</ins>: Specify the display margin separating each *column* of the plot

  > <ins>**Allowed values**</ins>: (float) Any non negative floating point value in the range $[0, 1]$, representing a display ratio.
    
  > <ins>**Default**</ins>: `0.02`

- ### <a name="axis_fontsize"></a>axis_fontsize
  > <ins>**Description**</ins>: Specify the font size of every axis label and annotation within the plot (in px.)

  > <ins>**Allowed values**</ins>: (integer) Any non negative integer in the range $[1, +\infty]$, representing a pixel size.

  > <ins>**Default**</ins>: `14`

- ### <a name="legend"></a>legend
  Categorizes plotting parameters related to the figure legend of the OCI-performance-plot.
  - #### <a name="size"></a>size
    > <ins>**Description**</ins>: Specify the font size of the legend labels (in px.)

    > <ins>**Allowed values**</ins>: (integer) Any non negative integer in the range $[1, +\infty]$, representing a pixel size.

    > <ins>**Default**</ins>: `10`

  - #### <a name="xpos"></a>xpos
    > <ins>**Description**</ins>: Specify the horizontal position of the figure legend.

    > <ins>**Allowed values**</ins>: (float) Any floating point value in the range $[-\infty, +\infty]$, representing a relative ratio, in relation to the left of the plotting surface.

    > <ins>**Default**</ins>: `-0.1`

  - #### <a name="ypos"></a>ypos
    > <ins>**Description**</ins>: Specify the vertical position of the figure legend.

    > <ins>**Allowed values**</ins>: (float) Any floating point value in the range $[-\infty, +\infty]$, representing a relative ratio, in relation to the bottom of the plotting surface.

    > <ins>**Default**</ins>: `-0.1`

- ### <a name="cm"></a>cm
  Categorizes plotting parameters related to the display of every individual confusion matrix found within the main OCI performance plot.
  - #### <a name="condense"></a>condense
    > <ins>**Description**</ins>: Hide columns containing no observations within the confusion matrix. i.e. hide prediction classes where no predictions were ever made by the method.

    > <ins>**Allowed values**</ins>: (boolean) Any value parseable as a boolean, according to the YAML format specifications (See: [YAML boolean](https://yaml.org/spec/1.2.2/#10212-boolean), for more details.)

    > <ins>**Default**</ins>: `no`

  - #### <a name="ratio"></a>ratio
    > <ins>**Description**</ins>: Specify whether the amount of observations within each class of the confusion matrix should be displayed as a within-class ratio.

    > <ins>**Allowed values**</ins>: (boolean) Any value parseable as a boolean, according to the YAML format specifications (See: [YAML boolean](https://yaml.org/spec/1.2.2/#10212-boolean), for more details.)

    > <ins>**Default**</ins>: `no`

  - #### <a name="colorscale"></a>colorscale
    > <ins>**Description**</ins>: 

    > <ins>**Allowed values**</ins>: (String) Any valid (ColorBrewer)[https://en.wikipedia.org/wiki/ColorBrewer], or [Viridis](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html) color palette label.

    > <ins>**Default**</ins>: `"Blues"`

  - #### <a name="ticklabels"></a>ticklabels
    > <ins>**Description**</ins>: Specify custom ticklabels for the rows and columns of every confusion matrix. Note that the provided list of labels must match the dimensions of the matrices. Setting this value to None (`~`) will instead use the expected relatedness coefficients found in the provided [`pedigree-codes`](#pedigree-codes) file.

    > <ins>**Allowed values**</ins>: (List[String]) A list of strings, representing labels for every row of a confusion matrix.

    > <ins>**Default**</ins>: `~`

    > <ins>**Examples**</ins>:  
    > `ticklabels: ~`  
    > `ticklabels: ["U", "3°", "2°", "1°", "S"]`  
    > `ticklabels: ["Unr", "Third", "Second", "First", "Self"]`  

  - #### <a name="tickfont"></a>tickfont
    > <ins>**Description**</ins>: Specify the font size of every confusion matrix' ticklabels. (in px.)

    > <ins>**Allowed values**</ins>: (integer) Any non negative value in the range $[1, +\infty]$, representing a pixel size.

    > <ins>**Default**</ins>: `10`

  - #### <a name="tickangle"></a>tickangle
    > <ins>**Description**</ins>: Specify the display angle of every confusion matrix' ticklabel (in degrees.)

    > <ins>**Allowed values**</ins>: (integer) Any integer value in the range $[0, 360]$, representing an angle. 

    > <ins>**Default**</ins>: `0`

  - #### <a name="fontsize"></a>fontsize
    > <ins>**Description**</ins>: Specify the size of every prediction count within the cells of the confusion matrices (in px.).

    > <ins>**Allowed values**</ins>: (integer) Any integer value in the range $[1, +\infty]$, representing a pixel size.

    > <ins>**Default**</ins>: `10`

  - #### <a name="show_xaxis"></a>show_xaxis
    > <ins>**Description**</ins>: Unilaterally choose to display or not the horizontal axes of the confusion matrices.

    > <ins>**Allowed values**</ins>: (boolean) Any value parseable as a boolean, according to the YAML format specifications (See: [YAML boolean](https://yaml.org/spec/1.2.2/#10212-boolean), for more details.)

    > <ins>**Default**</ins>: `no`

- ### <a name="scatter"></a>scatter
  Categorizes plotting parameters related to the display of the summary scatterplot of OCI performance values.
  - #### <a name="dash"></a>dash
    > <ins>**Description**</ins>: Specify a list of drawing style for every line found within the scatterplot. Note that values of the list will be recycled, if the number of lines to plot is greater than the number of elements within the list. 

    > <ins>**Allowed values**</ins>: (List[strings]) A list of plotly-compatible line dash styles. (i.e.: *`"solid"`*, *`"dot"`*, *`"dash"`*, *`"longdash"`*, *`"dashdot"`*, or *`"longdashdot"`*). For an extensive description of allowed values, see the corresponding plotly reference entry: [scatter line dash](https://plotly.com/r/reference/scatter/#scatter-line-dash)

    > <ins>**Default**</ins>: `"solid"`

    > <ins>**Examples**</ins>:  
    > `dash: "solid"`  
    > `dash: ["solid", "dash"]`  
    > `dash: ["solid", "dash", "dashdot", "dot"]`

  - #### <a name="mode"></a>mode
    > <ins>**Description**</ins>: Specify the general drawing mode of the scatterplot. 

    > <ins>**Allowed values**</ins>: (String) A plotly-compatible string, representing a drawing mode (i.e.: "lines", "markers", "line+markers"). See the corresponding plotly reference entry, for an extensive description of allowed values: [scatter mode](https://plotly.com/r/reference/scatter/#scatter-mode)

    > <ins>**Default**</ins>: *`"lines+markers"`*

  - #### <a name="yaxis"></a>yaxis
    Categorizes plotting parameters related to the display of the summary scatterplot's y-axis.

    - ##### <a name="range"></a>range
      > <ins>**Description**</ins>: Specify the display range of the y axis.

      > <ins>**Allowed values**</ins>: (List[float]) A list of two floating point values, representing the minimum and maximum displayed value of the yaxis, respectively.

      > <ins>**Default**</ins>: `[0.38, 1.02]`

    - ##### <a name="tickfont-1"></a>tickfont
      > <ins>**Description**</ins>: Specify the font size of the yaxis tick labels. (in px.)

      > <ins>**Allowed values**</ins>: (integer) Any non negative integer value in the range $[1, +\infty]$, representing a pixel size.

      > <ins>**Default**</ins>: `10`

    - ##### <a name="dtick"></a>dtick
      > <ins>**Description**</ins>: Specify the interval between displayed tick labels on the yaxis.

      > <ins>**Allowed values**</ins>: (float) Any floating point value in the range $[0, 1]$, representing a y-axis value.

      > <ins>**Default**</ins>: `0.1`

  - #### <a name="colors"></a>colors
    > <ins>**Description**</ins>: Specify a custom list of colors for every line of the scatterplot. If set to None (`~`), `badger.plots` will automatically generate a default color palette, using RColorBrewer's `Set2` palette.

    > <ins>**Allowed values**</ins>: (List<String>) A list of hexadecimal color codes, one for each plotted line.

    > <ins>**Default**</ins>: `~`

    > <ins>**Examples**</ins>:  
    > `colors: ~`  
    > `colors: ["#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"]`  
    > `colors: ["#00204D"]`  

- ### <a name="width"></a>width
  > <ins>**Description**</ins>: Specify the default width of the output plot (in px.)

  > <ins>**Allowed values**</ins>: (integer) Any non negative value in the range $[1, +\infty]$, representing a displayed width (in pixels).

  > <ins>**Default**</ins>: `1920`

- ### <a name="height"></a>height
  > <ins>**Description**</ins>: Specify the default height of the output plot (in px.)

  > <ins>**Allowed values**</ins>: (integer) Any non negative value in the range $[1, +\infty]$, representing a displayed height (in pixels).

  > <ins>**Default**</ins>: `1080`

## <a name="accuracy-plot"></a>accuracy-plot
Categorizes plotting parameters arguments for the normalized root mean square deviation (`nRMSD`) and mean bias error (`nMBE`) summary statistics plot. 
- ### <a name="filename-1"></a>filename
  > <ins>**Description**</ins>: Base name of the output plot. The specified filename should be provided *without* any file extension. The corresponding files will be created within the directory specified through [output-dir](#output-dir).

  > <ins>**Allowed values**</ins>: (String) Any valid string representing a filename. Use of special characters is not recommended. (Prefer using the following character set: `[a-zA-Z0-9-_.]`)
  
  > <ins>**Default**</ins>: *`"nRMSD-accuracy-plot"`


- ### <a name="transpose-1"></a>transpose

  > <ins>**Description**</ins>: Transpose the plotting order of benchmarked methods and biological conditions. When the value is set to `no`, every tested biological condition is plotted within a separate plot, while individual bars will correspond to a given kinship estimation method. Setting this value to `yes` will instead split the results of every kinship estimation method in a separate subplot, while individual bars will correspond to a given biological condition.

  > <ins>**Allowed values**</ins>: (boolean) Any value parseable as a boolean, according to the YAML format specifications (See: [YAML boolean](https://yaml.org/spec/1.2.2/#10212-boolean), for more details.)

  > <ins>**Default**</ins>: `no`

- ### <a name="fixed_axis"></a>fixed_axis
  > <ins>**Description**</ins>: Specify whether every `nRMSD` and `nMBE` subplot should share the same y-axis range. When set to `no`, the y-axis range of every subplot is tailored to the given values. When set to `yes`, the program will instead search for the maximum and minimum value across every subplot, and use these as a shared plotting range across all subplots.

  > <ins>**Allowed values**</ins>: (boolean) Any value parseable as a boolean, according to the YAML format specifications (See: [YAML boolean](https://yaml.org/spec/1.2.2/#10212-boolean), for more details.)

  > <ins>**Default**</ins>: `yes`

- ### <a name="title"></a>title
  > <ins>**Description**</ins>: Specify a main title for the plot.

  > <ins>**Allowed values**</ins>: (String) Any string representing a main title for the plot.

  > <ins>**Default**</ins>: `""`

- ### <a name="rmsd"></a>rmsd
  Categorizes plotting parameters shared across all `nRMSD` subplots
  - #### <a name="title-1"></a>title
    > <ins>**Description**</ins>: Specify a subtitle for the `nRMSD` subplots

    > <ins>**Allowed values**</ins>: (String) Any string representing a main subtitle for the `nRMSD` subplots.

    > <ins>**Default**</ins>: `nRMSD`

  - #### <a name="dtick-1"></a>dtick
    > <ins>**Description**</ins>: Specify the interval between displayed tick labels on the yaxis.

    > <ins>**Allowed values**</ins>: (float) Any floating point value in the range $[0, 1]$, representing a y-axis value.

    > <ins>**Default**</ins>: `0.2`

  - #### <a name="tickangle-1"></a>tickangle
    > <ins>**Description**</ins>: Specify the display angle of every displayed ticklabel (in degrees.)

    > <ins>**Allowed values**</ins>: (integer) Any integer value in the range $[0, 360]$, representing an angle. 

    > <ins>**Default**</ins>: `0`

- ### <a name="mbe"></a>mbe
  Categorizes plotting parameters shared across all `nMBE` subplots

  - #### <a name="title-2"></a>title
    > <ins>**Description**</ins>: Specify a subtitle for the `nMBE` subplots

    > <ins>**Allowed values**</ins>: (String) Any string representing a main subtitle for the `nMBE` subplots.

    > **Default:** `nMBE`
  
  - #### <a name="dtick-2"></a>dtick
    > <ins>**Description**</ins>: Specify the interval between displayed tick labels on the yaxis.

    > <ins>**Allowed values**</ins>: (float) Any floating point value in the range $[0, 1]$, representing a y-axis value.

    > <ins>**Default**</ins>: `0.2`

  - #### <a name="tickangle-2"></a>tickangle
    > <ins>**Description**</ins>: Specify the display angle of every displayed ticklabel (in degrees.)

    > <ins>**Allowed values**</ins>: (integer) Any integer value in the range $[0, 360]$, representing an angle. 

    > <ins>**Default**</ins>: `45`

  - #### <a name="scale_factor"></a>scale_factor
    > <ins>**Description**</ins>: Multiply the displayed values of `nMBE` by the provided scaling factor. Note that setting this value to anything other than `1` will force the program to display the scaling factor in the corresponding `nMBE` [title](#title-1).

    > <ins>**Allowed values**</ins>: (float) Any floating point value in the range $[-\infty, +\infty]$, representing a scaling factor for the y-axis values.

    > <ins>**Default**</ins>: `1`

- ### <a name="mbe_plot_ratio"></a>mbe_plot_ratio
  > <ins>**Description**</ins>: Specify the display ratio between the column of `nMBE` subplots and the column of `nRMSD` subplots. Note that increasing this value will have the effect of increasing the display width of every `nMBE` subplot.

  > <ins>**Allowed values**</ins>: (float) Any floating point value in the range $[0, 1]$, representing a display width ratio. 

  > <ins>**Default**</ins>: `0.3`

- ### <a name="ticksize"></a>ticksize
  Categorizes plotting parameters related to the display size of tick labels within the accuracy plot
  - #### <a name="xaxis"></a>xaxis
    > <ins>**Description**</ins>: Specify the font size of every label found on the x-axis (in px.).

    > <ins>**Allowed values**</ins>: (integer) Any non negative integer value in the range $[1, +\infty]$, representing a pixel size.

    > <ins>**Default**</ins>: `16`

  - #### <a name="yaxis-1"></a>yaxis
    > <ins>**Description**</ins>: Specify the font size of every label found on the y-axis (in px.).

    > <ins>**Allowed values**</ins>: (integer) Any non negative integer value in the range $[1, +\infty]$, representing a pixel size.

    > **Default:** `12`

- ### <a name="legend-1"></a>legend
  Categorizes plotting parameters related to the display of the legend within the main accuracy plot.
  - #### <a name="size-1"></a>size
    > <ins>**Description**</ins>: Specify the font size of every figure legend key. (in px.)

    > <ins>**Allowed values**</ins>: (integer) Any non negative integer value in the range $[1, +\infty]$, representing a pixel size.

    > <ins>**Default**</ins>: `12`

  - #### <a name="title-3"></a>title
    > <ins>**Description**</ins>: Specify a title categorizing every item within the legend

    > <ins>**Allowed values**</ins>: (string) Any string representing a legend title.

    > <ins>**Default**</ins>: `"Method"`

- ### <a name="marker"></a>marker
  Categorizes plotting parameters related to the display of every marker displayed within the plot.
  - #### <a name="colors-1"></a>colors
    > <ins>**Description**</ins>: Specify a custom list of colors for every grouped bar within the subplots. If set to None (`~`), `badger.plots` will automatically generate a default color palette. Depending on the value of [transpose](#transpose-1), the default color palette will either be RColorBrewer's `Set2`, or the `viridis` color palette.

    > <ins>**Allowed values**</ins>: (List<String>) A list of hexadecimal color codes, one for each grouped bar.

    > <ins>**Default**</ins>: `~`

    > <ins>**Examples**</ins>:  
    > `colors: ~`  
    > `colors: ["#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"]`  
    > `colors: ["#00204D"]`  

  - #### <a name="patterns"></a>patterns
    Categorizes arguments related to the fill pattern of bars within the plot. See the corresponding entry in plotly's reference: [bar marker pattern](https://plotly.com/r/reference/bar/#bar-marker-pattern)
    - ##### <a name="shape"></a>shape
      > <ins>**Description**</ins>: Specify the shape(s) of the hatching-pattern fill of every grouped bar.  Note that values of the list will be recycled, if the number of grouped bars is greater than the number of elements within the list.

      > <ins>**Allowed values**</ins>: (List[strings]) A list of plotly-compatible fill pattern labels (i.e.: *`""`*, *`"/"`*, *`"\"`*, *`"x"`*, *`"-"`*, *`"|"`*, *`"+"`*, *`"."`*), where each symbol corresponds to a different hatch-pattern. For an extensive description of allowed values, see the corresponding plotly reference entry: [pattern shape](https://plotly.com/r/reference/bar/#bar-marker-pattern-shape)

      > <ins>**Default**</ins>: `""`

      > <ins>**Examples**</ins>:  
      > `colors: ~`  
      > `colors: ["", "/"]`  
      > `colors: ["", "x", "/" "+"]`  

    - ##### <a name="solidity"></a>solidity
      > <ins>**Description**</ins>: Specify the fraction of the area filled by the provided patterns. Increasing this value will *increase* the filled area of the hatching pattern.

      > <ins>**Allowed values**</ins>: (float) Any floating point value in the range $[0, 1]$, representing a fraction of filled area.

      > <ins>**Default**</ins>: `0.5`

    - ##### <a name="size-2"></a>size
      > <ins>**Description**</ins>: Specify the size of the shapes used to generate the hatching pattern (in px.)

      > <ins>**Allowed values**</ins>: (integer) Any integer value in the range $[1, +\infty]$, representing a pixel size. 

      > <ins>**Default**</ins>: `5`

- ### <a name="width-1"></a>width
  > <ins>**Description**</ins>: Specify the default width of the output plot (in px.)

  > <ins>**Allowed values**</ins>: (integer) Any non negative value in the range $[1, +\infty]$, representing a displayed width (in pixels).

  > <ins>**Default**</ins>: `1120`

- ### <a name="height-1"></a>height
  > <ins>**Description**</ins>: Specify the default height of the output plot (in px.)

  > <ins>**Allowed values**</ins>: (integer) Any non negative value in the range $[1, +\infty]$, represending a displayed height (in pixels).

  > <ins>**Default**</ins>: `1190`
