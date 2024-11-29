# Creating custom pedigree configuration files for BADGER

## Pedigree definition files (.def) files

Pedigree definition files are provided through the `ped-sim['data']['definition']` keyword, within BADGER's main yaml configuration's file.

These files generally adhere to the proprietary format specifications of `ped-sim-1.4.2`, and extensive documentation on how to design these files may be found in the original documentation of `ped-sim` here: [def files](https://github.com/williamslab/ped-sim/blob/v1.4.2/README.md#def-file).

Here, BADGER currently comes with an additional restriction, in that the first three fields of the definition file (i.e., the `def`, `[name]`, and `[#copies]` field described in `ped-sim's` documentation) should ***always*** be written as follows:

```text
def ped {{ params.replicates }}
```
This current limitation resides in the facts that:
1. Snakemake makes extensive use of pattern matching to resolve job dependencies, explaining the strict requirement on setting the `[name]` field as `ped`
2. BADGER applies `jinja2` templating on these files to dynamically adjust the required number of pedigree replicates. Hence, the `{{ params.replicate }}` string for the `[#copies]` field.

Apart from this single limitation, the rest of the file strictly adheres to the format specification described in the documentation of `ped-sim`. 

Hence, the following pedigree definition file

```text
#def [name] [#copies] [#generations]
def ped {{ params.replicates }} 3
1 1 2
2 1 2 1:1_2 2:1_2
3 1 2 1:1 2:2
```

Will generate the following template pedigree within BADGER's simulations

<p align="center">

  | ![BADGER-pedigree-definition](assets/BADGER-pedigree-definition-01.svg) |
  |:---:|
  | ***Figure 1***: Example output of a simple pedigree definition file |

</p>


## Pedigree codes files

pedigree_codes.txt are provided through the `ped-sim['data']['codes']`keyword, within BADGER's main yaml configuration file.

The purpose of this file is simply to subset the number of investigated ties within the pedigree, and thus specifically target which relationships should be used by BADGER to apply its benchmark.
This file is a simple, unheaded and tab-separated txt file, containing four fields:

1. `<NAME>`: User-defined relationship label (e.g, *"Unrelated"*, *"Parent-Offspring"*, *"Cousins"*). Any valid alphanumeric string may be used here.
2. `ID1`   : `ped-sim` identifier of the first individual involved in the targeted relationship.
3. `ID2`   : `ped-sim` identifier of the *second* individual involved in the targeted relationship.
4. `Expected-r`: expected relatedness coefficient for this particular pair of individuals.

In other terms while a `pedigree.def` defines the nodes and topology of the pedigree, a `pedigree-codes.txt` file will define the edges. Hence, the following example pedigree-code file:

```text
Unrelated   g1-b1-i1    g1-b2-i1    0
Siblings    g2-b1-i1    g2-b2-i1    0.5
Avuncular   g3-b1-i1    g2-b2-i1    0.25
Cousins     g3-b1-i1    g3-b2-i1    0.125
```

Will request BADGER to target the following relationships when assessing classification performance:

<p align="center">

  | ![BADGER-pedigree-definition](assets/BADGER-pedigree-definition-02.svg) |
  |:---:|
  | ***Figure 2***: Example output of a pedigree-codes file |

</p>


## Simulating monozygotic twins within the pedigree.

BADGER currently provides with a simple way to generate monozygotic twins during the simulations. This can be useful to estimate the impact of such relationships on the accuracy and bias of kinship estimation methods.  

This feature can be requested by simply specifying a self-comparison within the `pedigree-codes` file. To distinguish original samples from their "doppelganger" counterparts, the `ID2` field ***must*** however be capitalized. 

Thus, providing BADGER with the following alternative `pedigree-codes.txt` file:

<pre><code>Unrelated   g1-b1-i1    g1-b2-i1    0 
Siblings    g2-b1-i1    g2-b2-i1    0.5
Avuncular   g3-b1-i1    g2-b2-i1    0.25
Cousins     g3-b1-i1    g3-b2-i1    0.125
<span style="color: gold"><strong><em>Twins       g3-b2-i1    G3-B2-I1    1</em></strong></span></code></pre>

Will both modify the tree topology, and investigate relationships as follows:
<p align="center">

  | ![BADGER-pedigree-definition](assets/BADGER-pedigree-definition-03.svg) |
  |:---:|
  | ***Figure 3***: Impact of specifying twins within a `pedigree-codes` file.|

</p>