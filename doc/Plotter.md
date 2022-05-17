# Operation Details

## Parameters
For an explanation of Plotter's parameters in the Run Config and by commandline, see [the parameters documentation](Parameters.md#plotter).

## Input Files
tinyRNA automatically handles file inputs when Plotter is called as a step in a [pipeline run](Pipeline.md). It is necessary to identify these files when running Plotter as a standalone step. The list of possible input files is lengthy, so you are free to specify only the subset of inputs that are required for your desired plot types. The dependencies per plot type are as follows:

| Plot Type                       | Input File Commandline Argument(s)                         | Source              |
|---------------------------------|------------------------------------------------------------|---------------------|
| len_dist                        | `--len-dist FILE FILE FILE ...`                            | Counter             |
| rule_charts                     | `--rule-counts FILE`                                       | Counter             |
| class_charts                    | `--raw-counts FILE`</br>`--summary-stats FILE`             | Counter</br>Counter |
| replicate_scatter               | `--norm-counts FILE`                                       | DESeq2              |
| sample_avg_scatter_by_dge       | `--norm-counts FILE`</br>`--dge-tables FILE FILE FILE ...` | DESeq2</br>DESeq2   |
| sample_avg_scatter_by_dge_class | `--norm-counts FILE`</br>`--dge-tables FILE FILE FILE ...` | DESeq2</br>DESeq2   |

# Plot Types

## len_dist
The distributions of 5' end nucleotides vs. sequence lengths can be used to assess the overall quality of your libraries. This can also be used for analyzing small RNA distributions in non-model organisms without annotations.

<p float="left" align="center">
    <img src="../images/plots/len_dist_short.jpeg" width="40%" alt="len_dist lengths 20-30"/>
    <img src="../images/plots/len_dist_long.jpeg" width="40%" alt="len_dist lengths 15-60"/>
</p>


#### Subtypes
Two plots are produced for each replicate:
- Distribution of _Mapped Reads_, which are counted for every alignment reported in Counter's input SAM files
- Distribution of _Assigned Reads_, which are counted at each alignment where at least one overlapping feature passed selection and was assigned a portion of the sequence's original counts

#### Length Bounds
Lengths are plotted over a continuous range, even if an intermediate length was not observed, and the bounds of this range can be assigned automatically or manually. Manual lengths can be assigned using [plot_len_dist_min and plot_len_dist_max](Parameters.md#bounds-for-len_dist-charts).

When Plotter is called as a step in a pipeline run, min and max bounds are determined independently in the following order of priority:
1. Manual assignment by `plot_len_dist_min` and / or `plot_len_dist_max` in the Run Config
2. The corresponding _optional_ entries for fastp (`length_required` and `length_limit`) in the Run Config
3. Automatic assignment from the data. Bounds are determined by considering the min/max lengths across all libraries such that all plots have the same bounds. This determination is performed separately for each plot subtype.

#### Non-Nucleotide Bases
Placeholder bases, e.g. N, will be reported if they are encountered at the 5' end. Otherwise only the 4 standard bases are reported.



## rule_charts
Counts are assigned only to the features that meet selection criteria at each alignment locus. It is useful to see how each selection rule contributed to the overall assignment of counts. The rule_charts plot type shows the percentage of _mapped reads_ that each rule contributed to the total _assigned reads_.

<p float="left" align="center">
    <img src="../images/plots/rule_charts.jpeg" width="60%" alt="rule_chart with 10 rules"/>
</p>

#### Rule Number
Rules are referred to by their row number in the Features Sheet and the first non-header row is considered rule 0. Rule **N** represents the percentage of mapped reads that were unassigned. Sources of unassigned reads include:
- A lack of features passing selection at alignment loci
- Alignments which do not overlap with any features

#### Rule Chart Styles
Percentage label darkness and bar colors reflect the magnitude of the rule's contribution. Magnitude is always considered on a 0-100% scale, rather than scaling down to the chart's view limits. These styles cannot be changed using a plot stylesheet.



## class_charts
Features can have multiple classes associated with them, so it is useful to see the proportions of counts by class. The class_charts plot type shows the percentage of _mapped_ reads that were assigned to features by class. Each feature's associated classes are determined by the `Class=` attribute in your GFF files.

<p float="left" align="center">
    <img src="../images/plots/class_charts.jpeg" width="60%" alt="rule_chart with 10 rules"/>
</p>

#### Class N
Class **N** represents the percentage of mapped reads that were unassigned. Sources of unassigned reads include:
- A lack of features passing selection at alignment loci
- Alignments which do not overlap with any features

#### Count Normalization
Features with multiple associated classes will have their counts split evenly across these classes before being grouped and summed.

#### Class Chart Styles
Proportions in rule_charts and class_charts are plotted using the same function. Styles are the same between the two. See [rule chart styles](#rule-chart-styles) for more info.


## replicate_scatter
Todo

## sample_avg_scatter_by_dge
Todo

## sample_avg_scatter_by_dge_class
Todo