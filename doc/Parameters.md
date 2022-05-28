# Options for tinyRNA Tools
This page provides an explanation of the parameters offered by each of our Python components. We represent each parameter with its Run Config key and its corresponding commandline argument. When running the tinyRNA pipeline (`tiny run/recount/replot`), you'll use the Run Config key to specify your preferences. The commandline argument is used when you run a tool as an individual, standalone step (`tiny-collapse`, `tiny-count`, or `tiny-plot`)

### Table of Contents
- [Collapser](#collapser)
- [Counter](#counter)
- [Plotter](#plotter)

## Collapser

### Threshold
| Run Config Key | Commandline Argument     |
|----------------|--------------------------|
| threshold:     | `--threshold THRESHOLD`  |

You can specify a minimum count threshold to determine which sequences are reported after unique sequence identification and counting is complete. Sequences with a count less than or equal to the threshold value are placed in a `low_count.fasta` file, and are not used in downstream analyses.

### Trimming
| Run Config Key | Commandline Argument |
|----------------|----------------------|
| 5p_trim:       | `--5p-trim LENGTH`   |
| 3p_trim:       | `--3p-trim LENGTH`   |

Bases can be trimmed from the 5' and/or 3' end of each sequence before it is evaluated for uniqueness and counted. This is useful for trimming randomized bases and UMIs (note: formal UMI deduplication is not offered at this time).

### Compression
| Run Config Key | Commandline Argument |
|----------------|----------------------|
| compress:      | `--compress`         |

Collapser outputs are often very large. You can save space by switching this option "on" so that outputs are gzipped before being written to disk.

## Counter

### All-Features
| Run Config Key         | Commandline Argument   |
|------------------------|------------------------|
| counter_all_features:  | `--all-features`       |

By default, Counter will only evaluate alignments to features which match a `Select for...` & `with value...` of at least one rule in your Features Sheet. It is this matching feature set, and only this set, which is included in `feature_counts.csv` and therefore available for analysis by DESeq2 and Plotter. Switching this option "on" will include all features in every input GFF file, regardless of attribute matches, for Counter and downstream steps.

### Normalize-by-Hits
 | Run Config Key             | Commandline Argument      |
|----------------------------|---------------------------|
| counter-normalize-by-hits: | `--normalize-by-hits T/F` |

By default, Counter will divide the number of counts associated with each sequence, twice, before they are assigned to a feature. Each unique sequence's count is determined by Collapser and is preserved through the alignment process. The original count is divided first by the number of loci that the sequence aligns to, and second by the number of features passing selection at each locus. Switching this option "off" disables the latter normalization step.

### Decollapse
 | Run Config Key      | Commandline Argument   |
|---------------------|------------------------|
| counter_decollapse: | `--decollapse`         |

The SAM files produced by the tinyRNA pipeline are collapsed by default; alignments sharing a SEQ field are strictly multi-alignments and do not reflect original sequence counts. If this option is switched "on", Counter will produce a decollapsed copy of each input SAM file. Each alignment in the decollapsed SAM will be duplicated by the sequence's original count. This is useful for browsing in IGV.

### Filters
 | Run Config Key                 | Commandline Argument                       |
|--------------------------------|--------------------------------------------|
| counter_source_filter: [ ]     | `--source-filter SOURCE SOURCE SOURCE ...` | 
| counter_type_filter: [ ]       | `--type-filter TYPE TYPE TYPE ...`         |

You can optionally filter features in your GFF files by specifying sources and/or types that are desired. Source and type refer to GFF columns 2 and 3 respectively. If source _and_ type filters are specified, each feature must match one of the sources _and_ one of the types in order to be included in the counting process. For both filters, an empty list is the same as "allow all."

### Is Pipeline
 | Run Config Key | Commandline Argument |
|----------------|----------------------|
|                | `--is-pipeline`      |

This commandline argument tells Counter that it is running as a workflow step rather than a standalone/manual run. Under these conditions Counter will look for all input files in the current working directory regardless of the paths defined in the Samples Sheet and Features Sheet.

### Report Diags
 | Run Config Key | Commandline Argument |
|----------------|----------------------|
| counter_diags: | `--report-diags`     |

Diagnostic information will include intermediate alignment files for each library and an additional stats table with information about counts that were not assigned to a feature. See [the description of these outputs](../README.md#Diagnostics) for details.

## Plotter

### Plot Requests
 | Run Config Key | Commandline Argument         |
|----------------|------------------------------|
| plot_requests: | `--plots PLOT PLOT PLOT ...` |

Plotter will only produce the list of plots requested.

### P value
 | Run Config Key | Commandline Argument |
|----------------|----------------------|
| plot_pval:     | `--p-value VALUE`    |

Feature expression levels are considered significant if their P value is less than this value, with a default of 0.05. Non-differentially expressed features are plotted as gray points, and in `sample_avg_scatter_by_dge_class`, these points are not colored by feature class.

### Style Sheet
 | Run Config Key | Paths File Key    | Commandline Argument     |
|----------------|-------------------|--------------------------|
|                | plot_style_sheet: | `--style-sheet MPLSTYLE` |

The plot style sheet can be used to override the default Matplotlib styles used by Plotter. Unlike the other parameters, this option is found in the Paths File. See the [Plot Stylesheet documentation](Configuration.md#plot-stylesheet-details) for more information.

### Vector Scatter
 | Run Config Key      | Commandline Argument |
|---------------------|----------------------|
| plot_vector_points: | `--vector-scatter`   |

The scatter plots produced by Plotter have rasterized points by default. This allows for faster plot generation, smaller file sizes, and files that are more easily handled by PDF readers. Plots are produced in 300 dpi by default, so in most cases this rasterization is seldom noticeable under normal zoom levels. Switching this option "on" will cause points to be vectorized allowing for zooming without pixelation.
>**Note**: only scatter points are rasterized with this option switched "off"; all other elements are vectorized in every plot type.

### Bounds for len_dist Charts
 | Run Config Key     | Commandline Argument   |
|--------------------|------------------------|
| plot_len_dist_min: | `--len-dist-min VALUE` | 
| plot_len_dist_max: | `--len-dist-max VALUE` |

The min and/or max bounds for plotted lengths can be set with this option. See [Plotter's documentation](Plotter.md#length-bounds) for more information about how these values are determined if they aren't set.