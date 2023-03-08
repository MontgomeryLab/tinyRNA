# Options for tinyRNA Tools
This page provides an explanation of the parameters offered by each of our Python components. We represent each parameter with its Run Config key and its corresponding commandline argument. When running the tinyRNA pipeline (`tiny run/recount/replot`), you'll use the Run Config key to specify your preferences. The commandline argument is used when you run a tool as an individual, standalone step (`tiny-collapse`, `tiny-count`, or `tiny-plot`)

### Table of Contents
- [tiny-collapse](#tiny-collapse)
- [tiny-count](#tiny-count)
- [tiny-deseq.r](#tiny-deseqr)
- [tiny-plot](#tiny-plot)

## tiny-collapse

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

tiny-collapse outputs are often very large. You can save space by switching this option "on" so that outputs are gzipped before being written to disk.

### Full tiny-collapse Help String
```
tiny-collapse -i FASTQFILE -o OUTPREFIX [-h] [-t THRESHOLD] [-c]
              [--5p-trim LENGTH] [--3p-trim LENGTH]

Collapse sequences from a fastq file to a fasta file. Headers in the output
fasta file will contain the number of times each sequence occurred in the
input fastq file, and an ID which indicates the relative order in which each
sequence was first encountered. Gzipped files are automatically supported for
fastq inputs, and compressed fasta outputs are available by request.

Required arguments:
  -i FASTQFILE, --input-file FASTQFILE
                        The input fastq(.gz) file to collapse
  -o OUTPREFIX, --out-prefix OUTPREFIX
                        The prefix for output files {prefix}_collapsed.fa and,
                        if counts fall below threshold,
                        {prefix}_collapsed_lowcounts.fa

Optional arguments:
  -h, --help            show this help message and exit
  -t THRESHOLD, --threshold THRESHOLD
                        Sequences <= THRESHOLD will be omitted from
                        {prefix}_collapsed.fa and will instead be placed in
                        {prefix}_collapsed_lowcounts.fa
  -c, --compress        Use gzip compression when writing fasta outputs
  --5p-trim LENGTH      Trim LENGTH bases from the 5' end of each sequence
  --3p-trim LENGTH      Trim LENGTH bases from the 3' end of each sequence
```

## tiny-count

### Get Templates
| Run Config Key | Commandline Argument |
|----------------|----------------------|
|                | `--get-templates`    |

Copies the template configuration files required by tiny-count into the current directory. This argument can't be combined with `--paths-file`. All other arguments are ignored when provided, and once the templates have been copied tiny-count exits.

### Normalize by Hits
| Run Config Key             | Commandline Argument      |
|----------------------------|---------------------------|
| counter-normalize-by-hits: | `--normalize-by-hits T/F` |

By default, tiny-count will divide the number of counts associated with each sequence, twice, before they are assigned to a feature. Each unique sequence's count is determined by tiny-collapse (or a compatible collapsing utility) and is preserved through the alignment process. The original count is divided first by the number of loci that the sequence aligns to, and second by the number of features passing selection at each locus. Switching this option "off" disables the latter normalization step.

### Decollapse
| Run Config Key      | Commandline Argument   |
|---------------------|------------------------|
| counter_decollapse: | `--decollapse`         |

The SAM files produced by the tinyRNA pipeline are collapsed by default; alignments sharing a SEQ field are strictly multi-alignments and do not reflect original sequence counts. If this option is switched "on", tiny-count will produce a decollapsed copy of each input SAM file. Each alignment in the decollapsed SAM will be duplicated by the sequence's original count. This is useful for browsing in IGV. If non-collapsed inputs are provided to tiny-count in standalone mode, this option will be ignored.

### StepVector
| Run Config Key     | Commandline Argument |
|--------------------|----------------------|
| counter_stepvector | `--stepvector`       |

A custom Cython implementation of HTSeq's StepVector is used for finding features that overlap each alignment interval. While the core C++ component of the StepVector is the same, we have found that our Cython implementation can result in runtimes up to 50% faster than HTSeq's implementation. This parameter allows you to use HTSeq's StepVector if you wish.

### Is Pipeline
| Run Config Key | Commandline Argument |
|----------------|----------------------|
|                | `--is-pipeline`      |

This commandline argument tells tiny-count that it is running as a workflow step rather than a standalone/manual run. Under these conditions tiny-count will look for all input files in the current working directory regardless of the paths defined in the Samples Sheet and Features Sheet.

### Report Diags
| Run Config Key | Commandline Argument |
|----------------|----------------------|
| counter_diags: | `--report-diags`     |

Diagnostic information will include intermediate alignment files for each library and an additional stats table with information about counts that were not assigned to a feature. See [the description of these outputs](../README.md#Diagnostics) for details.

### Full tiny-count Help String
```
tiny-count (-pf FILE | --get-templates) [-o PREFIX] [-nh T/F] [-dc]
           [-sv {Cython,HTSeq}] [-p] [-d]

tiny-count is a precision counting tool for hierarchical classification and
quantification of small RNA-seq reads

Required arguments:
  You must either provide a Paths File or request templates for detailing
  your configuration.

  -pf FILE, --paths-file FILE
                        your Paths File (default: None)
  --get-templates       Copies the template configuration files required by
                        tiny-count into the current directory. (default:
                        False)

Optional arguments:
  These options can be used in conjunction with the Paths File (-pf)
  argument mentioned above.

  -o PREFIX, --out-prefix PREFIX
                        The output prefix to use for file names. All
                        occurrences of the substring {timestamp} will be
                        replaced with the current date and time. (default:
                        tiny-count_{timestamp})
  -nh T/F, --normalize-by-hits T/F
                        If T/true, normalize counts by (selected) overlapping
                        feature counts. (default: T)
  -dc, --decollapse     Create a decollapsed copy of all SAM files listed in
                        your Samples Sheet. This option is ignored for non-
                        collapsed inputs. (default: False)
  -sv {Cython,HTSeq}, --stepvector {Cython,HTSeq}
                        Select which StepVector implementation is used to find
                        features overlapping an interval. (default: Cython)
  -p, --is-pipeline     Indicates that tiny-count was invoked as part of a
                        pipeline run and that input files should be sourced as
                        such. (default: False)
  -d, --report-diags    Produce diagnostic information about
                        uncounted/eliminated selection elements. (default:
                        False)
```
## tiny-deseq.r

### PCA Plot
| Run Config Key | Commandline Argument |
|----------------|----------------------|
| dge_pca_plot   | `--pca`              |

DESeq2 can produce a PCA plot for your samples if your experiment contains biological replicates. If this option is switched "on", the PCA plot will be produced and placed in the plots subdirectory for each run.

### Drop Zero
| Run Config Key | Commandline Argument |
|----------------|----------------------|
| dge_drop_zero  | `--drop-zero`        |

Features with zero counts across all libraries will be dropped before DGE analysis if this option is switched "on".


### Full tiny-deseq.r Help String
```
tiny-deseq.r --input-file COUNTFILE --outfile-prefix PREFIX [--control CONDITION] [--pca] [--drop-zero]

Required arguments:

    --input-file <count_file>
          A text file containing a table of features x samples of the run to
          process by DESeq2. The [...]feature_counts.csv output of tinyrna-count is expected here.
              
    --outfile-prefix <outfile>
          Name of the output files to write. These will be created:
              1. Normalized count table of all samples
              2. Differential gene expression table per comparison
              3. A PCA plot per comparison, if --pca is also provided.

Optional arguments:

    --control <control_condition>
          If the control condition is specified, comparisons will
          only be made between the control and experimental conditions.

    --pca
          This will produce principle component analysis plots
          using the DESeq2 library. Output files are PDF format.

    --drop-zero
          Prior to performing analysis, this will drop all
          rows/features which have a zero count in all samples."
```

## tiny-plot

### Plot Requests
| Run Config Key | Commandline Argument         |
|----------------|------------------------------|
| plot_requests: | `--plots PLOT PLOT PLOT ...` |

tiny-plot will only produce the list of plots requested.

### P value
| Run Config Key | Commandline Argument |
|----------------|----------------------|
| plot_pval:     | `--p-value VALUE`    |

Feature expression levels are considered significant if their P value is less than this value, with a default of 0.05. Non-differentially expressed features are plotted as gray points, and in `sample_avg_scatter_by_dge_class`, these points are not colored by feature class.

### Style Sheet
| Run Config Key | Paths File Key    | Commandline Argument     |
|----------------|-------------------|--------------------------|
|                | plot_style_sheet: | `--style-sheet MPLSTYLE` |

The plot style sheet can be used to override the default Matplotlib styles used by tiny-plot. Unlike the other parameters, this option is found in the Paths File. The expected value for this parameter is the path to your modified style sheet. See the [Plot Style Sheet documentation](Configuration.md#plot-style-sheet-details) for more information.

### Vector Scatter
| Run Config Key      | Commandline Argument |
|---------------------|----------------------|
| plot_vector_points: | `--vector-scatter`   |

The scatter plots produced by tiny-plot have rasterized points by default. This allows for faster plot generation, smaller file sizes, and files that are more easily handled by PDF readers. Plots are produced in 300 dpi by default, so in most cases this rasterization is seldom noticeable under normal zoom levels. Switching this option "on" will cause points to be vectorized allowing for zooming without pixelation.
>**Note**: only scatter points are rasterized with this option switched "off"; all other elements are vectorized in every plot type.

### Bounds for len_dist Charts
| Run Config Key     | Commandline Argument   |
|--------------------|------------------------|
| plot_len_dist_min: | `--len-dist-min VALUE` | 
| plot_len_dist_max: | `--len-dist-max VALUE` |

The min and/or max bounds for plotted lengths can be set with this option. See [tiny-plot's documentation](tiny-plot.md#length-bounds) for more information about how these values are determined if they aren't set.

### Bounds for scatter_by_dge Plots
| Run Config Key        | Commandline Argument |
|-----------------------|----------------------|
| plot_dge_scatter_min: | `--dge-min VALUE`    | 
| plot_dge_scatter_max: | `--dge-max VALUE`    |

The min and/or max bounds for DGE scatter plots can be set with this option. The value you provide should be a log2 count value and can be whole or fractional, e.g. `--dge-min 1.9` would produce a plot whose first tick mark is labeled 2 and would include points for feature counts as low as 3.74. Unspecified bounds are automatically calculated to fit the data, and will include the margin specified by the `axes.[x/y]margin` key in the [Plot Style Sheet](Configuration.md#plot-style-sheet-details).

### Labels for Class-related Plots
| Run Config Key         | Commandline Argument |
|------------------------|----------------------|
| plot_unknown_class:    | `--unknown-class`    | 
| plot_unassigned_class: | `--unassigned-class` |

The labels that should be used for special groups in `class_charts` and `sample_avg_scatter_by_dge_class` plots. The "unknown" class group represents counts which were assigned by a Features Sheet rule which lacked a "Classify as..." label. The "unassigned" class group represents counts which weren't assigned to a feature.

### Filtering Classes in DGE Class Scatter Plots
| Run Config Key             | Commandline Argument |
|----------------------------|----------------------|
| plot_class_scatter_filter: | `--classes-include`  |
|                            | `--classes-exclude`  |

If an inclusive filter is used, then only the classes in the list, if present, are shown. If an exclusive filter is used, then the listed classes are omitted from the plot. This behavior extends to features whose P value is above threshold. In the Run Config, the filter type can be set with the `style:` sub-key, and the desired list of classes for the filter can be provided between the brackets of the `classes:` sub-key

### Full tiny-plot Help String
```
tiny-plot [-rc RAW_COUNTS] [-nc NORM_COUNTS] [-uc RULE_COUNTS]
          [-ss STAT] [-dge COMPARISON [COMPARISON ...]]
          [-len 5P_LEN [5P_LEN ...]] [-o PREFIX] [-pv VALUE]
          [-s MPLSTYLE] [-v] [-ldi VALUE] [-lda VALUE] [-dgi VALUE]
          [-dga VALUE] [-una LABEL] [-unk LABEL] [-ic CLASS [CLASS ...]
          | -ec CLASS [CLASS ...]] -p PLOT [PLOT ...]

This script produces basic static plots for publication as part of the tinyRNA
workflow.

Input file requirements vary by plot type and you are free to supply only the
files necessary for your plot selections. If you are sourcing all of your
input files from the same run directory, you may find it easier to instead run
`tiny replot` within that run directory.

Required arguments:
  -p PLOT [PLOT ...], --plots PLOT [PLOT ...]
                        List of plots to create. Options:
                        • len_dist: A stacked barchart showing size & 5'
                          nucleotide distribution.
                        • rule_charts: A barchart showing percentages of
                          counts by matched rule.
                        • class_charts: A barchart showing percentages of
                          counts per classification.
                        • replicate_scatter: A scatter plot comparing
                          replicates for all count files given.
                        • sample_avg_scatter_by_dge: A scatter plot comparing
                          all sample groups, with differentially expressed
                          small RNAs highlighted based on P value cutoff.
                        • sample_avg_scatter_by_dge_class: A scatter plot
                          comparing all sample groups, with classes
                          highlighted for differentially expressed small RNAs
                          based on P value cutoff.

Input files produced by tiny-count:
  -rc RAW_COUNTS, --raw-counts RAW_COUNTS
                        The ...feature_counts.csv file
  -uc RULE_COUNTS, --rule-counts RULE_COUNTS
                        The ...counts-by-rule.csv file
  -ss STAT, --summary-stats STAT
                        The ...summary_stats.csv file
  -len 5P_LEN [5P_LEN ...], --len-dist 5P_LEN [5P_LEN ...]
                        The ...nt_len_dist.csv files

Input files produced by tiny-deseq.r:
  -nc NORM_COUNTS, --norm-counts NORM_COUNTS
                        The ...norm_counts.csv file
  -dge COMPARISON [COMPARISON ...], --dge-tables COMPARISON [COMPARISON ...]
                        The ...cond1...cond2...deseq.csv files

Optional arguments:
  -o PREFIX, --out-prefix PREFIX
                        Prefix to use for output filenames.
  -pv VALUE, --p-value VALUE
                        P value to use in DGE scatter plots.
  -s MPLSTYLE, --style-sheet MPLSTYLE
                        Optional matplotlib style sheet to use for plots.
  -v, --vector-scatter  Produce scatter plots with vectorized points (slower).
                        Note: only the points on scatter plots will be raster
                        if this option is not provided.
  -ldi VALUE, --len-dist-min VALUE
                        len_dist plots will start at this value
  -lda VALUE, --len-dist-max VALUE
                        len_dist plots will end at this value
  -dgi VALUE, --dge-min VALUE
                        scatter_by_dge plots will start at this log2 value
  -dga VALUE, --dge-max VALUE
                        scatter_by_dge plots will end at this log2 value
  -una LABEL, --unassigned-class LABEL
                        Use this label in class-related plots for unassigned
                        counts
  -unk LABEL, --unknown-class LABEL
                        Use this label in class-related plots for counts which
                        were assigned by rules lacking a "Classify as..."
                        value
  -ic CLASS [CLASS ...], --classes-include CLASS [CLASS ...]
                        Only include these classes, if present, in class
                        scatter plots (applies regardless of P value)
  -ec CLASS [CLASS ...], --classes-exclude CLASS [CLASS ...]
                        Omit these classes, if present, from class scatter
                        plots (applies regardless of P value)
```