#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: tiny-plot
stdout: console_output.log

inputs:

  raw_counts:
    type: File
    inputBinding:
      prefix: --raw-counts
    doc: "Raw, non-normalized feature counts from Counter"

  rule_counts:
    type: File
    inputBinding:
      prefix: --rule-counts
    doc: "Raw, non-normalized counts by matched rule from Counter"

  norm_counts:
    type: File?
    inputBinding:
      prefix: --norm-counts
    doc: "Normalized feature counts from DESeq"

  dge_tables:
    type: File[]?
    inputBinding:
      prefix: --dge-tables
    doc: "Sample comparison tables from DESeq"

  summ_stats:
    type: File
    inputBinding:
      prefix: --summary-stats
    doc: "The summary stats csv from Counter"

  len_dist_tables:
    type: File[]
    inputBinding:
      prefix: --len-dist
    doc: "5' end nucleotide vs. length matrices from Counter"

  dge_pval:
    type: float?
    inputBinding:
      prefix: --p-value
    doc: "The p-value to use for DGE scatter plots (default: 0.05)"

  style_sheet:
    type: File?
    inputBinding:
      prefix: --style-sheet
    doc: "A .mplstyle sheet to use instead of tinyrna default styles"

  vector_scatter:
    type: boolean?
    inputBinding:
      prefix: --vector-scatter
    doc: "If provided, scatter plots will have vectorized points (slower)"

  len_dist_min:
    type: int?
    inputBinding:
      prefix: --len-dist-min
    doc: "The first length to plot in the range for len_dist plots"

  len_dist_max:
    type: int?
    inputBinding:
      prefix: --len-dist-max
    doc: "The last length to plot in the range for len_dist plots"

  dge_min:
    type: double?
    inputBinding:
      prefix: --dge-min
    doc: "The log2 lower view limit in DGE scatter plots"

  dge_max:
    type: double?
    inputBinding:
      prefix: --dge-max
    doc: "The log2 upper view limit in DGE scatter plots"

  unknown_class_label:
    type: string?
    inputBinding:
      prefix: --unknown-class
    doc: \
      'Use this label in class-related plots for counts which were '
      'assigned by rules lacking a "Classify as..." value'

  unassigned_class_label:
    type: string?
    inputBinding:
      prefix: --unassigned-class
    doc: 'Use this label in class-related plots for unassigned counts'

  classes_include:
    type: string[]?
    inputBinding:
      prefix: --classes-include
    doc: \
      'Only include these classes, if present, in class scatter '
      'plots (applies regardless of P value)'

  classes_exclude:
    type: string[]?
    inputBinding:
      prefix: --classes-exclude
    doc: \
      'Omit these classes, if present, from class scatter plots '
      '(applies regardless of P value)'

  out_prefix:
    type: string?
    inputBinding:
      prefix: --out-prefix
    doc: "The prefix to use when naming output files (optional)"

  plot_requests:
    type: string[]
    inputBinding:
      prefix: --plots
    doc: "A list of desired plot types to produce"

outputs:

  len_dist:
    type: Directory?
    outputBinding:
      glob: len_dist

  rule_chart:
    type: Directory?
    outputBinding:
      glob: rule_chart

  class_chart:
    type: Directory?
    outputBinding:
      glob: class_chart

  replicate_scatter:
    type: Directory?
    outputBinding:
      glob: replicate_scatter

  sample_avg_scatter_by_dge:
    type: Directory?
    outputBinding:
      glob: scatter_by_dge

  sample_avg_scatter_by_dge_class:
    type: Directory?
    outputBinding:
      glob: scatter_by_dge_class

  console_output:
    type: stdout