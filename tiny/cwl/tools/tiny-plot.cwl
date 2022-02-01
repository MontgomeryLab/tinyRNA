#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: tiny-plot
stdout: console_output.log

inputs:

  raw_counts:
    type: File
    inputBinding:
      prefix: -rc
    doc: "Raw, non-normalized feature counts from Counter"

  norm_counts:
    type: File?
    inputBinding:
      prefix: -nc
    doc: "Normalized feature counts from DESeq"

  dge_tables:
    type: File[]?
    inputBinding:
      prefix: -dge
    doc: "Sample comparison tables from DESeq"

  summ_stats:
    type: File
    inputBinding:
      prefix: -ss
    doc: "The summary stats csv from Counter"

  len_dist:
    type: File[]
    inputBinding:
      prefix: -len
    doc: "5' end nucleotide vs. length matrices from Counter"

  dge_pval:
    type: float?
    inputBinding:
      prefix: -pv
    doc: "The p-value to use for DGE scatter plots (default: 0.05)"

  style_sheet:
    type: File?
    inputBinding:
      prefix: -s
    doc: "A .mplstyle sheet to use instead of tinyrna default styles"

  vector_scatter:
    type: boolean?
    inputBinding:
      prefix: -v
    doc: "If provided, scatter plots will have vectorized points (slower)"

  out_prefix:
    type: string?
    inputBinding:
      prefix: -o
    doc: "The prefix to use when naming output files (optional)"

  plot_requests:
    type: string[]
    inputBinding:
      prefix: -p
    doc: "A list of desired plot types to produce"

outputs:
  plots:
    type: File[]?
    outputBinding:
      glob: "*.pdf"

  console_output:
    type: stdout