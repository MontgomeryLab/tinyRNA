#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: tiny-plot

inputs:
  raw_counts:
    type: File
    inputBinding:
      prefix: -rc
    doc: "Raw, non-normalized feature counts from Counter"

  norm_counts:
    type: File
    inputBinding:
      prefix: -nc
    doc: "Normalized feature counts from DESeq"

  dge_tables:
    type: File[]
    inputBinding:
      prefix: -dge
    doc: "Sample comparison tables from DESeq"

  len_dist:
    type: File[]
    inputBinding:
      prefix: -len
    doc: "5' end nucleotide vs. length matrices from Counter"

  style_sheet:
    type: File?
    inputBinding:
      prefix: -s
    doc: "A .mplstyle sheet to use instead of tinyrna default styles"

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

