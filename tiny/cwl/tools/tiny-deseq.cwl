#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: tiny-deseq.r
stdout: console_output.log
stderr: console_output.log

inputs:
  input_file:
    type: File
    inputBinding:
      prefix: --input-file
    doc: The merged count table output of tiny-count

  outfile_prefix:
    type: string
    inputBinding:
      prefix: --outfile-prefix
    doc: The prefix for naming output files

  control:
    type: string?
    inputBinding:
      prefix: --control
    doc: If specified, comparisons will only be made against the control condition

  plots:
    type: boolean?
    inputBinding:
      prefix: --pca
    doc: Produce PCA plots for each library comparison

  drop_zero:
    type: boolean?
    inputBinding:
      prefix: --drop-zero
    doc: Drop features which have a zero count across all samples before analysis

outputs:
  norm_counts:
    type: File
    outputBinding:
      glob: $(inputs.outfile_prefix)_norm_counts.csv

  comparisons:
    type: File[]
    outputBinding:
      glob: $(inputs.outfile_prefix)_*_*_deseq_table.csv

  pca_plots:
    type: File[]?
    outputBinding:
      glob: $(inputs.outfile_prefix)_pca_plot.pdf

  console_output:
    type: stdout