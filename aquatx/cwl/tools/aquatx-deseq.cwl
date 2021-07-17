#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: aquatx-deseq.r

inputs:
  input_file:
    type: File
    inputBinding:
      prefix: --input-file
    doc: The merged count table output of aquatx-count

  outfile_prefix:
    type: string
    inputBinding:
      prefix: --outfile-prefix
    doc: The prefix for naming output files

  pca_plots:
    type: boolean?
    inputBinding:
      prefix: --pca
    doc: Produce PCA plots for each library comparison

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
      glob: $(inputs.outfile_prefix)_*_*_pca_plot.pdf
