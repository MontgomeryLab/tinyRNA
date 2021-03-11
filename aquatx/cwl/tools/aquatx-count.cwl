#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: aquatx-count

inputs:
  input_files:
    type: File[]
    inputBinding:
      prefix: -i
      itemSeparator: ','
      position: 0

  config_file:
    type: File
    inputBinding:
      prefix: -c
      position: 1

  out_prefix:
    type: string[]
    inputBinding:
      position: 2
      prefix: -o

  intermed_file:
    type: boolean?
    inputBinding:
      position: 3
      prefix: -t

outputs:
  feature_counts:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_out_feature_counts.txt

  other_counts:
    type: File[]
    outputBinding:
      glob:
        - $(inputs.out_prefix)_out_nt_len_dist.csv
        - $(inputs.out_prefix)_out_class_counts.csv

  stats_file:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_stats.txt

  intermed_out_file:
    type: File[]?
    outputBinding:
      glob: $(*_out_aln_table.txt)
