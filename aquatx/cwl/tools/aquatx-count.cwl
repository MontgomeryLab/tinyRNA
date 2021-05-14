#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing: $(inputs.gff_files)

baseCommand: aquatx-count

inputs:
  samples_csv:
    type: File
    inputBinding:
      prefix: -i
      position: 0

  config_csv:
    type: File
    inputBinding:
      prefix: -c
      position: 1

  out_prefix:
    type: string
    inputBinding:
      position: 2
      prefix: -o

  intermed_file:
    type: boolean?
    inputBinding:
      position: 3
      prefix: -t

  is_pipeline:
    type: boolean?
    inputBinding:
      position: 4
      prefix: -p

  # Specifies the GFF files defined in features.csv
  gff_files:  # This optional input is for pipeline execution.
    type: File[]?

  # These optional inputs are for producing pipeline summary statistics at the conclusion of counting
  fastp_logs:
    type: File[]?

  collapsed_fa:
    type: File[]?

outputs:
  feature_counts:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_feature_counts.csv

  other_counts:
    type: File[]
    outputBinding:
      glob: $(inputs.out_prefix)_nt_len_dist.csv

  stats_file:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_alignment_stats.csv

  intermed_out_file:
    type: File[]?
    outputBinding:
      glob: "*_aln_table.txt"
