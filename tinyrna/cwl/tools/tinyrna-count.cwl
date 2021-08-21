#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing: [
      $(inputs.gff_files),
      $(inputs.aligned_seqs),
      $(inputs.fastp_logs),
      $(inputs.collapsed_fa)
    ]

baseCommand: tinyrna-count

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
    default: false
    inputBinding:
      position: 4
      prefix: -p

  diagnostics:
    type: boolean?
    default: false
    inputBinding:
      position: 5
      prefix: -d

  # The following optional inputs are for staging InitialWorkingDir files for pipeline execution

  # Specifies the GFF files defined in features.csv
  gff_files:
    type: File[]?

  # Specifies the fastp json logs for producing pipeline Summary Statistics
  fastp_logs:
    type: File[]?

  # Specifies the collapsed Fasta files produced by tinyrna-collapse for producing pipeline Summary Statistics
  collapsed_fa:
    type: File[]?

  # Specifies the .sam files produced by the bowtie step of the pipeline
  aligned_seqs:
    type: File[]?

outputs:
  feature_counts:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_feature_counts.csv

  other_counts:
    type: File[]
    outputBinding:
      glob: $(inputs.out_prefix)*_nt_len_dist.csv

  alignment_stats:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_alignment_stats.csv

  summary_stats:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_summary_stats.csv

  intermed_out_files:
    type: File[]?
    outputBinding:
      glob: "*_aln_table.txt"

  alignment_diags:
    type: File?
    outputBinding:
      glob: $(inputs.out_prefix)_alignment_diags.csv

  selection_diags:
    type: File?
    outputBinding:
      glob: $(inputs.out_prefix)_selection_diags.txt