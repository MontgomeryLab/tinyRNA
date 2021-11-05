#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InitialWorkDirRequirement
    listing: [
      $(inputs.gff_files),
      $(inputs.aligned_seqs),
      $(inputs.fastp_logs),
      $(inputs.collapsed_fa)
    ]

baseCommand: tiny-count
stdout: console_output.log

inputs:
  samples_csv:
    type: File
    inputBinding:
      prefix: -i

  config_csv:
    type: File
    inputBinding:
      prefix: -c

  out_prefix:
    type: string
    inputBinding:
      prefix: -o

  # Optional inputs

  source_filter:
    type: string[]?
    inputBinding:
      prefix: -sf

  type_filter:
    type: string[]?
    inputBinding:
      prefix: -tf

  all_features:
    type: boolean?
    inputBinding:
      prefix: -a

  is_pipeline:
    type: boolean?
    inputBinding:
      prefix: -p

  diagnostics:
    type: boolean?
    inputBinding:
      prefix: -d

  # The following optional inputs are for staging InitialWorkingDir files for pipeline execution

  # Specifies the GFF files defined in features.csv
  gff_files:
    type: File[]?

  # Specifies the fastp json logs for producing pipeline Summary Statistics
  fastp_logs:
    type: File[]?

  # Specifies the collapsed Fasta files produced by tiny-collapse for producing pipeline Summary Statistics
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

  console_output:
    type: stdout