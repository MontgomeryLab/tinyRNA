#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InitialWorkDirRequirement
    listing: [
      $(inputs.samples_csv),
      $(inputs.features_csv),
      $(inputs.gff_files),
      $(inputs.aligned_seqs),
      $(inputs.fastp_logs),
      $(inputs.collapsed_fa)
    ]

baseCommand: tiny-count
stdout: console_output.log

inputs:
  paths_file:
    type: File
    inputBinding:
      prefix: --paths-file

  out_prefix:
    type: string
    inputBinding:
      prefix: --out-prefix

  # Optional inputs

  normalize_by_feature_hits:
    type: string?
    inputBinding:
      prefix: --normalize-by-feature-hits

  normalize_by_genomic_hits:
    type: string?
    inputBinding:
      prefix: --normalize-by-genomic-hits

  decollapse:
    type: boolean?
    inputBinding:
      prefix: --decollapse

  stepvector:
    type: string?
    inputBinding:
      prefix: --stepvector

  in_pipeline:
    type: boolean?
    inputBinding:
      prefix: --in-pipeline

  diagnostics:
    type: boolean?
    inputBinding:
      prefix: --report-diags

  # The following optional inputs are for staging InitialWorkingDir files for pipeline execution

  samples_csv:
    type: File

  features_csv:
    type: File

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

  norm_counts:
    type: File[]?
    outputBinding:
      glob: $(inputs.out_prefix)_norm_counts.csv

  rule_counts:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_counts_by_rule.csv

  mapped_nt_len_dist:
    type: File[]
    outputBinding:
      glob: $(inputs.out_prefix)*_mapped_nt_len_dist.csv

  assigned_nt_len_dist:
    type: File[]
    outputBinding:
      glob: $(inputs.out_prefix)*_assigned_nt_len_dist.csv

  alignment_stats:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_alignment_stats.csv

  summary_stats:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_summary_stats.csv

  decollapsed_sams:
    type: File[]?
    outputBinding:
      glob: "*_decollapsed.sam"

  alignment_tables:
    type: File[]?
    outputBinding:
      glob: "*_alignment_table.csv"

  assignment_diags:
    type: File?
    outputBinding:
      glob: $(inputs.out_prefix)_assignment_diags.csv

  selection_diags:
    type: File?
    outputBinding:
      glob: $(inputs.out_prefix)_selection_diags.txt

  stats_check:
    type: File?
    outputBinding:
      glob: "*_stats_check.csv"

  console_output:
    type: stdout