#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
  - class: ScatterFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:

  # multi input
  threads: int?

  # fastp inputs
  in_fq: File[]
  out_fq: string[]
  fp_phred64: boolean?
  compression: int?
  dont_overwrite: boolean?
  disable_adapter_trimming: boolean?
  adapter_sequence: string?
  trim_poly_x: boolean?
  poly_x_min_len: int?
  disable_quality_filtering: boolean?
  qualified_quality_phred: int?
  unqualified_percent_limit: int?
  n_base_limit: int?
  disable_length_filtering: boolean?
  length_required: int?
  length_limit: int?
  overrepresentation_analysis: boolean?
  overrepresentation_sampling: int?
  json: string[]
  html: string[]
  report_title: string[]

  # collapser inputs
  uniq_seq_prefix: string[]
  threshold: int?
  compress: boolean?

  # bowtie inputs
  bt_index_files: File[]
  ebwt: string
  outfile: string[]
  fastq: boolean?
  fasta: boolean?
  trim5: int?
  trim3: int?
  bt_phred64: boolean?
  solexa: boolean?
  solexa13: boolean?
  end_to_end: int?
  nofw: boolean?
  norc: boolean?
  k_aln: int?
  all_aln: boolean?
  no_unal: boolean?
  un: string[]
  sam: boolean?
  seed: int?
  shared_memory: boolean?

  #counter inputs
  output_prefix: string
  samples_csv: File
  features_csv: File
  intermed_file: boolean?
  gff_files: File[]?
  aligned_seqs: File[]?
  is_pipeline: boolean?

steps:

  counts-prep:
    run: per-library.cwl
    scatter: [in_fq, out_fq, json, html, report_title, uniq_seq_prefix, outfile, un]
    scatterMethod: dotproduct
    in:
      # fastp
      in_fq: in_fq
      out_fq: out_fq
      fp_phred64: fp_phred64
      compression: compression
      dont_overwrite: dont_overwrite
      disable_adapter_trimming: disable_adapter_trimming
      adapter_sequence: adapter_sequence
      trim_poly_x: trim_poly_x
      poly_x_min_len: poly_x_min_len
      disable_quality_filtering: disable_quality_filtering
      qualified_quality_phred: qualified_quality_phred
      unqualified_percent_limit: unqualified_percent_limit
      n_base_limit: n_base_limit
      disable_length_filtering: disable_length_filtering
      length_required: length_required
      length_limit: length_limit
      overrepresentation_analysis: overrepresentation_analysis
      overrepresentation_sampling: overrepresentation_sampling
      json: json
      html: html
      report_title: report_title

      # Collapser
      uniq_seq_prefix: uniq_seq_prefix
      threshold: threshold
      compress: compress

      # Bowtie
      bt_index_files: bt_index_files
      ebwt: ebwt
      outfile: outfile
      fastq: fastq
      fasta: fasta
      trim5: trim5
      trim3: trim3
      bt_phred64: bt_phred64
      solexa: solexa
      solexa13: solexa13
      end_to_end: end_to_end
      nofw: nofw
      k_aln: k_aln
      all_aln: all_aln
      no_unal: no_unal
      un: un
      sam: sam
      threads: threads
      shared_memory: shared_memory
      seed: seed
    out: [fastq_clean, html_report_file, json_report_file, uniq_seqs, aln_seqs, uniq_seqs_low]

  counts:
    run: ../tools/aquatx-count.cwl
    in:
      aligned_seqs: counts-prep/aln_seqs
      gff_files: gff_files
      samples_csv: samples_csv
      config_csv: features_csv
      out_prefix: output_prefix
      intermed_file: intermed_file
      fastp_logs: counts-prep/json_report_file
      collapsed_fa: counts-prep/uniq_seqs
      is_pipeline: {default: true}
    out: [feature_counts, other_counts, alignment_stats, summary_stats, intermed_out_files]

  deseq2:
    run: ../tools/aquatx-deseq.cwl
    in:
      input_file: counts/feature_counts
      outfile_prefix: output_prefix
    out: [norm_counts, comparisons]

outputs:
  # Per-library outputs
  fastq_clean:
    type: File[]
    outputSource: counts-prep/fastq_clean

  html_report_file:
    type: File[]
    outputSource: counts-prep/html_report_file

  json_report_file:
    type: File[]
    outputSource: counts-prep/json_report_file

  uniq_seqs:
    type: File[]
    outputSource: counts-prep/uniq_seqs

  aln_seqs:
    type: File[]
    outputSource: counts-prep/aln_seqs

  other_count_files:
    type: File[]
    outputSource: counts/other_counts

  # Pipeline summary outputs from aquatx-count
  feat_count_file:
    type: File
    outputSource: counts/feature_counts

  alignment_stats:
    type: File
    outputSource: counts/alignment_stats

  summary_stats:
    type: File
    outputSource: counts/summary_stats

  deseq_normed:
    type: File
    outputSource: deseq2/norm_counts

  deseq_tables:
    type: File[]
    outputSource: deseq2/comparisons

  # Optional outputs
  aln_tables:
    type: File[]?
    outputSource: counts/intermed_out_files