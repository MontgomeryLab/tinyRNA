#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow

requirements:
  - class: ScatterFeatureRequirement
  - class: SubworkflowFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:

  # multi input
  threads: int?
  run_name: string

  # bowtie build
  run_bowtie_build: boolean
  reference_genome_files: File[]
  offrate: int?
  ntoa: boolean?
  noref: boolean?
  ftabchars: int?

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
  logfile: string[]
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

  # counter inputs
  samples_csv: File
  features_csv: File
  intermed_file: boolean?
  gff_files: File[]?
  aligned_seqs: File[]?
  is_pipeline: boolean?

  # plotter options
  plot_input_files: string[]
  plot_data_types: string[]
  plots: string[]

steps:

  bt_build_optional:
    run: ../tools/bowtie-build.cwl
    when: $(inputs.run_bowtie_build)
    in:
      run_bowtie_build: run_bowtie_build
      ref_in: reference_genome_files
      ebwt_base: ebwt
      offrate: offrate
      ntoa: ntoa
      noref: noref
      ftabchars: ftabchars
      threads: threads
    out: [index_files]

  counts-prep:
    run: per-library.cwl
    scatter: [in_fq, out_fq, json, html, report_title, uniq_seq_prefix, outfile, logfile, un]
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
      bt_index_files:
        source: [bt_build_optional/index_files, bt_index_files]
        pickValue: first_non_null
      ebwt: ebwt
      outfile: outfile
      logfile: logfile
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
    out: [fastq_clean, html_report_file, json_report_file, uniq_seqs, uniq_seqs_low, aln_seqs, unal_seqs, bowtie_log]

  counts:
    run: ../tools/aquatx-count.cwl
    in:
      aligned_seqs: counts-prep/aln_seqs
      gff_files: gff_files
      samples_csv: samples_csv
      config_csv: features_csv
      out_prefix: run_name
      intermed_file: intermed_file
      fastp_logs: counts-prep/json_report_file
      collapsed_fa: counts-prep/uniq_seqs
      is_pipeline: {default: true}
    out: [feature_counts, other_counts, alignment_stats, summary_stats, intermed_out_files]

  subdirs:
    run: organize-outputs.cwl
    in:
      bt_build_indexes: bt_build_optional/index_files
      run_bowtie_build: run_bowtie_build

      fastp_cleaned_fastq: counts-prep/fastq_clean
      fastp_html_report: counts-prep/html_report_file
      fastp_json_report: counts-prep/json_report_file

      collapser_uniq: counts-prep/uniq_seqs
      collapser_low:
        # If no output in scatter, this value will actually be an
        #  array of nulls rather than just null. Can't use default.
        source: counts-prep/uniq_seqs_low
        pickValue: all_non_null

      bowtie_sam: counts-prep/aln_seqs
      bowtie_log: counts-prep/bowtie_log
      bowtie_unal:
        source: counts-prep/unal_seqs
        default: []

      counter_features: counts/feature_counts
      counter_other: counts/other_counts
      counter_alignment_stats: counts/alignment_stats
      counter_summary_stats: counts/summary_stats
      counter_intermed:
        source: counts/intermed_out_files
        default: []
    out: [ bt_build_dir, fastp_dir, collapser_dir, bowtie_dir, counter_dir ]

  differential-expression:
    run: ../tools/aquatx-deseq.cwl
    in:
      input_file: counts/feature_counts
      outfile_prefix: run_name
    out: [ norm_counts, comparisons ]

  plotter:
    run: ../tools/aquatx-plot.cwl
    in:
      input_files: plot_input_files
      data_types: plot_data_types
      out_prefix: output_prefix
      reference: ref_annotations
      plots: plots
    out: [plots]

outputs:

  # Subdirectory outputs
  bt_build_out_dir:
    type: Directory?
    outputSource: subdirs/bt_build_dir

  fastp_out_dir:
    type: Directory
    outputSource: subdirs/fastp_dir

  collapser_out_dir:
    type: Directory
    outputSource: subdirs/collapser_dir

  bowtie_out_dir:
    type: Directory
    outputSource: subdirs/bowtie_dir

  counter_out_dir:
    type: Directory
    outputSource: subdirs/counter_dir

  # Differential expression outputs
  diffex_norm_counts:
    type: File
    outputSource: differential-expression/norm_counts

  diffex_comparisons:
    type: File[]
    outputSource: differential-expression/comparisons

  # Plotter outputs
  plots:
    type: File[]
    outputSource: plotter/plots
