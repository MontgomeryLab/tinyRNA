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
  sample_basenames: string[]

  # bowtie build
  run_bowtie_build: boolean
  reference_genome_files: File[]
  offrate: int?
  ntoa: boolean?
  noref: boolean?
  ftabchars: int?

  # fastp inputs
  in_fq: File[]
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
  fastp_report_titles: string[]

  # collapser inputs
  threshold: int?
  compress: boolean?

  # bowtie inputs
  bt_index_files: File[]
  ebwt: string
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
  sam: boolean?
  seed: int?
  shared_memory: boolean?

  # counter inputs
  samples_csv: File
  features_csv: File
  gff_files: File[]?
  aligned_seqs: File[]?
  is_pipeline: boolean?
  counter_diags: boolean?
  counter_all_features: boolean?
  counter_type_filter: string[]?
  counter_source_filter: string[]?

  # deseq inputs
  control_condition: string?
  dge_pca_plots: boolean?
  dge_drop_zero: boolean?

  # plotter options
  plot_requests: string[]
  plot_style_sheet: File?

  # output directory names
  dir_name_bt_build: string
  dir_name_fastp: string
  dir_name_collapser: string
  dir_name_bowtie: string
  dir_name_counter: string
  dir_name_dge: string
  dir_name_plotter: string

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
    out: [index_files, console_output]

  counter-prep:
    run: per-library.cwl
    scatter: [in_fq, sample_basename, fastp_report_title]
    scatterMethod: dotproduct
    in:
      sample_basename: sample_basenames
      # fastp
      in_fq: in_fq
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
      fastp_report_title: fastp_report_titles

      # Collapser
      threshold: threshold
      compress: compress

      # Bowtie
      bt_index_files:
        source: [bt_build_optional/index_files, bt_index_files]
        pickValue: first_non_null
        default: bt_index_files  # To appease the workflow validator
      ebwt: ebwt
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
      sam: sam
      threads: threads
      shared_memory: shared_memory
      seed: seed
    out: [fastq_clean, html_report_file, json_report_file, fastp_console, uniq_seqs, uniq_seqs_low, collapser_console,
          aln_seqs, unal_seqs, bowtie_console]
    
  counter-prep-subdirs:
    run: organize-outputs.cwl
    in:
      bt_build_name: dir_name_bt_build
      bt_build_indexes: bt_build_optional/index_files
      bt_build_console: bt_build_optional/console_output
      run_bowtie_build: run_bowtie_build

      fastp_name: dir_name_fastp
      fastp_cleaned_fastq: counter-prep/fastq_clean
      fastp_html_report: counter-prep/html_report_file
      fastp_json_report: counter-prep/json_report_file
      fastp_console: counter-prep/fastp_console

      collapser_name: dir_name_collapser
      collapser_uniq: counter-prep/uniq_seqs
      collapser_low: counter-prep/uniq_seqs_low
      collapser_console: counter-prep/collapser_console

      bowtie_name: dir_name_bowtie
      bowtie_sam: counter-prep/aln_seqs
      bowtie_unal: counter-prep/unal_seqs
      bowtie_console: counter-prep/bowtie_console
    out: [ bt_build_dir, fastp_dir, collapser_dir, bowtie_dir ]

  counter:
    run: ../tools/tiny-count.cwl
    in:
      samples_csv: samples_csv
      config_csv: features_csv
      aligned_seqs: counter-prep/aln_seqs
      gff_files: gff_files
      out_prefix: run_name
      all_features: counter_all_features
      source_filter: counter_source_filter
      type_filter: counter_type_filter
      is_pipeline: {default: true}
      diagnostics: counter_diags
      fastp_logs: counter-prep/json_report_file
      collapsed_fa: counter-prep/uniq_seqs
    out: [feature_counts, other_counts, alignment_stats, summary_stats, intermed_out_files,
          alignment_diags, selection_diags, console_output]

  counter-subdir:
    run: organize-outputs.cwl
    in:
      counter_name: dir_name_counter
      counter_features: counter/feature_counts
      counter_other: counter/other_counts
      counter_alignment_stats: counter/alignment_stats
      counter_summary_stats: counter/summary_stats
      counter_intermed: counter/intermed_out_files
      counter_aln_diag: counter/alignment_diags
      counter_selection_diag: counter/selection_diags
      counter_console: counter/console_output
      features_csv: features_csv
    out: [counter_dir]

  dge:
    run: ../tools/tiny-deseq.cwl
    in:
      input_file: counter/feature_counts
      outfile_prefix: run_name
      control: control_condition
      plots: dge_pca_plots
      drop_zero: dge_drop_zero
    out: [ norm_counts, comparisons, pca_plots, console_output ]

  dge-subdir:
    run: organize-outputs.cwl
    in:
      dge_name: dir_name_dge
      dge_norm: dge/norm_counts
      dge_comparisons: dge/comparisons
      dge_pca: dge/pca_plots
      dge_console: dge/console_output
    out: [dge_dir]

  plotter:
    run: ../tools/tiny-plot.cwl
    in:
      norm_counts: dge/norm_counts
      dge_tables: dge/comparisons
      len_dist: counter/other_counts
      style_sheet: plot_style_sheet
      out_prefix: run_name
      plot_requests: plot_requests
    out: [plots, console_output]

  plotter-subdir:
    run: organize-outputs.cwl
    in:
      plotter_name: dir_name_plotter
      plotter_plots: plotter/plots
      plotter_console: plotter/console_output
    out: [plotter_dir]

outputs:

  # Subdirectory outputs
  bt_build_out_dir:
    type: Directory?
    outputSource: counter-prep-subdirs/bt_build_dir

  fastp_out_dir:
    type: Directory?
    outputSource: counter-prep-subdirs/fastp_dir

  collapser_out_dir:
    type: Directory?
    outputSource: counter-prep-subdirs/collapser_dir

  bowtie_out_dir:
    type: Directory?
    outputSource: counter-prep-subdirs/bowtie_dir

  counter_out_dir:
    type: Directory?
    outputSource: counter-subdir/counter_dir

  dge_out_dir:
    type: Directory?
    outputSource: dge-subdir/dge_dir

  plotter_dir:
    type: Directory?
    outputSource: plotter-subdir/plotter_dir

