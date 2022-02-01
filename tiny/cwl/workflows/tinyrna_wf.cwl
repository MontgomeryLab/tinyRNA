#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow

requirements:
  - class: ScatterFeatureRequirement
  - class: SubworkflowFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

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
  adapter_fasta: File?
  trim_front1: int?
  trim_tail1: int?

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
  counter_decollapse: boolean?
  counter_no_normalize: boolean?
  counter_all_features: boolean?
  counter_type_filter: string[]?
  counter_source_filter: string[]?

  # deseq inputs
  run_deseq: boolean
  control_condition: string?
  dge_pca_plot: boolean?
  dge_drop_zero: boolean?

  # plotter options
  plot_requests: string[]
  plot_vector_points: boolean?
  plot_style_sheet: File?
  plot_pval: float?

  # output directory names
  dir_name_bt_build: string
  dir_name_fastp: string
  dir_name_collapser: string
  dir_name_bowtie: string
  dir_name_counter: string
  dir_name_dge: string
  dir_name_plotter: string

steps:

  preprocessing:
    run: preprocessing.cwl
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
      adapter_fasta: adapter_fasta
      trim_front1: trim_front1
      trim_tail1: trim_tail1

      # Collapser
      threshold: threshold
      compress: compress
    out: [fastq_clean, html_report_file, json_report_file, fastp_console, uniq_seqs, uniq_seqs_low, collapser_console]

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
    out: [ index_files, console_output ]

  bowtie:
    run: ../tools/bowtie.cwl
    scatter: [ reads, sample_basename ]
    scatterMethod: dotproduct
    in:
      reads: preprocessing/uniq_seqs
      sample_basename: sample_basenames
      bt_index_files:
        source: [ bt_build_optional/index_files, bt_index_files ]
        pickValue: first_non_null
        default: bt_index_files  # To appease the workflow validator
      ebwt: ebwt
      outfile: { valueFrom: $(inputs.sample_basename + "_aligned_seqs.sam") }
      logfile: { valueFrom: $(inputs.sample_basename + "_console_output.log") }
      fastq: fastq
      fasta: fasta
      trim5: trim5
      trim3: trim3
      phred64: bt_phred64
      solexa: solexa
      solexa13: solexa13
      end_to_end: end_to_end
      nofw: nofw
      k_aln: k_aln
      all_aln: all_aln
      no_unal: no_unal
      un: { valueFrom: $(inputs.sample_basename + "_unaligned_seqs.fa") }
      sam: sam
      threads: threads
      shared_memory: shared_memory
      seed: seed
    out: [ sam_out, unal_seqs, console_output ]

  counter:
    run: ../tools/tiny-count.cwl
    in:
      samples_csv: samples_csv
      config_csv: features_csv
      aligned_seqs: bowtie/sam_out
      gff_files: gff_files
      out_prefix: run_name
      all_features: counter_all_features
      source_filter: counter_source_filter
      no_normalize: counter_no_normalize
      decollapse: counter_decollapse
      type_filter: counter_type_filter
      is_pipeline: {default: true}
      diagnostics: counter_diags
      fastp_logs: preprocessing/json_report_file
      collapsed_fa: preprocessing/uniq_seqs
    out: [ feature_counts, other_counts, alignment_stats, summary_stats, console_output,
           decollapsed_sams, intermed_out_files, alignment_diags, selection_diags ]

  dge:
    run: ../tools/tiny-deseq.cwl
    when: $(inputs.run_deseq)
    in:
      run_deseq: run_deseq
      input_file: counter/feature_counts
      outfile_prefix: run_name
      control: control_condition
      pca: dge_pca_plot
      drop_zero: dge_drop_zero
    out: [ norm_counts, comparisons, pca_plot, console_output ]

  plotter:
    run: ../tools/tiny-plot.cwl
    in:
      norm_counts: dge/norm_counts
      dge_tables: dge/comparisons
      raw_counts: counter/feature_counts
      summ_stats: counter/summary_stats
      len_dist: counter/other_counts
      dge_pval: plot_pval
      style_sheet: plot_style_sheet
      out_prefix: run_name
      plot_requests: plot_requests
      vector_scatter: plot_vector_points
    out: [plots, console_output]

  organize_bt_indexes:
    run: ../tools/make-subdir.cwl
    when: $(inputs.run_bowtie_build)
    in:
      run_bowtie_build: run_bowtie_build
      dir_files:
        source: [ bt_build_optional/index_files, bt_build_optional/console_output ]
        pickValue: all_non_null
      dir_name: dir_name_bt_build
    out: [ subdir ]

  organize_fastp:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ preprocessing/fastq_clean, preprocessing/html_report_file, preprocessing/json_report_file,
                  preprocessing/fastp_console ]
      dir_name: dir_name_fastp
    out: [ subdir ]

  organize_collapser:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ preprocessing/uniq_seqs, preprocessing/uniq_seqs_low, preprocessing/collapser_console ]
      dir_name: dir_name_collapser
    out: [ subdir ]

  organize_bowtie:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ bowtie/sam_out, bowtie/unal_seqs, bowtie/console_output ]
      dir_name: dir_name_bowtie
    out: [ subdir ]

  organize_counter:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ counter/feature_counts, counter/other_counts, counter/alignment_stats, counter/summary_stats,
                  counter/intermed_out_files, counter/alignment_diags, counter/selection_diags, counter/console_output,
                  counter/decollapsed_sams, features_csv ]
      dir_name: dir_name_counter
    out: [ subdir ]

  organize_dge:
    run: ../tools/make-subdir.cwl
    when: $(inputs.run_deseq)
    in:
      run_deseq: run_deseq
      dir_files:
        source: [ dge/norm_counts, dge/comparisons, dge/console_output ]
        pickValue: all_non_null
      dir_name: dir_name_dge
    out: [ subdir ]

  organize_plotter:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ plotter/plots, plotter/console_output, dge/pca_plot ]
        pickValue: all_non_null
      dir_name: dir_name_plotter
    out: [ subdir ]

outputs:

  # Subdirectory outputs
  bt_build_out_dir:
    type: Directory?
    outputSource: organize_bt_indexes/subdir

  fastp_out_dir:
    type: Directory
    outputSource: organize_fastp/subdir

  collapser_out_dir:
    type: Directory
    outputSource: organize_collapser/subdir

  bowtie_out_dir:
    type: Directory
    outputSource: organize_bowtie/subdir

  counter_out_dir:
    type: Directory
    outputSource: organize_counter/subdir

  dge_out_dir:
    type: Directory?
    outputSource: organize_dge/subdir

  plotter_out_dir:
    type: Directory
    outputSource: organize_plotter/subdir
