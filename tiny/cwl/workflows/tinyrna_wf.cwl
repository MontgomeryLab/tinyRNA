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

  # tiny-collapse inputs
  threshold: int?
  compress: boolean?
  5p_trim: int?
  3p_trim: int?

  # bowtie inputs
  bt_index_files: File[]
  ebwt: string
  trim5: int?
  trim3: int?
  bt_phred64: boolean?
  solexa: boolean?
  solexa13: boolean?
  end_to_end: int?
  nofw: boolean?
  norc: boolean?
  seedmms: int?
  seedlen: int?
  best: boolean?
  strata: boolean?
  suppress_aln: int?
  k_aln: int?
  all_aln: boolean?
  no_unal: boolean?
  seed: int?
  shared_memory: boolean?

  # tiny-count inputs
  paths_file: File
  samples_csv: File
  features_csv: File
  gff_files: File[]?
  aligned_seqs: File[]?
  is_pipeline: boolean?
  counter_diags: boolean?
  counter_decollapse: boolean?
  counter_stepvector: string?
  counter_all_features: boolean?
  counter_normalize_by_feature_hits: boolean?
  counter_normalize_by_genomic_hits: boolean?

  # tiny-deseq inputs
  run_deseq: boolean
  control_condition: string?
  dge_pca_plot: boolean?
  dge_drop_zero: boolean?

  # tiny-plot options
  plot_requests: string[]
  plot_vector_points: boolean?
  plot_len_dist_min: int?
  plot_len_dist_max: int?
  plot_dge_scatter_min: double?
  plot_dge_scatter_max: double?
  plot_style_sheet: File?
  plot_pval: float?
  plot_unknown_class: string?
  plot_unassigned_class: string?
  plot_class_scatter_filter_include: string[]?
  plot_class_scatter_filter_exclude: string[]?

  # output directory names
  dir_name_bt_build: string
  dir_name_fastp: string
  dir_name_tiny-collapse: string
  dir_name_bowtie: string
  dir_name_tiny-count: string
  dir_name_tiny-deseq: string
  dir_name_tiny-plot: string

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
      threads: threads

      # tiny-collapse
      threshold: threshold
      compress: compress
      5p_trim: 5p_trim
      3p_trim: 3p_trim
    out: [fastq_clean, html_report_file, json_report_file, uniq_seqs, uniq_seqs_low]

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
    out: [ index_files ]

  bowtie:
    run: ../tools/bowtie.cwl
    scatter: [ reads, sample_basename ]
    scatterMethod: dotproduct
    in:
      ebwt: ebwt
      bt_index_files:
        source: [ bt_build_optional/index_files, bt_index_files ]
        pickValue: first_non_null
        default: bt_index_files  # To appease the workflow validator
      reads: preprocessing/uniq_seqs
      sample_basename: sample_basenames
      trim5: trim5
      trim3: trim3
      phred64: bt_phred64
      solexa: solexa
      solexa13: solexa13
      end_to_end: end_to_end
      seedmms: seedmms
      seedlen: seedlen
      nofw: nofw
      norc: norc
      k_aln: k_aln
      all_aln: all_aln
      suppress_aln: suppress_aln
      best: best
      strata: strata
      no_unal: no_unal
      threads: threads
      shared_memory: shared_memory
      seed: seed
    out: [ sam_out, unal_seqs ]

  tiny-count:
    run: ../tools/tiny-count.cwl
    in:
      paths_file: paths_file
      samples_csv: samples_csv
      features_csv: features_csv
      aligned_seqs: bowtie/sam_out
      gff_files: gff_files
      out_prefix: run_name
      all_features: counter_all_features
      normalize_by_feature_hits:
        source: counter_normalize_by_feature_hits
        valueFrom: $(String(self))  # convert boolean -> string
      normalize_by_genomic_hits:
        source: counter_normalize_by_genomic_hits
        valueFrom: $(String(self))  # convert boolean -> string
      decollapse: counter_decollapse
      stepvector: counter_stepvector
      is_pipeline: {default: true}
      diagnostics: counter_diags
      fastp_logs: preprocessing/json_report_file
      collapsed_fa: preprocessing/uniq_seqs
    out: [ feature_counts, rule_counts, norm_counts, mapped_nt_len_dist, assigned_nt_len_dist,
           alignment_stats, summary_stats, decollapsed_sams, alignment_tables,
           assignment_diags, selection_diags, stats_check ]

  tiny-deseq:
    run: ../tools/tiny-deseq.cwl
    when: $(inputs.run_deseq)
    in:
      run_deseq: run_deseq
      input_file: tiny-count/feature_counts
      outfile_prefix: run_name
      control: control_condition
      pca: dge_pca_plot
      drop_zero: dge_drop_zero
    out: [ norm_counts, comparisons, pca_plot ]

  tiny-plot:
    run: ../tools/tiny-plot.cwl
    in:
      norm_counts: tiny-deseq/norm_counts
      dge_tables: tiny-deseq/comparisons
      raw_counts: tiny-count/feature_counts
      rule_counts: tiny-count/rule_counts
      summ_stats: tiny-count/summary_stats
      len_dist_tables:
        source: [ tiny-count/mapped_nt_len_dist, tiny-count/assigned_nt_len_dist ]
        linkMerge: merge_flattened
      len_dist_min:
        source: [ plot_len_dist_min, length_required ]
        pickValue: all_non_null
        valueFrom: |
          $(self.length ? self[0] : null)
      len_dist_max:
        source: [ plot_len_dist_max, length_limit ]
        pickValue: all_non_null
        valueFrom: |
          $(self.length ? self[0] : null)
      dge_min: plot_dge_scatter_min
      dge_max: plot_dge_scatter_max
      unknown_class_label: plot_unknown_class
      unassigned_class_label: plot_unassigned_class
      classes_include: plot_class_scatter_filter_include
      classes_exclude: plot_class_scatter_filter_exclude
      dge_pval: plot_pval
      style_sheet: plot_style_sheet
      out_prefix: run_name
      plot_requests: plot_requests
      vector_scatter: plot_vector_points
    out:
      - len_dist
      - rule_chart
      - class_chart
      - replicate_scatter
      - sample_avg_scatter_by_dge
      - sample_avg_scatter_by_dge_class

  organize_bt_indexes:
    run: ../tools/make-subdir.cwl
    when: $(inputs.run_bowtie_build)
    in:
      run_bowtie_build: run_bowtie_build
      dir_files:
        source: [ bt_build_optional/index_files ]
        pickValue: all_non_null
      dir_name: dir_name_bt_build
    out: [ subdir ]

  organize_fastp:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ preprocessing/fastq_clean, preprocessing/html_report_file, preprocessing/json_report_file ]
      dir_name: dir_name_fastp
    out: [ subdir ]

  organize_tiny-collapse:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ preprocessing/uniq_seqs, preprocessing/uniq_seqs_low ]
      dir_name: dir_name_tiny-collapse
    out: [ subdir ]

  organize_bowtie:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ bowtie/sam_out, bowtie/unal_seqs ]
      dir_name: dir_name_bowtie
    out: [ subdir ]

  organize_tiny-count:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ tiny-count/feature_counts, tiny-count/rule_counts, tiny-count/norm_counts,
                  tiny-count/mapped_nt_len_dist, tiny-count/assigned_nt_len_dist,
                  tiny-count/alignment_stats, tiny-count/summary_stats, tiny-count/decollapsed_sams,
                  tiny-count/assignment_diags, tiny-count/selection_diags, tiny-count/alignment_tables,
                  tiny-count/stats_check, features_csv ]
      dir_name: dir_name_tiny-count
    out: [ subdir ]

  organize_tiny-deseq:
    run: ../tools/make-subdir.cwl
    when: $(inputs.run_deseq)
    in:
      run_deseq: run_deseq
      dir_files:
        source: [ tiny-deseq/norm_counts, tiny-deseq/comparisons ]
        pickValue: all_non_null
      dir_name: dir_name_tiny-deseq
    out: [ subdir ]

  organize_tiny-plot:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [tiny-plot/len_dist, tiny-plot/rule_chart, tiny-plot/class_chart, tiny-plot/replicate_scatter,
                 tiny-plot/sample_avg_scatter_by_dge, tiny-plot/sample_avg_scatter_by_dge_class,
                 tiny-deseq/pca_plot ]
        pickValue: all_non_null
      dir_name: dir_name_tiny-plot
    out: [ subdir ]

outputs:

  # Subdirectory outputs
  bt_build_out_dir:
    type: Directory?
    outputSource: organize_bt_indexes/subdir

  fastp_out_dir:
    type: Directory
    outputSource: organize_fastp/subdir

  tiny-collapse_out_dir:
    type: Directory
    outputSource: organize_tiny-collapse/subdir

  bowtie_out_dir:
    type: Directory
    outputSource: organize_bowtie/subdir

  tiny-count_out_dir:
    type: Directory
    outputSource: organize_tiny-count/subdir

  tiny-deseq_out_dir:
    type: Directory?
    outputSource: organize_tiny-deseq/subdir

  tiny-plot_out_dir:
    type: Directory
    outputSource: organize_tiny-plot/subdir
