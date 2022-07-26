#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow

requirements:
 - class: StepInputExpressionRequirement
 - class: InlineJavascriptRequirement
 - class: MultipleInputFeatureRequirement

inputs:
  # multi input
  threads: int?
  sample_basename: string # unscatter

  # fastp inputs
  in_fq: File # unscatter
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
  fastp_report_title: string # unscatter
  adapter_fasta: File?
  trim_front1: int?
  trim_tail1: int?

  # collapser inputs
  run_collapser: boolean
  threshold: int?
  compress: boolean?
  5p_trim: int?
  3p_trim: int?

steps:
  fastp:
    run: ../tools/fastp.cwl
    in:
      thread: threads
      in1: in_fq
      sample_basename: sample_basename
      out1: {valueFrom: $(inputs.sample_basename + "_cleaned.fastq")}
      json: {valueFrom: $(inputs.sample_basename + "_qc.json")}
      html: {valueFrom: $(inputs.sample_basename + "_qc.html")}
      phred64: fp_phred64
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
      report_title: fastp_report_title
      adapter_fasta: adapter_fasta
      trim_front1: trim_front1
      trim_tail1: trim_tail1
    out: [fastq1, report_json, report_html, console_output]

  collapse:
    run: ../tools/tiny-collapse.cwl
    when: $(inputs.run_collapser)
    in:
      run_collapser: run_collapser
      input_file: fastp/fastq1
      sample_basename: sample_basename
      out_prefix: { valueFrom: $(inputs.sample_basename) }
      threshold: threshold
      compress: compress
      5p_trim: 5p_trim
      3p_trim: 3p_trim
    out: [collapsed_fa, low_counts_fa, console_output]

outputs:

  fastq_clean:
    type: File # unscatter
    outputSource: fastp/fastq1

  html_report_file:
    type: File # unscatter
    outputSource: fastp/report_html

  json_report_file:
    type: File # unscatter
    outputSource: fastp/report_json

  fastp_console:
    type: File # unscatter
    outputSource: fastp/console_output

  # Optional outputs
  uniq_seqs:
    type: File? # unscatter
    outputSource: collapse/collapsed_fa

  collapser_console:
    type: File? # unscatter
    outputSource: collapse/console_output

  uniq_seqs_low:
    type: File? # unscatter
    outputSource: collapse/low_counts_fa