#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
 - class: StepInputExpressionRequirement
 - class: InlineJavascriptRequirement

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
    out: [fastq1, report_json, report_html, console_output]

  collapse:
    run: ../tools/tiny-collapse.cwl
    in:
      input_file: fastp/fastq1
      sample_basename: sample_basename
      out_prefix: { valueFrom: $(inputs.sample_basename) }
      threshold: threshold
      compress: compress
    out: [collapsed_fa, low_counts_fa, console_output]

  bowtie:
    run: ../tools/bowtie.cwl
    in:
      reads: collapse/collapsed_fa
      sample_basename: sample_basename
      bt_index_files: bt_index_files
      ebwt: ebwt
      outfile: {valueFrom: $(inputs.sample_basename + "_aligned_seqs.sam")}
      logfile: {valueFrom: $(inputs.sample_basename + "_console_output.log")}
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
      un: {valueFrom: $(inputs.sample_basename + "_unaligned_seqs.fa")}
      sam: sam
      threads: threads
      shared_memory: shared_memory
      seed: seed
    out: [sam_out, unal_seqs, console_output]

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

  uniq_seqs:
    type: File # unscatter
    outputSource: collapse/collapsed_fa

  collapser_console:
    type: File # unscatter
    outputSource: collapse/console_output

  aln_seqs:
    type: File # unscatter
    outputSource: bowtie/sam_out

  bowtie_console:
    type: File # unscatter
    outputSource: bowtie/console_output

  # Optional outputs
  unal_seqs:
    type: File?
    outputSource: bowtie/unal_seqs

  uniq_seqs_low:
    type: File? # unscatter
    outputSource: collapse/low_counts_fa