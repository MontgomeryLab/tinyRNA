#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

inputs:
  # multi input
  threads: int?

  # fastp inputs
  in_fq: File # unscatter
  out_fq: string  # unscatter
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
  json: string # unscatter
  html: string # unscatter
  report_title: string # unscatter

  # collapser inputs
  uniq_seq_prefix: string # unscatter
  threshold: int?
  compress: boolean?

  # bowtie inputs
  bt_index_files: File[]
  ebwt: string
  logfile: string # unscatter
  outfile: string # unscatter
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
  un: string # unscatter
  sam: boolean?
  seed: int?
  shared_memory: boolean?

steps:
  fastp:
    run: ../tools/fastp.cwl
    in:
      thread: threads
      in1: in_fq
      out1: out_fq
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
      json: json
      html: html
      report_title: report_title
    out: [fastq1, report_json, report_html]

  collapse:
    run: ../tools/tiny-collapse.cwl
    in:
      input_file: fastp/fastq1
      out_prefix: uniq_seq_prefix
      threshold: threshold
      compress: compress
    out: [collapsed_fa, low_counts_fa]

  bowtie:
    run: ../tools/bowtie.cwl
    in:
      bt_index_files: bt_index_files
      ebwt: ebwt
      reads: collapse/collapsed_fa
      outfile: outfile
      logfile: logfile
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
      un: un
      sam: sam
      threads: threads
      shared_memory: shared_memory
      seed: seed
    out: [sam_out, unal_seqs, bowtie_log]

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

  uniq_seqs:
    type: File # unscatter
    outputSource: collapse/collapsed_fa

  aln_seqs:
    type: File # unscatter
    outputSource: bowtie/sam_out

  bowtie_log:
    type: File # unscatter
    outputSource: bowtie/bowtie_log

  # Optional outputs
  unal_seqs:
    type: File?
    outputSource: bowtie/unal_seqs

  uniq_seqs_low:
    type: File? # unscatter
    outputSource: collapse/low_counts_fa