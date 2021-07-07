#!/usr/bin/env cwl-runner

######-------------------------------------------------------------------------------######
#
# This Workflow gathers outputs from all steps and organizes them into subdirectories.
#
######-------------------------------------------------------------------------------######

cwlVersion: v1.2
class: Workflow

requirements:
  - class: MultipleInputFeatureRequirement

inputs:

  bt_build_name: string
  bt_build_indexes: File[]?
  run_bowtie_build: boolean

  fastp_name: string
  fastp_cleaned_fastq: File[]
  fastp_html_report: File[]
  fastp_json_report: File[]

  collapser_name: string
  collapser_uniq: File[]
  collapser_low: File[]?

  bowtie_name: string
  bowtie_sam: File[]
  bowtie_log: File[]
  bowtie_unal: File[]?

  counter_name: string
  counter_features: File
  counter_other: File[]
  counter_alignment_stats: File
  counter_summary_stats: File
  counter_intermed: File[]?

  dge_name: string
  dge_norm: File
  dge_comparisons: File[]

steps:

  organize_bt_indexes:
    run: ../tools/make-subdir.cwl
    when: $(inputs.run_bowtie_build)
    in:
      run_bowtie_build: run_bowtie_build
      dir_files: bt_build_indexes
      dir_name: bt_build_name
    out: [ subdir ]

  organize_fastp:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ fastp_cleaned_fastq, fastp_html_report, fastp_json_report ]
        linkMerge: merge_flattened
      dir_name: fastp_name
    out: [ subdir ]

  organize_collapser:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ collapser_uniq, collapser_low ]
        linkMerge: merge_flattened
      dir_name: collapser_name
    out: [ subdir ]

  organize_bowtie:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ bowtie_sam, bowtie_unal, bowtie_log ]
        linkMerge: merge_flattened
      dir_name: bowtie_name
    out: [ subdir ]

  organize_counter:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ counter_features, counter_other, counter_alignment_stats, counter_summary_stats, counter_intermed ]
        linkMerge: merge_flattened
      dir_name: counter_name
    out: [ subdir ]

  organize_dge:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ dge_norm, dge_comparisons ]
        linkMerge: merge_flattened
      dir_name: dge_name
    out: [ subdir ]

outputs:

  bt_build_dir:
    type: Directory?
    outputSource: organize_bt_indexes/subdir

  fastp_dir:
    type: Directory
    outputSource: organize_fastp/subdir

  collapser_dir:
    type: Directory
    outputSource: organize_collapser/subdir

  bowtie_dir:
    type: Directory
    outputSource: organize_bowtie/subdir

  counter_dir:
    type: Directory
    outputSource: organize_counter/subdir

  dge_dir:
    type: Directory
    outputSource: organize_dge/subdir