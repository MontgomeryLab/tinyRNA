#!/usr/bin/env cwl-runner

######-------------------------------------------------------------------------------######
#
# This Workflow gathers outputs from the preceding counts-prep SubWorkflow into
# corresponding subdirectories. At the WorkflowOutput level within aquatx_wf.cwl,
# these subdirectories are the final output for all non-counts steps.
#
######-------------------------------------------------------------------------------######

cwlVersion: v1.2
class: Workflow

requirements:
  - class: MultipleInputFeatureRequirement

inputs:

  bt_build_indexes: File[]?
  run_bowtie_build: boolean

  fastp_cleaned_fastq: File[]
  fastp_html_report: File[]
  fastp_json_report: File[]

  collapser_uniq: File[]
  collapser_low: File[]?

  bowtie_sam: File[]
  bowtie_unal: File[]?

steps:

  organize_bt_indexes:
    run: ../tools/make-subdir.cwl
    when: $(inputs.run_bowtie_build)
    in:
      run_bowtie_build: run_bowtie_build
      dir_files: bt_build_indexes
      dir_name: { default: "bowtie-build" }
    out: [ subdir ]

  organize_fastp:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ fastp_cleaned_fastq, fastp_html_report, fastp_json_report ]
        linkMerge: merge_flattened
      dir_name: { default: "fastp" }
    out: [ subdir ]

  organize_collapser:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ collapser_uniq, collapser_low ]
        linkMerge: merge_flattened
      dir_name: { default: "collapser" }
    out: [ subdir ]

  organize_bowtie:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [ bowtie_sam, bowtie_unal ]
        linkMerge: merge_flattened
      dir_name: { default: "bowtie" }
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