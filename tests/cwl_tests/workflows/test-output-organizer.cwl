#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow

requirements:
  - class: MultipleInputFeatureRequirement

inputs:
  ls_dir: Directory
  outdir_name: string

steps:
  get_files:
    run: ../tools/ls.cwl
    in:
      dir: ls_dir
    out: [ls_files]

  organize_output:
    run: ../../../tiny/cwl/tools/make-subdir.cwl
    in:
      dir_files: get_files/ls_files
      dir_name: outdir_name
    out: [subdir]

outputs:
  type: Directory
  outputSource: organize_output/subdir