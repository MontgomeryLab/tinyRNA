#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: aquatx-merge

inputs:
  input_files:
    type: File[]
    inputBinding:
      prefix: -i
    doc: input files to merge

  sample_names:
    type: string[]
    inputBinding:
      prefix: -s
    doc: sample names associated with input files

  mode:
    type: string
    inputBinding:
      prefix: -m
    doc: mode to run. One of [counts, stats]

  output_file:
    type: string
    inputBinding:
      prefix: -o
    doc: name of the final merged file

outputs:
  merged_file:
    type: File
    outputBinding:
      glob: $(inputs.output_file)
