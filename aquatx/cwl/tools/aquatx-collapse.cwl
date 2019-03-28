#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: aquatx-collapse

inputs:
  input_file:
    type: File
    inputBinding:
      position: 2
      prefix: -i

  out_file:
    type: string
    inputBinding:
      position: 2
      prefix: -o

  threshold:
    type: int?
    inputBinding:
      position: 2
      prefix: -t

  keep_low_counts:
    type: string?
    inputBinding:
      position: 2
      prefix: -k

outputs:
  collapsed_fa:
    type: File
    outputBinding:
      glob: $(inputs.out_file)

  low_counts_fa:
    type: File?
    outputBinding:
      glob: $(inputs.keep_low_counts)
