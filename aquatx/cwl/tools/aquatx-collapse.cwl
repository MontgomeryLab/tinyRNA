#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: aquatx-collapse

inputs:
  input_file:
    type: File
    inputBinding:
      position: 0
      prefix: -i

  out_prefix:
    type: string
    inputBinding:
      position: 1
      prefix: -o

  threshold:
    type: int?
    default: 0
    inputBinding:
      position: 2
      prefix: -t

  compress:
    type: boolean?
    default: false
    inputBinding:
      position: 3
      prefix: -c

outputs:
  collapsed_fa:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_collapsed.fa*

  low_counts_fa:
    type: File?
    outputBinding:
      glob: $(inputs.out_prefix)_collapsed_lowcounts.fa*
