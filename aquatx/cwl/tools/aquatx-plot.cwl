#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: aquatx-plot

inputs:
  input_files:
    type: File[]
    inputBinding:
      position: 2
      prefix: -i

  out_prefix:
    type: string?
    inputBinding:
      position: 2
      prefix: -o

  data_types:
    type: string[]
    inputBinding:
      position: 2
      prefix: -d

  references:
    type: File[]?
    inputBinding:
      position: 2
      prefix: -r

  plots:
    type: string[]?
    inputBinding:
      position: 2
      prefix: -p

outputs:
  plots:
    type: File[]
    outputBinding:
      glob: *.pdf

