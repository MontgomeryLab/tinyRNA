#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
 - class: InlineJavascriptRequirement

baseCommand: aquatx-count

inputs:
  input_file:
    type: File
    inputBinding:
      position: 2
      prefix: -i

  ref_annotations:
    type: File[]
    inputBinding:
      position: 2
      prefix: -r

  mask_annotations:
    type: File[]?
    inputBinding:
      position: 2
      prefix: -m

  antisense:
    type: string[]?
    inputBinding:
      position: 2
      prefix: -a

  out_prefix:
    type: string
    inputBinding:
      position: 2
      prefix: -o
      valueFrom: |
        ${ 
          if (self) {
            return self;
          } else {
            return inputs.input_file.basename;
          }
        }

  intermed_file:
    type: boolean?
    inputBinding:
      position: 2 
      prefix: -t

outputs:
  feature_counts:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_out_feature_counts.txt

  other_counts:
    type: File[]
    outputBinding:
      glob:
        - $(inputs.out_prefix)_out_nt_len_dist.csv
        - $(inputs.out_prefix)_out_class_counts.csv

  stats_file:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_stats.txt

  intermed_out_file:
    type: File?
    outputBinding:
      glob: $(inputs.out_prefix)_out_aln_table.txt
