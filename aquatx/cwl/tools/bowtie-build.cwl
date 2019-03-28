#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: bowtie2-build

inputs:
  # references
  ref_in:
    type: File[]
    inputBinding:
      position: 2

  ebwt_base:
    type: string
    inputBinding:
      position: 3

  # options
  fasta:
    type: boolean
    inputBinding: 
      position: 1
      prefix: -f
    default: true

  ftabchars:
    type: int?
    inputBinding:
      position: 1
      prefix: -t
    doc: num of chars consumed in initial lookup (default 10)
  
  ntoa:
    type: boolean?
    inputBinding: 
      position: 1
      prefix: --ntoa
    doc: convert Ns in reference to As
  
  noref:
    type: boolean?
    inputBinding: 
      position: 1
      prefix: -r

  offrate:
    type: int?
    inputBinding: 
      position: 1
      prefix: -o

  seed:
    type: int?
    inputBinding: 
      position: 1
      prefix: --seed

  threads:
    type: int?
    inputBinding: 
      position: 1
      prefix: --threads

outputs:
  index_files:
    type: File[]
    outputBinding:
      glob: $(inputs.ebwt_base).*.ebwt
    
