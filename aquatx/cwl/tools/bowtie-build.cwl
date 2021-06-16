#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: bowtie-build

inputs:

  ### REFERENCES ###

  ref_in:
    type: File[]
    inputBinding:
      itemSeparator: ","
      position: 2
    doc: "comma-separated list of files with ref sequences"

  ebwt_base:
    type: string
    inputBinding:
      position: 3
    doc: "write Ebwt data to files with this dir/basename"

  ### OPTIONS ###

  fasta:
    type: boolean?
    default: true
    inputBinding:
      position: 1
      prefix: -f
    doc: "reference files are fasta format"

  noref:
    type: boolean?
    default: false
    inputBinding:
      position: 12
      prefix: --noref
    doc: "don't build .3/.4 index files"

  offrate:
    type: int?
    default: 5
    inputBinding:
      position: 14
      prefix: --offrate
    doc: "SA is sampled every 2^<int> BWT chars (default: 5)"

  ftabchars:
    type: int?
    default: 10
    inputBinding:
      position: 15
      prefix: --ftabchars
    doc: "# of chars consumed in initial lookup (default: 10)"

  ntoa:
    type: boolean?
    default: false
    inputBinding:
      position: 16
      prefix: --ntoa
    doc: "convert Ns in reference to As"

  seed:
    type: int?
    inputBinding:
      position: 17
      prefix: --seed
    doc: "seed for random number generator"

outputs:
  index_files:
    type: File[]
    outputBinding:
      glob: $(inputs.ebwt_base).*.ebwt
