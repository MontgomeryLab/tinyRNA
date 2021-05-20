#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
 - class: InitialWorkDirRequirement
   listing: $(inputs.bt_index_files)

baseCommand: bowtie

# Only inputs relevant to the srna pipeline are listed
inputs:

  ebwt:
    type: string
    inputBinding:
      position: 17
    doc: "The basename of the index to be searched."

  # Only used by InitialWorkDirRequirement
  bt_index_files:
    type: File[]
    doc: "Index files for bowtie alignment."

  reads:
    type: File
    inputBinding:
      itemSeparator: ","
      position: 18
    doc: "Comma-separated list of files containing unpaired reads"

  outfile:
    type: string
    inputBinding:
      position: 19
    doc: "File to write hits to"

  fastq:
    type: boolean?
    inputBinding:
      prefix: -q
      position: 0
    doc: "query input files are FASTQ .fq/.fastq"

  fasta:
    type: boolean?
    inputBinding:
      prefix: -f
      position: 1
    default: true
    doc: "query input files are (multi-)FASTA .fa/.mfa"

  trim5:
    type: int?
    inputBinding:
      prefix: --trim5
      position: 2
    doc: "trim <int> bases from 5' (left) end of reads"

  trim3:
    type: int?
    inputBinding:
      prefix: --trim3
      position: 3
    doc: "trim <int> bases from 3' (right) end of reads"

  phred64:
    type: boolean?
    inputBinding:
      prefix: --phred64-quals
      position: 4
    doc: "input quals are Phred+64 (same as --solexa1.3-quals)"

  solexa:
    type: boolean?
    inputBinding:
      prefix: --solexa-quals
      position: 5
    doc: "input quals are from GA Pipeline ver. < 1.3"

  solexa13:
    type: boolean?
    inputBinding:
      prefix: --solexa1.3-quals
      position: 6
    doc: "input quals are from GA Pipeline ver. >= 1.3"

  ### ALIGNMENT ###

  end_to_end:
    type: int?
    inputBinding:
      prefix: -v
      position: 7
    default: 0
    doc: "report end-to-end hits w/ <=v mismatches; ignore qualities"

  nofw:
    type: boolean?
    inputBinding:
      prefix: --nofw
      position: 8
    doc: "do not align to forward/reverse-complement reference strand"

  ### REPORTING ###

  k_aln:
    type: int?
    inputBinding:
      prefix: -k
      position: 9
    doc: "report up to <int> good alignments per read (default: 1)"

  all_aln:
    type: boolean?
    inputBinding:
      prefix: --all
      position: 10
    default: true
    doc: "report all alignments per read (much slower than low -k)"

  ### OUTPUT ###

  time:
    type: boolean?
    inputBinding:
      prefix: -t
      position: 11
    default: true
    doc: "print wall-clock time taken by search phases"

  un:
    type: string?
    inputBinding:
      prefix: --un
      position: 12
    doc: "write unaligned reads/pairs to file(s) <fname>"

  no_unal:
    type: boolean?
    inputBinding:
      prefix: --no-unal
      position: 12
    default: true
    doc: "suppress SAM records for unaligned reads"

  ### SAM ###

  sam:
    type: boolean?
    inputBinding:
      prefix: --sam
      position: 13
    default: true
    doc: "write hits in SAM format"

  ### PERFORMANCE ###

  threads:
    type: int?
    inputBinding:
      prefix: --threads
      position: 14
    doc: "number of alignment threads to launch (default: 1)"

  shared_memory:
    type: boolean?
    inputBinding:
      prefix: --shmem
      position: 15
    doc: "use shared mem for index; many bowtie's can share"

  ### OTHER ###

  seed:
    type: int?
    inputBinding:
      prefix: --seed
      position: 16
    doc: "seed for random number generator"

outputs:
  sam_out:
    type: File
    outputBinding:
      glob: $(inputs.outfile)

  unal_seqs:
    type: File?
    outputBinding:
      glob: $(inputs.un)
