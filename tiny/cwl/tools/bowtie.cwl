#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
 - class: InitialWorkDirRequirement
   listing: $(inputs.bt_index_files)

baseCommand: bowtie
stdout: $(inputs.reads.basename + "_console_output.log")
stderr: $(inputs.reads.basename + "_console_output.log")

# Only inputs relevant to the srna pipeline are listed
inputs:

  ebwt:
    type: string
    inputBinding:
      prefix: -x
      position: 23
    doc: "The basename of the index to be searched."

  # Only used by InitialWorkDirRequirement
  bt_index_files:
    type: File[]
    doc: "Index files for bowtie alignment."

  reads:
    type: File
    inputBinding:
      position: 24
    doc: "File containing unpaired reads"

  outfile:
    type: string
    inputBinding:
      position: 25
    doc: "File to write hits to"

  logfile:
    type: string
    doc: "File to write Bowtie's stdout and stderr streams to"

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
    doc: "do not align to forward reference strand"

  norc:
    type: boolean?
    inputBinding:
      prefix: --norc
      position: 9
    doc: "do not align to reverse-complement reference strand"

  seedmms:
    type: int?
    inputBinding:
      prefix: --seedmms
      position: 10
    doc: "max mismatches in seed (can be 0-3, default: -n 2)"

  seedlen:
    type: int?
    inputBinding:
      prefix: --seedlen
      position: 11
    doc: "seed length for --seedmms (default: 28)"

  ### REPORTING ###

  best:
    type: boolean?
    inputBinding:
      prefix: --best
      position: 12
    doc: "hits guaranteed best stratum; ties broken by quality"

  strata:
    type: boolean?
    inputBinding:
      prefix: --strata
      position: 13
    doc: "hits in sub-optimal strata aren't reported (requires --best)"

  k_aln:
    type: int?
    inputBinding:
      prefix: -k
      position: 14
    doc: "report up to <int> good alignments per read (default: 1)"

  all_aln:
    type: boolean?
    inputBinding:
      prefix: --all
      position: 15
    default: true
    doc: "report all alignments per read (much slower than low -k)"

  ### OUTPUT ###

  time:
    type: boolean?
    inputBinding:
      prefix: -t
      position: 16
    default: true
    doc: "print wall-clock time taken by search phases"

  un:
    type: string?
    inputBinding:
      prefix: --un
      position: 17
    doc: "write unaligned reads/pairs to file(s) <fname>"

  no_unal:
    type: boolean?
    inputBinding:
      prefix: --no-unal
      position: 18
    default: true
    doc: "suppress SAM records for unaligned reads"

  ### SAM ###

  sam:
    type: boolean?
    inputBinding:
      prefix: --sam
      position: 19
    default: true
    doc: "write hits in SAM format"

  ### PERFORMANCE ###

  threads:
    type: int?
    inputBinding:
      prefix: --threads
      position: 20
    doc: "number of alignment threads to launch (default: 1)"

  shared_memory:
    type: boolean?
    inputBinding:
      prefix: --shmem
      position: 21
    doc: "use shared mem for index; many bowtie's can share"

  ### OTHER ###

  seed:
    type: int?
    inputBinding:
      prefix: --seed
      position: 22
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

  console_output:
    type: stdout
