#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
 - class: InitialWorkDirRequirement
   listing: $(inputs.bt_index_files)
 - class: InlineJavascriptRequirement

baseCommand: bowtie

inputs:
  # Only inputs relevant to the srna pipeline are listed
  ebwt:
    type: string
    inputBinding:
      position: 2
    doc: "The basename of the index to be searched."
  
  bt_index_files:
    type: File[]
    doc: "Index files for bowtie2 alignment."

  outfile:
    type: string
    inputBinding:
      position: 4
    doc: "File to write hits to"

  reads:
    type: File
    inputBinding:
      itemSeparator: ","
      position: 3
    doc: "Comma-separated list of files containing unpaired reads"

  fastq:
    type: boolean?
    inputBinding:
      prefix: -q
      position: 1
    doc: "query input files are FASTQ .fq/.fastq"

  fasta:
    type: boolean
    inputBinding:
      prefix: -f
      position: 1
    default: true
    doc: "query input files are (multi-)FASTA .fa/.mfa"
  
  trim5:
    type: int?
    inputBinding:
      prefix: --trim5
      position: 1
    doc: "trim <int> bases from 5' (left) end of reads"
  
  trim3: 
    type: int?
    inputBinding:
      prefix: --trim3
      position: 1
    doc: "trim <int> bases from 3' (right) end of reads"

  phred64:
    type: boolean?
    inputBinding:
      prefix: --phred64-quals
      position: 1
    doc: "input quals are Phred+64 (same as --solexa1.3-quals)"

  solexa:
    type: boolean?
    inputBinding:
      prefix: --solexa-quals
      position: 1
    doc: "input quals are from GA Pipeline ver. < 1.3"

  solexa13:
    type: boolean?
    inputBinding:
      prefix: --solexa1.3-quals
      position: 1
    doc: "input quals are from GA Pipeline ver. >= 1.3"

  end_to_end:
    type: int
    inputBinding:
      prefix: -v
      position: 1
    default: 0
    doc: "report end-to-end hits w/ <=v mismatches; ignore qualities"

  nofw:
    type: boolean?
    inputBinding:
      prefix: --nofw
      position: 1
    doc: "do not align to forward reference strand"

  norc:
    type: boolean?
    inputBinding:
      prefix: --norc
      position: 1
    doc: "do not align to reverse-complement reference strand"
  
  # Reporting
  k_aln:
    type: int?
    inputBinding:
      prefix: -k
      position: 1
    doc: "report up to <int> good alignments per read (default: 1)"

  all:
    type: boolean
    inputBinding:
      prefix: --all
      position: 1
    default: true
    doc: "report all alignments per read (much slower than low -k)"

  # Output
  no_unal:
    type: boolean?
    inputBinding:
      prefix: --no-unal
      position: 1
    default: true   
    doc: "suppress SAM records for unaligned reads"

  un:
    type: string?
    inputBinding:
      prefix: --un
    doc: "write unaligned reads/pairs to file(s)"

  # SAM
  sam:
    type: boolean?
    inputBinding:
      prefix: --sam
      position: 1
    default: true
    doc: "write hits in SAM format"
  
  # Performance inputs
  threads:
    type: int?
    inputBinding:
      prefix: --threads
      position: 1
    doc: "number of alignment threads to launch (default: 1)"

  seed:
    type: int?
    inputBinding: 
      prefix: --seed  
      position: 1
    doc: "seed for random number generator"

outputs:
  sam_out:
    type: File
    outputBinding:
      glob: $(inputs.outfile)

  unal_seqs:
    type: File?
    outputBinding:
      glob: |
        ${
          if (inputs.un) {
            return inputs.un;
          } else {
            return [];
          }
        }
