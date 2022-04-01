#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: tiny-collapse
stdout: $(inputs.input_file.basename + "_console_output.log")

inputs:
  # Fastq files
  input_file:
    type: File
    inputBinding:
      position: 0
      prefix: -i
    doc: "The optionally gzipped fastq files to collapse"

  # Collapsed fasta name
  out_prefix:
    type: string
    inputBinding:
      position: 1
      prefix: -o
    doc: "The prefix for output files {prefix}_collapsed.fa and, if
      counts fall below threshold, {prefix}_collapsed_lowcounts.fa"

  # Count filtering
  threshold:
    type: int?
    default: 0
    inputBinding:
      position: 2
      prefix: -t
    doc: "Sequences <= THRESHOLD will be omitted from {prefix}_collapsed.fa
      and will instead be placed in {prefix}_collapsed_lowcounts.fa"

  # Gzip outputs
  compress:
    type: boolean?
    default: false
    inputBinding:
      position: 3
      prefix: -c
    doc: "Use gzip compression when writing fasta outputs"

  5p_trim:
    type: int?
    default: 0
    inputBinding:
      position: 4
      prefix: --5p-trim
    doc: "Trim the specified number of bases from the 5' end of each sequence"

  3p_trim:
    type: int?
    default: 0
    inputBinding:
      position: 5
      prefix: --3p-trim
    doc: "Trim the specified number of bases from the 3' end of each sequence"

outputs:
  collapsed_fa:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)_collapsed.fa*

  low_counts_fa:
    type: File?
    outputBinding:
      glob: $(inputs.out_prefix)_collapsed_lowcounts.fa*

  console_output:
      type: stdout