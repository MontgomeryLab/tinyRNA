#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: | 
  CWL wrapper for fastp - a tool designed to provide fast all-in-one 
  preprocessing for FastQ files. This tool is developed in C++ with 
  multithreading supported to afford high performance.

requirements:
 - class: InlineJavascriptRequirement

baseCommand: fastp

inputs:
  # File I/O options
  in1: 
    type: File
    inputBinding:
      position: 1
      prefix: --in1
    doc: |
      read1 input file name (string)

  out1:
    type: string
    inputBinding:
      position: 2
      prefix: --out1
    doc: |
      read1 output file name (string [=])

  phred64:
    type: boolean?
    inputBinding:
      prefix: --phred64
    doc: |
       indicate the input is using phred64 scoring
       (it'll be converted to phred33, so the output
       will still be phred33)
 
  compression:
    type: int?
    default: 4
    inputBinding:
      prefix: --compression
    doc: |
      compression level for gzip output (1 ~ 9). 
      1 is fastest, 9 is smallest, default is 4.
      (int [=4])

  dont_overwrite:
    type: boolean?
    inputBinding:
      prefix: --dont_overwrite
    doc: |
      don't overwrite existing files. Overwritting is
      allowed by default.
  
  # Adapter trimming options
  disable_adapter_trimming:
    type: boolean?
    inputBinding:
      prefix: --disable_adapter_trimming 
    doc: |
      adapter trimming is enabled by default. If this
      option is specified, adapter trimming is disabled

  adapter_sequence:
    type: string?
    default: 'auto'
    inputBinding:
      prefix: --adapter_sequence
    doc: |
      the adapter for read1. For SE data, if not specified,
      the adapter will be auto-detected. For PE data, this
      is used if R1/R2 are found not overlapped. 
      (string [=auto])
  
  # Trimming options regardless of quality
  trim_poly_x:
    type: boolean?
    inputBinding:
      prefix: --trim_poly_x
    doc: |
      enable polyX trimming in 3' ends.

  poly_x_min_len:
    type: int?
    inputBinding:
      prefix: --poly_x_min_len
    doc: |
      the minimum length to detect polyX in the read tail.
      10 by default. (int [=10]) 

  # Quality score trimming and filtering options
  disable_quality_filtering:
    type: boolean?
    inputBinding:
      prefix: --disable_quality_filtering
    doc: |
      quality filtering is enabled by default. If this
      option is specified, quality filtering is disabled

  qualified_quality_phred:
    type: int?
    default: 15
    inputBinding:
      prefix: --qualified_quality_phred
    doc: |
      the quality value that a base is qualified. Default
      15 means phred quality >=Q15 is qualified. (int [=15])

  unqualified_percent_limit:
    type: int?
    default: 0
    inputBinding:
      prefix: --unqualified_percent_limit
    doc: |
      how many percents of bases are allowed to be unqualified
      (0~100). Default 40 means 40% (int [=40])

  n_base_limit:
    type: int?
    default: 1
    inputBinding:
      prefix: --n_base_limit
    doc: |
      if one read's number of N base is >n_base_limit, then
      this read/pair is discarded. Default is 5 (int [=5])

  # Sequence length filtering
  disable_length_filtering:
    type: boolean?
    inputBinding:
      prefix: --disable_length_filtering
    doc: |
      length filtering is enabled by default. If this option
      is specified, length filtering is disabled
  
  length_required:
    type: int?
    default: 15
    inputBinding:
      prefix: --length_required
    doc: |
      reads shorter than length_required will be discarded,
      default is 15. (int [=15]) 

  length_limit:
    type: int?
    default: 30
    inputBinding:
      prefix: --length_limit
    doc: |
      reads longer than length_limit will be discarded,
      default 0 means no limitation. (int [=0])
  
  # Over-representation options
  overrepresentation_analysis:
    type: boolean?
    inputBinding:
      prefix: --overrepresentation_analysis
    doc: |
      enable overrepresented sequence analysis.

  overrepresentation_sampling:
    type: int?
    inputBinding:
      prefix: --overrepresentation_sampling
    doc: |
      one in (--overrepresentation_sampling) reads will be computed
      for overrepresentation analysis (1~10000), smaller is slower, 
      default is 20. (int [=20]) 
  
  # Output report options
  json:
    type: string
    inputBinding:
      prefix: --json
    doc: |
      the json format report file name (string [=fastp.json])

  html:
    type: string
    inputBinding:
      prefix: --html
    doc: |
      the html format report file name (string [=fastp.html]) 

  report_title: 
    type: string?
    inputBinding:
      prefix: --report_title
    doc: |
      should be quoted with '' or "", default is "fastp report" 
      (string [=fastp report])
  
  # Parallel processing options
  thread:
    type: int?
    inputBinding:
      prefix: --thread
    default: 2
    doc: |
      worker thread number, default is 2 (int [=2])
  
outputs:
  fastq1:
    type: File
    outputBinding:
      glob: $(inputs.out1)
  
  report_json:
    type: File
    outputBinding:
      glob: $(inputs.json)
  
  report_html:
    type: File
    outputBinding:
      glob: $(inputs.html)
