#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

requirements:
  - class: ShellCommandRequirement
  - class: InitialWorkDirRequirement
    listing: $([inputs.dir])


# This is effectively a shell command no-op since CWL is doing all the work here
arguments: [":"]

inputs:
  dir:
    type: Directory

outputs:
  ls_files:
    type: File[]
    outputBinding:
      glob: $(inputs.dir.basename)/*