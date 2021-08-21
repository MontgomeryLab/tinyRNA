#!/usr/bin/env cwl-runner

######-------------------------------------------------------------------------------######
# This CommandLineTool converts an array of Files into a Directory output
# It allows us to gather multiple files from a workflow step into a single subdirectory
# No external tools are utilized. This is pure CWL.
######-------------------------------------------------------------------------------######

cwlVersion: v1.2
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: InitialWorkDirRequirement
    listing: $(inputs.dir_files)

# This is effectively a no-op shell command
arguments: [":"]

inputs:
  dir_files: File[]
  dir_name: {type: string?, default: "default_dirname"}

outputs:
  subdir:
    type: Directory
    outputBinding:
      # This glob pattern returns the enclosing temp dir
      glob: "."
      # Then we simply rename the enclosing temp dir and return it
      outputEval: |
        ${self[0].basename = inputs.dir_name; return self;}