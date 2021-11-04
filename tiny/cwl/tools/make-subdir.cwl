#!/usr/bin/env cwl-runner

######-------------------------------------------------------------------------------######
# This ExpressionTool returns a directory containing the files/directories provided.
# The dir_files input is extremely flexible. It can contain a single file or directory,
# or lists of files and/or directories, any of which may be or contain null values. Nested
# lists are unpacked recursively. Directories are added as themselves rather than their
# individual file listings.
######-------------------------------------------------------------------------------######

cwlVersion: v1.2
class: ExpressionTool

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing: $(inputs.dir_files)

inputs:
  dir_files: Any?
  dir_name: {type: string?, default: "default_dirname"}

outputs:
  subdir: Directory

expression: |
  ${
    function add_array(files){
      var listing = [];
      for (var i in files){
        var item = files[i]
        if (item == null) continue;
        if (Array.isArray(item)) listing.push.apply(listing, add_array(item));
        if (["File", "Directory"].includes(item["class"])) listing.push(item);
      }
      return listing;
    }

    return {"subdir": {
      "class": "Directory",
      "basename": inputs.dir_name,
      "listing": add_array([inputs.dir_files])
    }}
  }