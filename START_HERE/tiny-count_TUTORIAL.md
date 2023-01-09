# Getting Started With tiny-count

tiny-count is a counting utility that allows for hierarchical assignment of small RNA reads to features based on user-defined selection rules. This tutorial offers an introductory procedure for setting up and running tiny-count using your own data files. 

If you instead want to use the tinyRNA workflow, where tiny-count execution is handled automatically, please see the [other tutorial](tinyRNA_TUTORIAL).

## Installation
Standalone installation requires [conda](https://docs.conda.io/en/main/miniconda.html). If conda is already installed, you can install tiny-count from the bioconda channel. See the [tiny-count installation section in the README](../README.md#tiny-count-standalone-installation) for instructions.

Alternatively, if you have already installed tinyRNA, you can use the `tiny-count` command within the tinyrna conda environment.

## Your Data Files
Gather the following files for the analysis:
1. **SAM files** containing small RNA reads aligned to a reference genome, one file per sample
2. **GFF3 or GFF2/GTF file(s)** containing annotations for features that you want to assign reads to

## Configuration Files
First, you'll need to obtain template copies of the configuration files. Start by activating the conda environment where tiny-count is installed, then run the following command:

```
tiny-count --get-templates
```

Next, fill out the configuration files that were copied:

### 1. The Samples Sheet (samples.csv)
Edit this file to add the paths to your SAM files, and to define the group name, replicate number, etc. for each sample.

### 2. The Paths File (paths.yml)
Edit this file to add the paths to your GFF annotation(s) under the `gff_files` key. You can leave the `alias` key as-is for now. All other keys in this file are used in the tinyRNA workflow.

### 3. The Features Sheet (features.csv)
Edit this file to define the selection rules for assigning reads to features. For now, we'll add a fully permissive rule:

| Select for...  | with value... | Classify as... | Source Filter | Type Filter | Hierarchy | Strand | 5' End Nucleotide | Length | Overlap |
|----------------|---------------|----------------|---------------|-------------|-----------|--------|-------------------|--------|---------|
| Any            | Any           | Any            |               |             | 0         | Both   | Any               | Any    | Partial |

## First Run
Now you're ready to run tiny-count. Make sure you've activated the conda environment where tiny-count is installed, then run the following command:

```
tiny-count --paths-file paths.yml
```

## Outputs
The primary output is feature_counts.csv, a table of classified counts per feature. You can read about the other file outputs in the [Counts and Pipeline Statistics section of the README](../README.md#counts-and-pipeline-statistics).

## Next Steps
Now that you've run tiny-count, you can edit the configuration files to customize the analysis. For example, you can increase the specificity of your selection rule, or add more selection rules with similar or different hierarchy values, or add more GFF files to the Paths File. You can also add more samples to the Samples Sheet, and run tiny-count again to add them to the output.

### What to read next:
- [Feature selection rules and the selection process](../doc/tiny-count.md#feature-selection)
- [GFF aliases in the Features Sheet](../doc/Configuration.md#gff-files)
- [Command line options](../doc/Parameters.md#tiny-count)