Shortcuts:
- [Run Config](#run-config-details)
  - [The Processed Run Config](#the-processed-run-config)
- [Paths File](#paths-file-details)
  - [Building Bowtie Indexes](#building-bowtie-indexes)
- [Samples Sheet](#samples-sheet-details)
  - [Assigning the Control Group](#assigning-the-control-group)
  - [Applying Custom Normalization](#applying-custom-normalization)
- [Features Sheet](#features-sheet-details)

# Configuration Files
The pipeline requires that you identify:
- Your preferences for the pipeline and its steps via the *Run Config*
- The location of your config files and other file inputs via the *Paths File*
- Your selection preferences for feature counting via the *Features Sheet*
- Your samples via the *Samples Sheet*

<img src="../images/config-files.png" width="50%" alt="[Configuration File Diagram]">

The `START_HERE` directory demonstrates a working configuration using these files. You can also get a copy of them (and other optional template files) with:
```
tiny get-template
```

## Overview

>**Tip**: Each of the following will allow you to map out paths to your input files for analysis. You can use either relative or absolute paths to do so. **Relative paths will be evaluated relative to the file in which they are defined.** This allows you to flexibly organize and share configurations between projects.

#### Run Config

The overall behavior of the pipeline and its steps is determined by the Run Config file (`run_config.yml`). This YAML file can be edited using a simple text editor. Within it you must specify the location of your Paths file (`paths.yml`). All other settings are optional.

#### Paths File

The locations of pipeline file inputs are defined in the Paths file (`paths.yml`). This YAML file includes paths to your Samples and Features Sheets, in addition to your bowtie index prefix (optional) and the final run directory name. The final run directory will contain all pipeline outputs. The directory name is prepended with the `run_name` and current date and time to keep outputs separate.

#### Samples Sheet

To make it simple to specify your fastq files and their locations, along with associated sample names and replicate numbers, we provide a CSV file (`samples.csv`) which can be modified in a spreadsheet editor such as Microsoft Excel or LibreOffice Calc.

#### Features Sheet

Small RNAs can often be classified by sequence characteristics, such as length, strandedness, and 5' nucleotide. We provide a Features Sheet (`features.csv`) in which you can define selection rules to more accurately capture counts for the small RNAs of interest.

#### Plot Stylesheet

Matplotlib styles can be optionally overridden using a matplotlib stylesheet. The stylesheet provided by `tiny get-template` contains the defaults utilized by Plotter. Please keep in mind that Plotter overrides these defaults for a few specific elements of certain plots. Feel free to reach out if there is a plot style you wish to override but find you are unable to.

## Run Config Details

### The processed Run Config
When the pipeline starts up, tinyRNA will process the Run Config based on the contents of it and your other configuration files, and the processed copy will be saved to your run directory. The processed configuration is what ultimately determines the behavior of the workflow. This provides auto-documentation for all of your pipeline runs.

## Paths File Details

### Building Bowtie Indexes
If you don't have bowtie indexes already built for your reference genome, tinyRNA can build them for you at the beginning of an end-to-end run and reuse them on subsequent runs with the same Paths File.

To build bowtie indexes:
1. Open your Run Config in a text editor and find the `run_bowtie_build` key. Set its value to `true` and save it.
2. Open your Paths File in a text editor and find the `reference_genome_files` key. Add your reference genome file(s) under this key, one per line with a `- ` in front.
3. Find the `ebwt` key and delete its value.
4. Execute an end-to-end pipeline run.

Once your indexes have been built, your Paths File will be modified such that `ebwt` points to their location (prefix) within your Run Directory. This means that indexes will not be unnecessarily rebuilt on subsequent runs as long as the same Paths File is used. If you need them rebuilt, simply repeat steps 3 and 4 above.

## Samples Sheet Details

### Assigning the Control Group
Assigning the control group allows the proper DGE comparisons to be made and plotted. The Control column is where you'll make this indication by writing `true` on any corresponding row. Regardless of the number of replicates in each group, only one associated row needs to have this indication. Do not write `false` or anything else for the other groups; this column should only be used to indicate the affirmative.

### Applying Custom Normalization
Custom normalization can be applied at the conclusion of feature counting using the Normalization column. Unlike the Control column, values in the Normalization column apply to the specific library that they share a row with.

Supported values are:
- **Blank or 1**: no normalization is applied to the corresponding library
- **Any number**: the corresponding library's counts are divided by this number
- **RPM or rpm**: the corresponding library's counts are divided by (its mapped read count / 1,000,000)

>**NOTE**: These normalizations operate independently of Counter's --normalize-by-hits commandline option. The former is concerned with per-library normalization, whereas the latter is concerned with normalization by selected feature count at each locus ([more info](Counter.md#count-normalization)). The commandline option does not enable or disable the normalizations detailed above.

## Features Sheet Details
The Features Sheet allows you to define selection rules that determine how features are chosen when multiple features are found overlap an alignment locus. Selected features are "assigned" a portion of the reads associated with the alignment.

Rules apply to features parsed from **all** Feature Sources, with the exception of "Alias by..." which only applies to the Feature Source on the same row. Selection first takes place against feature attributes (GFF column 9), and is directed by defining the attribute you want to be considered (Select for...) and the acceptable values for that attribute (with value...). 

Rules that match features in the first stage of selection will be used in a second stage which performs elimination by hierarchy and interval overlap characteristics. Remaining candidates pass to the third and final stage of selection which examines characteristics of the alignment itself: strand relative to the feature of interest, 5' end nucleotide, and length. 

Since this configuration file is primarily used by Counter, we recommend you see [Counter's documentation](Counter.md) for a more thorough explanation.

>**Tip**: Don't worry about having duplicate Feature Source entries. Each GFF file is parsed only once.