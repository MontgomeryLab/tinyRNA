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

The overall behavior of the pipeline and its steps is determined by the Run Config file (`run_config.yml`). This YAML file can be edited using a simple text editor. Within it you must specify the location of your Paths file (`paths.yml`). All other settings are optional. [More info](#run-config-details).

#### Paths File

The locations of pipeline file inputs are defined in the Paths file (`paths.yml`). This YAML file includes paths to your Samples and Features Sheets, in addition to your bowtie index prefix (optional) and the final run directory name. The final run directory will contain all pipeline outputs. The directory name is prepended with the `run_name` and current date and time to keep outputs separate. [More info](#paths-file-details).

#### Samples Sheet

To make it simple to specify your fastq files and their locations, along with associated sample names and replicate numbers, we provide a CSV file (`samples.csv`) which can be modified in a spreadsheet editor such as Microsoft Excel or LibreOffice Calc. [More info](#samples-sheet-details).

#### Features Sheet

Small RNAs can often be classified by sequence characteristics, such as length, strandedness, and 5' nucleotide. We provide a Features Sheet (`features.csv`) in which you can define selection rules to more accurately capture counts for the small RNAs of interest. [More info](#features-sheet-details).

#### Plot Stylesheet

Plot styles can be optionally overridden using a matplotlibrc stylesheet. [More info](#plot-stylesheet-details).

## Editing YAML Files
The Run Config and Paths File are YAML formatted files that can be edited with a text editor. Changing values in these files is pretty straight forward, but it is useful to know a little about YAML syntax.

#### Comments
```yaml
# Lines beginning with a # are comments.
# We've used block comments to designate sections
# and to provide documentation locality
```

#### Key-Value Pairs
```yaml
this-is-a-key: "and this is the key's value"
```

#### Empty / Unassigned Values
```yaml
empty-1:
empty-2: ''
empty-3: ""
empty-4: ~
```

#### Lists
```yaml
# Notice: dashes are at the same indentation level; each item begins with dash space
reference_genome_files:
- ../relative/path/genome1.fasta
- /absolute/path/genome2.fasta

# Brackets can be used for more compact lists
source_filter: [source1, source2, source3]
```

#### Strings
```yaml
quoted-string-key: "this is a string"
unquoted-scalar: this is also a valid string
# Beware: unquoted scalars cannot contain a colon followed by a space
```

#### Booleans
```yaml
valid-1: true
valid-2: TrUe
invalid: 'true'  # <-- This is a string, not a boolean
```

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
|  _Column:_ | Input FASTQ Files   | Sample/Group Name | Replicate Number | Control | Normalization |
|-----------:|---------------------|-------------------|------------------|---------|---------------|
| _Example:_ | cond1_rep1.fastq.gz | condition1        | 1                | True    | RPM           |

### Assigning the Control Group
Assigning the control group allows the proper DGE comparisons to be made and plotted. The Control column is where you'll make this indication by writing `true` on any corresponding row. Regardless of the number of replicates in each group, only one associated row needs to have this indication. Do not write `false` or anything else for the other groups; this column should only be used to indicate the affirmative.

### Applying Custom Normalization
Custom normalization can be applied at the conclusion of feature counting using the Normalization column. Unlike the Control column, values in the Normalization column apply to the specific library that they share a row with.

Supported values are:
- **Blank or 1**: no normalization is applied to the corresponding library
- **Any number**: the corresponding library's counts are divided by this number
- **RPM or rpm**: the corresponding library's counts are divided by (its mapped read count / 1,000,000)

>**NOTE**: These normalizations operate independently of Counter's --normalize-by-hits commandline option. The former is concerned with per-library normalization, whereas the latter is concerned with normalization by selected feature count at each locus ([more info](Counter.md#count-normalization)). The commandline option does not enable or disable the normalizations detailed above.

### Low DF Experiments
DESeq2 requires that your experiment design has at least one degree of freedom. If your experiment doesn't include at least one sample group with more than one replicate, DESeq2 will be skipped and DGE related plots will not be produced.

## Features Sheet Details
| _Column:_ | Select for... | with value... | Alias by... | Tag | Hierarchy | Strand | 5' End Nucleotide | Length | Overlap | Feature Source |
|--------------| --- | --- |-------------| --- | --- | --- | --- | --- |---------| --- |
| _Example:_ | Class | CSR | Name        | tag1 | 1 | both | G | 22 | Exact   | ram1.gff3 |

The Features Sheet allows you to define selection rules that determine how features are chosen when multiple features are found overlap an alignment locus. Selected features are "assigned" a portion of the reads associated with the alignment.

Rules apply to features parsed from **all** Feature Sources, with the exception of "Alias by..." which only applies to the Feature Source on the same row. Selection first takes place against feature attributes (GFF column 9), and is directed by defining the attribute you want to be considered (Select for...) and the acceptable values for that attribute (with value...). 

Rules that match features in the first stage of selection will be used in a second stage which performs elimination by hierarchy and interval overlap characteristics. Remaining candidates pass to the third and final stage of selection which examines characteristics of the alignment itself: strand relative to the feature of interest, 5' end nucleotide, and length. 

See [Counter's documentation](Counter.md) for an explanation of each column.

>**Tip**: Don't worry about having duplicate Feature Source entries. Each GFF file is parsed only once.

## Plot Stylesheet Details
Matplotlib uses key-value "rc parameters" to allow for customization of its properties and styles, and one way these parameters can be specified is with a [matplotlibrc file](https://matplotlib.org/3.4.3/tutorials/introductory/customizing.html#a-sample-matplotlibrc-file), which we simply refer to as the Plot Stylesheet. You can obtain a copy of the default stylesheet used by Plotter with the command `tiny get-template`. Please keep in mind that Plotter overrides these defaults for a few specific elements of certain plots. Feel free to reach out if there is a plot style you wish to override but find you are unable to.