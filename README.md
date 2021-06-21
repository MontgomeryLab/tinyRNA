# AQuATx for small RNA sequencing data
:warning: **Under Development & Testing** :warning:

*This repository is being actively developed and tested, and is thus incomplete and only recommended for testing until a release is made. Feedback, suggestions, and bug reports are welcome under the [issues tab](https://github.com/biokcb/aquatx-srna/issues). Thank you!*

[![Build Status](https://travis-ci.com/MontgomeryLab/aquatx-srna.svg?branch=master)](https://travis-ci.com/MontgomeryLab/aquatx-srna)

- [Getting Started](#getting-started)
  - [Prerequisites & Installation](#prerequisites--installation)
    - [1. Install conda](#1-install-conda)
    - [2. Install R & DESeq2](#2-install-r--deseq2)
    - [3. Install AQuATx](#3-install-aquatx)
- [Usage](#usage)
  - [Configuration Files](#configuration-files)
    - [Run Config](#run-config)
    - [Paths File](#paths-file)
    - [Samples Sheet](#samples-sheet)
    - [Features Sheet](#features-sheet)
  - [Input formats](#input-formats)
  - [Running the end-to-end analysis](#running-the-end-to-end-analysis)
  - [Running individual steps](#running-individual-steps)
      - [Create Workflow](#create-workflow)
      - [Collapser](#collapser)
      - [Counter](#counter)
      - [fastp, bowtie-build, and bowtie](#fastp-bowtie-build-and-bowtie)
  - [Using a different workflow runner](#using-a-different-workflow-runner)
- [Outputs](#outputs)
  - [Fastq file quality analysis](#fastq-file-quality-analysis)
  - [Counts and pipeline statistics](#counts-and-pipeline-statistics)
    - [Feature Counts](#feature-counts)
    - [Alignment Statistics](#alignment-statistics)
    - [Summary Statistics](#summary-statistics)
    - [5'nt vs. length matrix](#5nt-vs-length-matrix)
- [Contributing](#contributing)
- [Authors](#authors)
- [License](#license)


AQuATx (Automated QUantitative Analysis of Transcript eXpression) is a set of tools to simplify the analysis of next-generation sequencing data. The goal of this specific repository is to provide an entire workflow for processing small RNA sequencing data with options for advanced hierarchical feature selection.

The current workflow is as follows:

![AQuATx basic pipeline](images/aquatx-workflow_current.png)

Steps outlined in orange are Python scripts we've created for this task, optional steps are teal, and purple steps are part of the end-to-end standard workflow. Outputs are shown in green.

## Getting Started

### Prerequisites & Installation

A conda environment file has been provided for easy installation of the AQuATx environment and its dependencies.

#### 1. Install conda
To install `conda` you can download and follow instructions for either:
- [Anaconda](https://www.anaconda.com/distribution/), which includes many packages, a few of our dependencies, and additional tools that you may find useful in the context of data science, or
- [Miniconda3](https://docs.conda.io/en/latest/miniconda.html), which contains only conda and its dependencies, so that its installation requires significantly less time and disk space.

See the [conda user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda) if you need additional help choosing.

#### 2. Install R & DESeq2
Rather than installing R via conda, we recommend you install it yourself first from [CRAN](https://www.r-project.org/), then install DESeq2 following [instructions](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) in Bioconductor. 

#### 3. Install AQuATx

To install AQuATx and its remaining dependencies:
```
# Clone the repository into a local directory
git clone https://github.com/MontgomeryLab/aquatx-srna.git
cd aquatx-srna

# Install the aquatx-srna environment and dependencies
conda env create -f environment.yml
```

## Usage
If you'd like to jump right in and start using AQuATx, see our [tutorial](./START_HERE/TUTORIAL.md).

You may execute the workflow in its entirety for a full end-to-end analysis pipeline, or you may choose to execute individual steps on their own. In most cases you will use the `aquatx` command for pipeline level operations, including running the pipeline.

### Configuration Files
The pipeline requires that you specify your input library files (Samples Sheet), your selection rules for feature counting (Features Sheet), the paths to your configuration files and other file inputs (Paths), and your preferences for the overall pipeline configuration (Run Config).

![AQuATx basic pipeline](images/config-files.png)

You may obtain a template copy of these files with the command:
```
aquatx get-template
```
:warning: You may use either relative or absolute paths in your configuration files. **Relative paths will be evaluated relative to the file in which they are defined.** This means that you may store your configuration files each in a different location, allowing for flexibility in your project organization. :warning:

#### Run Config

The pipeline revolves around a configuration file to make it easy to set up and run. This `YAML` Run Config file can be edited using a simple text editor. Within it you must define the location of your Paths file, and you may optionally define your preferences for the pipeline and its individual steps. During the setup phase of pipeline execution, AQuATx will further process this configuration file based on its contents and the contents of your Paths, Samples, and Features files. The processed configuration is what ultimately determines the behavior of the workflow. A copy will be saved in the final run directory specified in your Paths file.

#### Paths File

The locations of pipeline file inputs are defined in the Paths file. This `YAML` file includes your Samples and Features Sheets, in addition to your bowtie index prefix (optional) and the final run directory name. The final run directory will contain all pipeline outputs, and is therefore recreated and prepended with the `run_name` and current date and time of each run to keep outputs separate.

#### Samples Sheet

To make it simple to specify your sample files, along with associated sample names and replicate numbers, we provide a `csv` file which can be edited in a spreadsheet editor such as Microsoft Excel or LibreOffice Calc.

#### Features Sheet

sRNA research is often confounded by the fact that loci of features of interest may overlap with the loci of overabundant features, such as rRNA coding genes. These cases can lead to artificial inflation or deflation of feature counts of interest. Thus, we provide a Features Sheet in which you may define selection rules and obtain a more honest representation of a library's sRNA feature counts. These rules may be defined, per GFF file, on the basis of any attribute key/value pair, strand, 5' end nucleotide, and length, with options for full or partial interval matching with the target sequence.

Selection takes place for every feature associated with every alignment of every read. It occurs in two phases:
1. Against the candidate feature's attribute key/value pairs, as defined in your reference annotation files
2. Against the target sequence's attributes (strand, 5' end nucleotide, and length)

Each rule must be assigned a hierarchy value. A lower value indicates higher selection preference and multiple rules may share the same value. We utilize this value only during the first phase of selection; if multiple features match the attribute key/value pairs defined in your rules, then only the feature(s) with the lowest hierarchy values move to the second selection phase. The remaining features are discarded for the given read alignment.

Any candidate features which pass selection will receive a normalized count increment. Read counts are normalized twice before being assigned to a feature. A read's count is divided: 
1. By the number of loci it aligns to in the genome
2. By the number of selected features for each of its alignments.

You may wish for your final feature counts to be expressed in terms of a particular feature attribute key, e.g. by sequence name, and we have provided the option to do so with the `ID Attribute` column. 

### Input formats

The pipeline accepts the following formats:
  1. Reference annotations must be in GFF3
     - Each feature must have an attributes column which defines its `ID` and `Class` (case-sensitive)
     - Each feature's `ID` attribute must be unique
     - All features must be stranded
     - Attribute values which contain commas will be parsed as lists
  2. Library files must be FASTQ(.gz) <sup>*</sup>
  3. Reference genome files must be FASTA
  4. Bowtie indexes must be small indexes (.ebwt)

<sup>*</sup> `aquatx-count` accepts SAM files in your Samples Sheet only when invoked as an individual step. If these entries are present during an end-to-end run, an error will be thrown at pipeline start.

### Running the end-to-end analysis
In most cases you will use this toolset as an end-to-end pipeline. This will run a full, standard small RNA sequencing data analysis according to your configuration files, from raw fastq files to feature counts and summary statistics.

```
aquatx run --config path/to/Run_Config.yml
```

### Running individual steps
The process for running individual steps differs depending on whether the step is an AQuATx Python component, or a CWL wrapped third party tool.

The following steps are Python components. Their corresponding commands may be invoked at the command line:
##### Create Workflow
```
aquatx-config -i CONFIG

  required arguments:
    -i CONFIG, --input-file CONFIG
                          The Run Config file to be processed
```
##### Collapser
```
aquatx-collapse [-h] -i FASTQFILE -o OUTPREFIX [-t THRESHOLD] [-c]

  optional arguments:
    -h, --help            show this help message and exit
    -t THRESHOLD, --threshold THRESHOLD
                          Sequences <= THRESHOLD will be omitted from
                          {prefix}_collapsed.fa and will instead be placed in
                          {prefix}_collapsed_lowcounts.fa
    -c, --compress        Use gzip compression when writing fasta outputs
  
  required arguments:
    -i FASTQFILE, --input-file FASTQFILE
                          The input fastq(.gz) file to collapse
    -o OUTPREFIX, --out-prefix OUTPREFIX
                          The prefix for output files {prefix}_collapsed.fa and,
                          if counts fall below threshold,
                          {prefix}_collapsed_lowcounts.fa
```
##### Counter
```
aquatx-count [-h] -i SAMPLES -c CONFIGFILE -o OUTPUTPREFIX [-t] [-p]

  optional arguments:
    -h, --help            show this help message and exit
    -t, --intermed-file   Save the intermediate file containing all alignments
                          and associated features.
    -p, --is-pipeline     Indicates that counter was invoked from the aquatx
                          pipeline and that input files should be sources as
                          such.
  
  required arguments:
    -i SAMPLES, --input-csv SAMPLES
                          the csv samples file/library list
    -c CONFIGFILE, --config CONFIGFILE
                          the csv features configuration file
    -o OUTPUTPREFIX, --out-prefix OUTPUTPREFIX
                          output prefix to use for file names
```
##### fastp, bowtie-build, and bowtie
These are CWL wrapped third party tools.
1. Copy the workflow CWL folder to your current working directory with the command `aquatx setup-cwl --config none`
2. Within `./CWL/tools` find the file for the step you wish to run. Navigate to this folder in terminal (or copy your target .cwl file to a more convenient location)
3. Run `cwltool --make-template step-file.cwl > step-config.YML`. This will produce a `YML` configuration file specific to this step. Optional arguments will be indicated as such; if you do not wish to set a value for an optional argument, best practice is to remove it from the file
4. Fill in your preferences and inputs in this step configuration file and save it
5. Execute the tool with the command `cwltool step-file.cwl step-config.YML`

### Using a different workflow runner

We have used CWL to define the workflow for scalability and interoperability. The default runner, or interpreter, utilized by AQuATx is `cwltool`. You may use a different CWL runner if you would like, and in order to do so you will need the workflow CWL and your **processed** Run Config file. The following will copy these files to your current working directory:

```
aquatx setup-cwl --config <path/to/Run_Config.yml>
```

If you don't have a Run Config file or do not wish to obtain a processed copy, you may instead use "None" or "none" in the `--config` argument:
```
aquatx setup-cwl --config none
```

## Outputs

The files produced by each pipeline step will be included in the final run directory by default. These intermediate files are organized into subdirectories by step.

### Fastq file quality analysis

[Fastp](https://github.com/OpenGene/fastp) produces a quality filtered FASTQ file along with a summary and quality statistics report file for each library. These report files determine the pipeline summary statistics for total reads and retained reads.

### Counts and pipeline statistics
The Counter step produces final analysis tables for feature counts, alignment statistics, pipeline summary statistics, and 5'nt vs. length matrices for each library.

#### Feature Counts
A single table of feature counts includes columns for each library analyzed. A gene's _Feature ID_ and _Feature Class_ are simply the values of its `ID` and `class` attributes. However, we don't often think of genes by their ID attribute, so we have included a _Feature Name_ column to display a more human friendly alias of your choice. In the Features Sheet you will find the _ID Attribute_ column. Here you may choose another attribute to represent each feature in the _Feature Name_ column of this table.

For example, say your Features Sheet has a rule which specifies an ID Attribute of `sequence_name`, which has a value of "abc123,def456,123.456" for gene1. Let's also say gene1 is both ALG and CSR class. The GFF entry for this feature would have the following attributes column:
```
... ID=gene1;sequence_name=abc123,def456,123.456;Class=ALG,CSR; ...
```
The row for this feature in the feature counts table would read:

| Feature ID | Feature Name | Feature Class | Group1_rep_1 | Group1_rep_2 | ... |
|------------|--------------|---------------|-----------|-----------|-----|
| gene1 | abc123,def456,123.456 | ALG, CSR | 1234 | 999 | ... |


#### Alignment Statistics
A single table of alignment statistics includes columns for each library and the following rows:
- Total Assigned Reads
- Total Assigned Sequences
- Assigned Single-Mapping Reads
- Assigned Multi-Mapping Reads
- Reads Assigned to Single Feature
- Sequences Assigned to Single Feature
- Reads Assigned to Multiple Features
- Sequences Assigned to Multiple Features
- Total Unassigned Reads
- Total Unassigned Sequences

#### Summary Statistics
A single table of summary statistics includes columns for each library and the following rows:
- Total Reads
- Retained Reads
- Unique Sequences
- Mapped Sequences
- Aligned Reads

#### 5'nt vs. length matrix

After alignment, a size and 5'nt distribution table is created for each library. The distribution of lengths and 5'nt can be used to assess the overall quality of your libraries. This can also be used for analyzing small RNA distributions in non-model organisms without annotations.

## Contributing

See the [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines. To see what is actively being worked or planned go to the [projects tab](https://github.com/MontgomeryLab/aquatx-srna/projects) or the [issues tab](https://github.com/MontgomeryLab/aquatx-srna/issues).

## Authors

* **Kristen Brown** - 2018-2019 - Colorado State University - [biokcb](https://github.com/biokcb)
* **Alex Tate** - 01/2021-present - Colorado State University - [AlexTate](https://github.com/AlexTate)

See also the list of [contributors](https://github.com/MontgomeryLab/aquatx-srna/contributors) who participated in this project.

## License

This project is licensed under the GPLv3 license (along with HTSeq, bowtie). License - see the [LICENSE.md](LICENSE.md) file for details
