# tinyRNA: comprehensive analysis of small RNA high-throughput sequencing data
:warning: **Under Development & Testing** :warning:

*This repository is being actively developed and tested, and is thus incomplete and only recommended for testing until a release is made. Feedback, suggestions, and bug reports are welcome under the [issues tab](https://github.com/MontgomeryLab/tinyrna/issues). See our [changelog](CHANGELOG.md) for recent updates. Thank you!*

- [Installation](#installation)
- [Usage](#usage)
  - [Configuration Files](#configuration-files)
    - [Run Config](#run-config)
    - [Paths File](#paths-file)
    - [Samples Sheet](#samples-sheet)
    - [Features Sheet](#features-sheet)
  - [User-Provided Input Files](#user-provided-input-file-requirements)
  - [Running the End-to-End Analysis](#running-the-end-to-end-analysis)
  - [Running Individual Steps](#running-individual-steps)
  - [Using a Different Workflow Runner](#using-a-different-workflow-runner)
- [Outputs](#outputs)
  - [Data Pre-Processing](#data-pre-processing)
  - [Collapsed FASTA files](#collapsed-fasta-files)
  - [Counts and Pipeline Statistics](#counts-and-pipeline-statistics)
    - [Feature Counts](#feature-counts)
    - [Alignment Statistics](#alignment-statistics)
    - [Summary Statistics](#summary-statistics)
    - [5'nt vs. Length Matrix](#5nt-vs-length-matrix)
  - [Differential Expression Analysis](#differential-expression-analysis)
  - [Plots](#plots)
- [Contributing](#contributing)
- [Authors](#authors)
- [License](#license)


tinyRNA is a set of tools to simplify the analysis of next-generation sequencing data. The goal of this specific repository is to provide an entire workflow for processing small RNA sequencing data with options for advanced hierarchical feature selection.

The current workflow is as follows:

![tinyRNA basic pipeline](images/tinyrna-workflow_current.png)

## Installation

A setup script has been provided for easy installation of tinyRNA. The project and its dependencies will be installed in a conda environment called `tinyrna`.

```
# Clone the repository into a local directory
 git clone https://github.com/MontgomeryLab/tinyrna.git
 cd tinyrna

# Install the tinyrna environment and dependencies
 ./setup.sh

# Activate the tinyrna environment
 conda activate tinyrna

# When you are done running tinyRNA, you can deactivate the conda environment
 conda deactivate
```

If the installation script runs the Miniconda installer:
- Press "q" if you find yourself trapped on the license page
- We recommend answering "yes" to running `conda init`

## Usage
If you'd like to jump right in and start using tinyRNA, see our [tutorial](./START_HERE/TUTORIAL.md).

You can execute the workflow in its entirety for a full end-to-end analysis pipeline, or you can execute individual steps on their own. In most cases you will use the command `tiny` for pipeline level operations, including running the pipeline.

### Configuration Files
The pipeline requires that you identify:
- Your samples via the *Samples Sheet*
- Your selection preferences for feature counting via the *Features Sheet*
- The location of your config files and other file inputs via the *Paths File*
- Your preferences for the pipeline and its steps via the *Run Config*

<img src="images/config-files.png" width="50%" alt="[Configuration File Diagram]">

The `START_HERE` directory demonstrates a working configuration using these files. You can also get a copy of them (and other optional template files) with:
```
tiny get-template
```
You can use either relative or absolute paths in your configuration files. **Relative paths will be evaluated relative to the file in which they are defined.** This gives you the freedom to organize your project as you see fit.

#### Run Config

The overall behavior of the pipeline and its steps is determined by the Run Config file (`run_config.yml`). This YAML file can be edited using a simple text editor such as BBEDIT for Mac or Notepad++ for Windows. Within it you must specify the location of your Paths file (`paths.yml`). All other settings are optional. When the pipeline starts up, tinyRNA will process the Run Config based on the contents of it and your other configuration files, and the processed copy will be saved to your run directory. The processed configuration is what ultimately determines the behavior of the workflow. This provides auto-documentation for all of your pipeline runs.

#### Paths File

The locations of pipeline file inputs are defined in the Paths file (`paths.yml`). This YAML file includes paths to your Samples and Features Sheets, in addition to your bowtie index prefix (optional) and the final run directory name. The final run directory will contain all pipeline outputs. The directory name is prepended with the `run_name` and current date and time to keep outputs separate.

#### Samples Sheet

To make it simple to specify your fastq files and their locations, along with associated sample names and replicate numbers, we provide a `csv` file (`samples.csv`) which can be edited in a spreadsheet editor such as Microsoft Excel or LibreOffice Calc.

#### Features Sheet

Small RNAs can often be classified by sequence characteristics, such as length, strandedness, and 5' nucleotide. We provide a Features Sheet (`features.csv`) in which you can define selection rules to more accurately capture counts for the small RNAs of interest.  Rules apply to features parsed from **all** Feature Sources, with the exception of "Alias by..." which only applies to the Feature Source on the same row. These rules are used to select among overlapping features at each alignment locus. Selection first takes place against feature attributes (GFF3 column 9), and is directed by defining the attribute you want to be considered (Select by...) and the acceptable values for that attribute (with value...). Rules that match features at this stage will undergo hierarchical elimination and pass to the second stage of selection which examines characteristics of the alignment itself: strand relative to the feature of interest, 5' end nucleotide, and length, with options for full or partial interval matching with the target sequence. See [Counter](docs/Counter.md) for more information.

>**Tip**: Don't worry about having duplicate Feature Source entries. Each GFF3 file is parsed only once.

### User-Provided Input File Requirements

  1. A GFF3 formatted file with with genomic coordinates for your features of interest, such as miRNAs (see c_elegans_WS279_chr1.gff3 in the reference_data folder for an example)
     - Each feature must have an attributes column (column 9) which defines its `ID` and `Class` (case-sensitive).
       - For example: `chrI	.	miRNA	100	121	.	+	.	ID=miR-1;Class=miRNA`
     - All features must be stranded.
     - Attribute values which contain commas will be parsed as lists.
  2. FASTQ(.gz) <sup>*</sup> formatted files with your small RNA high-throughput sequencing data (files must be demultiplexed).
  3. A reference genome file in FASTA format (be sure that chromosome identifiers are identical between your reference annotations and genome sequence files).
  4. Optional: Bowtie indexes (must be small indexes (.ebwt)). Bowtie indexes will be created if the `ebwt` setting is empty in your Paths File.

*`tiny-count` accepts SAM files in your Samples Sheet only when invoked as an individual step. Because genome alignments are done after collapsing reads, the pipeline does not currently support SAM files from other sources.

### Running the End-to-End Analysis
In most cases you will use this toolset as an end-to-end pipeline. This will run a full, standard small RNA sequencing data analysis according to your configuration file. Before starting, you will need the following:

1. High-throughput sequencing data in fastq format. 
2. The genome sequence of interest in fasta format.
3. Genome coordinates of small RNA features of interest in gff3 format.
4. A completed Samples Sheet (`samples.csv`) with paths to the fastq files.
5. A completed Features Sheet (`features.csv`) with paths to the gff3 file(s).
6. An updated Paths File (`paths.yml`) with the path to the genome sequence.
7. A Run Config file (`run_config.yml`) located in your working directory or the path to the file. The template provided does not need to be updated if you wish to use the default settings.

To run an end-to-end analysis, be sure that you're working within the conda tinyrna environment (instructions above) in your terminal and optionally set your working directory that contains the Run Config file. Then, simply enter the following code into your terminal (if you are not working in the directory containing `run_config.yml`, provide the path before the name of the file - `path/to/run_config.yml`:

```
tiny run --config run_config.yml
```

### Running Individual Steps
The process for running individual steps differs depending on whether the step is a tinyRNA Python component, or a CWL wrapped third party tool.

The following steps are Python components that can be run from the command line within the tinyRNA conda environment:
##### Create Workflow
```
tiny-config -i CONFIG

  required arguments:
    -i CONFIG, --input-file CONFIG
                          The Run Config file to be processed
```
##### Collapser
```
tiny-collapse [-h] -i FASTQFILE -o OUTPREFIX [-t THRESHOLD] [-c]

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
tiny-count [-h] -i SAMPLES -c CONFIGFILE -o OUTPUTPREFIX
                  [-sf [SOURCE [SOURCE ...]]] [-tf [TYPE [TYPE ...]]] [-nn]
                  [-a] [-p] [-d]

This submodule assigns feature counts for SAM alignments using a Feature Sheet
ruleset. If you find that you are sourcing all of your input files from a prior
run, we recommend that you instead run `tiny recount` within that run's
directory.

optional arguments:
  -h, --help            show this help message and exit
  -sf [SOURCE [SOURCE ...]], --source-filter [SOURCE [SOURCE ...]]
                        Only produce counts for features whose GFF column 2
                        matches the source(s) listed
  -tf [TYPE [TYPE ...]], --type-filter [TYPE [TYPE ...]]
                        Only produce counts for features whose GFF column 3
                        matches the type(s) listed
  -nn, --no-normalize   Do not normalize counts by genomic hits and (selected)
                        overlapping feature counts.
  -a, --all-features    Represent all features in output counts table,
                        regardless of counts or identity rules.
  -p, --is-pipeline     Indicates that counter was invoked as part of a
                        pipeline run and that input files should be sourced as
                        such.
  -d, --report-diags    Produce diagnostic information about
                        uncounted/eliminated selection elements.

required arguments:
  -i SAMPLES, --input-csv SAMPLES
                        your Samples Sheet
  -c CONFIGFILE, --config CONFIGFILE
                        your Features Sheet
  -o OUTPUTPREFIX, --out-prefix OUTPUTPREFIX
                        output prefix to use for file names
```
##### Deseq2
```
tiny-deseq.r --input-file COUNTFILE --outfile-prefix PREFIX [--control CONDITION] [--pca] [--drop-zero]

    --input-file <count_file>
          A text file containing a table of features x samples of the run to
          process by DESeq2. The [...]feature_counts.csv output of tinyrna-count is expected here.
              
    --outfile-prefix <outfile>
          Name of the output files to write. These will be created:
              1. Normalized count table of all samples
              2. Differential gene expression table per comparison
              3. A PCA plot per comparison, if --pca is also provided.

    --control <control_condition>
          Opional. If the control condition is specified, comparisons will
          only be made between the control and experimental conditions.

    --pca
          Optional. This will produce principle component analysis plots
          using the DESeq2 library. Output files are PDF format.

    --drop-zero
          Optional. Prior to performing analysis, this will drop all
          rows/features which have a zero count in all samples."
```
##### Plotter
```
tiny-plot [-nc NORM_COUNTS] [-dge COMPARISON [COMPARISON ...]]
                 [-len 5P_LEN [5P_LEN ...]] [-h] [-o PREFIX] [-pv VALUE]
                 [-s MPLSTYLE] [-v] -p PLOT [PLOT ...]

This script produces basic static plots for publication as part of the tinyRNA
workflow. Input file requirements vary by plot type and you are free to supply
only the files necessary for your plot selections. If you are sourcing all of
your input files from the same run directory, you may find it easier to instead
run `tiny replot` within that run directory.

Required arguments:
  -p PLOT [PLOT ...], --plots PLOT [PLOT ...]
                        List of plots to create. Options:
                        • len_dist: A stacked barplot showing size & 5'
                          nucleotide distribution.
                        • class_charts: A pie and barchart showing proportions
                          of counts per class.
                        • replicate_scatter: A scatter plot comparing
                          replicates for all count files given.
                        • sample_avg_scatter_by_dge: A scatter plot comparing
                          all sample groups, with significantly different
                          genes highlighted. P-value can be set using
                          --p-value. Default: 0.05.
                        • sample_avg_scatter_by_dge_class: A scatter plot
                          comparing all sample groups, with classes and
                          significantly different genes highlighted.

Input files produced by Counter:
  -len 5P_LEN [5P_LEN ...], --len-dist 5P_LEN [5P_LEN ...]
                        The ...nt_len_dist.csv files produced by tiny-count

Input files produced by DGE:
  -nc NORM_COUNTS, --norm-counts NORM_COUNTS
                        The ...norm_counts.csv file produced by tiny-deseq.r
  -dge COMPARISON [COMPARISON ...], --dge-tables COMPARISON [COMPARISON ...]
                        The ...cond1...cond2...deseq.csv files produced by
                        tiny-deseq.r

Optional arguments:
  -h, --help            show this help message and exit
  -o PREFIX, --out-prefix PREFIX
                        Prefix to use for output filenames.
  -pv VALUE, --p-value VALUE
                        P-value to use in DGE scatter plots.
  -s MPLSTYLE, --style-sheet MPLSTYLE
                        Optional matplotlib style sheet to use for plots.
  -v, --vector-scatter  Produce scatter plots with vectorized points (slower).
                        Note: only the points on scatter plots will be raster
                        if this option is not provided.


```

##### fastp, bowtie-build, and bowtie
These are CWL wrapped third party tools.
1. Copy the workflow CWL folder to your current working directory with the command `tiny setup-cwl --config none`
2. Within `./CWL/tools` find the file for the step you wish to run. Navigate to this folder in terminal (or copy your target .cwl file to a more convenient location)
3. Run `cwltool --make-template step-file.cwl > step-config.YML`. This will produce a `YML` configuration file specific to this step. Optional arguments will be indicated as such; if you do not wish to set a value for an optional argument, best practice is to remove it from the file
4. Fill in your preferences and inputs in this step configuration file and save it
5. Execute the tool with the command `cwltool step-file.cwl step-config.YML`

### Using a Different Workflow Runner

We have used CWL to define the workflow for scalability and interoperability. The default runner, or interpreter, utilized by tinyRNA is `cwltool`. You may use a different CWL runner if you would like, and in order to do so you will need the workflow CWL and your **processed** Run Config file. The following will copy these files to your current working directory:

```
tiny setup-cwl --config <path/to/Run_Config.yml>
```

If you don't have a Run Config file or do not wish to obtain a processed copy, you may instead use "None" or "none" in the `--config` argument:
```
tiny setup-cwl --config none
```

## Outputs

The files produced by each pipeline step will be included in the final run directory by default. These intermediate files are organized into subdirectories by step.

### Data Pre-Processing

[fastp](https://github.com/OpenGene/fastp) is used to trim adapters and remove poor quality reads from FASTQ input files. Summary and quality statistics reports are generated for each library. These reports are used to calculate the pipeline summary statistics for total reads and retained reads.

### Collapsed FASTA files
A "collapsed" FASTA contains unique reads found in fastp's quality filtered FASTQ files. Each header indicates the number of times that sequence occurred in the input. This allows for faster bowtie alignments while preserving counts for downstream analysis.

### Counts and Pipeline Statistics
The counter step produces a variety of outputs

#### Feature Counts
Custom Python scripts and HTSeq are used to generate a single table of feature counts that includes columns for each library analyzed. A feature's _Feature ID_ and _Feature Class_ are simply the values of its `ID` and `Class` attributes. We have also included a _Feature Name_ column which displays aliases of your choice, as specified in the _Alias by..._ column of the Features Sheet. If _Alias by..._ is set to`ID`, the _Feature Name_ column is left empty.

For example, if your Features Sheet has a rule which specifies an ID Attribute of `sequence_name` and the GFF entry for this feature has the following attributes column:
```
... ID=406904;sequence_name=mir-1,hsa-miR-1;Class=miRNA; ...
```
The row for this feature in the feature counts table would read:

| Feature ID | Feature Name | Feature Class | Group1_rep_1 | Group1_rep_2 | ... |
|------------|--------------|---------------|-----------|-----------|-----|
| 406904 | mir-1,hsa-miR-1 | miRNA | 1234 | 999 | ... |


#### Alignment Statistics
A single table of alignment statistics includes columns for each library and the following rows:
- Total Assigned Reads (i.e. counts from sequences that aligned to at least one feature in your Features Sheet)
- Total Assigned Sequences (i.e. unique sequences that aligned to at least one feature in your Features Sheet)
- Assigned Single-Mapping Reads (i.e. counts from sequences mapping to a single genomic locus and aligning to at least one feature in your Features Sheet)
- Assigned Multi-Mapping Reads (i.e. counts from sequences mapping to multiple genomic loci that aligned to at leats one feature in your Features Sheet)
- Reads Assigned to Single Feature (i.e. counts from sequences that aligned to a single feature in your Features Sheet)
- Sequences Assigned to Single Feature (i.e. unique sequences that aligned to a single feature in your Features Sheet)
- Reads Assigned to Multiple Features (i.e. counts from sequences that aligned to multiple features in your Features Sheet)
- Sequences Assigned to Multiple Features (i.e. unique sequences that aligned to multiple features in your Features Sheet)
- Total Unassigned Reads (i.e. total counts for sequences that didn't align to any features in your Features Sheet)
- Total Unassigned Sequences (i.e. total unique sequences that didn't align to any features in your Features Sheet)

#### Summary Statistics
A single table of summary statistics includes columns for each library and the following rows:
- Total Reads (i.e. total reads represented in FASTQ input files)
- Retained Reads (i.e. total reads passing quality filtering)
- Unique Sequences (i.e. total unique sequences passing quality filtering)
- Mapped Sequences (i.e. total genome-mapping sequences passing quality filtering)
- Assigned Reads (i.e. total genome-mapping reads passing quality filtering that were aligned to at least one feature in your Features Sheet)

#### 5'nt vs. Length Matrix

After alignment, a size and 5' nt distribution table is created for each library. The distribution of lengths and 5' nt can be used to assess the overall quality of your libraries. This can also be used for analyzing small RNA distributions in non-model organisms without annotations.

### Differential Expression Analysis
DGE is performed using the `DESeq2` R package. It reports differential expression tables for your experiment design, and a table of normalized feature counts. If your control condition is indicated in your Samples Sheet then pairwise comparisons will be  made against the control. If a control condition is not indicated then all possible bidirectional pairwise comparisons are made. 

### Plots
Simple static plots are generated from the outputs of Counter and Deseq2. These plots are useful for assessing the quality of your experiment design and the quality of your libraries. The available plots are:
- len_dist: A stacked barplot showing size & 5' nucleotide distribution.
- class_charts: A pie and barchart showing proportions of counts per class.
- replicate_scatter: A scatter plot comparing replicates for all count files given.
- sample_avg_scatter_by_dge: A scatter plot comparing all sample groups, with significantly different genes highlighted.
- sample_avg_scatter_by_dge_class: A scatter plot comparing all sample groups, with classes and significantly different genes highlighted.


## Contributing

See the [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines. To see what is actively being worked or planned go to the [projects tab](https://github.com/MontgomeryLab/tinyrna/projects) or the [issues tab](https://github.com/MontgomeryLab/tinyrna/issues).

## Authors

* **Kristen Brown** - 2018-2019 - Colorado State University - [biokcb](https://github.com/biokcb)
* **Alex Tate** - 01/2021-present - Colorado State University - [AlexTate](https://github.com/AlexTate)

See also the list of [contributors](https://github.com/MontgomeryLab/tinyrna/contributors) who participated in this project.

## License

This project is licensed under the GPLv3 license (along with HTSeq, bowtie). License - see the [LICENSE.md](LICENSE.md) file for details
