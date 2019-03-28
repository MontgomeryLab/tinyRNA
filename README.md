# AQuATx for small RNA sequencing data
:warning: **Under Development & Testing** :warning:

*This repository is being actively developed and tested, thus incomplete and only recommended for testing and playing around with until a release is made. Feedback, suggestions, and bug reports are welcome under the [issues tab](https://github.com/biokcb/aquatx-srna/issues). Thank you!*


AQuATx (Automated QUantitative Analysis of Transcript eXpression) is a set of tools to simplify the analysis of next-generation sequencing data. The goal of this specific repository is to provide an entire workflow for processing small RNA sequencing data. 

The current workflow is as follows:

![AQuATx basic pipeline](images/aquatx-workflow_current.png)

All tools are wrapped with CWL. The optional steps are teal, purple are part of the end-to-end standard workflow, and steps outlined in orange are Python scripts specifically created for this tool. 

## Getting Started

### Prerequisites & Installation

The basic AQuATx pipeline for small RNA analysis depends on Python (tested for 3.6), numpy, pandas, matplotlib, HTSeq, fastp, bowtie, R, and DESeq2. To install the pipeline it is also easiest to install with `conda`, otherwise you'll likely need to install each tool separately. 

#### 1. Install conda
To install `conda` you can either download and follow instructions for:

[Anaconda](https://www.anaconda.com/distribution/), which comes packaged with a few dependencies already or 
[miniconda3](https://docs.conda.io/en/latest/miniconda.html), which doesn't install unnecessary packages. 

#### 2. Install R & DESeq2
Since there are some known conflicts with installing R via conda, I recommend you install it yourself first from [CRAN](https://www.r-project.org/) and DESeq2 following [instructions](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) in Bioconductor. 

#### 3. Install AQuATx

To install AQuATx & the remaining dependencies:
```
# Clone the repository into a local directory
git clone https://github.com/MontgomeryLab/aquatx-srna.git
cd aquatx-srna

# Install dependencies using conda
conda install -c bioconda --file requirements.txt

# Run the aquatx setup script
python setup.py install
```

## Usage

There are several ways to use the pipeline, the easiest of which is to run an end-to-end analysis with `aquatx`. You can also access the individual scripts for collapsing reads, counting, differential expression, and plotting separately. More detailed documentation to be added to the Wiki soon.

The inputs to the pipeline are a sample sheet, a reference sheet, and a configuration file. To obtain the templates for these inputs use the following command:

```
aquatx get-template
```

### Configuration file:

The pipeline revolves around a configuration file to make it easy to set up and run. The configuration file is a `YAML` that can be edited in a simple text editor. This template is then used alongside the sample and references sheet to produce your workflow and input configuration file for the run. The updated configuration file is saved in the final run directory, under `config_files` with the date and time of the run. You can modify either just a few lines in the input or modify many of the options per tool. 

:warning: Not all items in the configuration are implemented just yet. bowtie indexing, turning on/off steps/intermediate file saving, and sorting into output directories are unavailable, but are next to be developed. See the [projects tab](https://github.com/MontgomeryLab/aquatx-srna/projects) for more information on what is being worked on.

### Sample sheet:

To make it simple to specify your sample files, the associated sample names and replicate numbers, we provide a `.csv` file that can be edited with a program like Excel. Simply be sure to save as a `.csv`, not to add commas in your cells, and to fill out the information appropriately. 

### Reference sheet:

You might have multiple reference annotation files and annotations you want to count differently than others based on what class of small RNA you are looking at. We allow you to specify counting both sense and antisense reads, masking locations from being counted toward your reference, and specifying subclass information. Class information should be specified in the 3rd column of your gff3 file.

### Input formats:

The accepted inputs to the analysis workflow are:
  1. Reference annotations must be in GFF3 format
  2. FastQ/A files
  3. bowtie index

### Running the end-to-end analysis

The main option for using the workflow is to run an end-to-end analysis. This will automatically take your input configuration and run a full, standard small RNA sequencing data analysis pipeline from raw fastq files to DEGs and plots. This utilizes `cwltool` to run a CWL based workflow. This will first run the configuration setup based on your input files.

```
aquatx run --config <path/to/config.yml>
```

### Creating a configuration to run separately

If you would instead like to create a configuration file and a CWL workflow to run the analysis using a different implementation (such as `CWLEXEC`, on DNANexus, using `Toil`, etc), you can just create the workflow and input files using:

```
aquatx setup-cwl --config <path/to/config.yml>
```

## Outputs

The pipeline will produce all intermediate files by default. We provide summary statitics of the run itself at each step, a differential gene expression table, and multiple visualizations in a vector-editable format. 

### Fastq file quality analysis

`fastp` produces summary and quality statistics for each of your raw fastq files and we collect that information to output into summary HTML reports. 

### Size distributions

After alignment, a size and 5'nt distribution table and histogram is created. The distribution of lengths and 5'nt can be used to assess the overall quality of your libraries. This can also be used for analyzing small RNA distributions in non-model organisms without annotations.

### Class analysis

One of the main parts of small RNA analysis is grouping small RNAs by class. `aquatx` produces pie and bar charts to display the class distribution of your aligned reads. 

### Differential gene expression

For runs with replicates we also perform differential gene expression analysis. The output produced includes a table of all small rnas, their fold change among comparisons, and associated p-values. 

## Contributing

See the [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines. To see what is actively being worked or planned go to the [projects tab](https://github.com/MontgomeryLab/aquatx-srna/projects) or the [issues tab](https://github.com/MontgomeryLab/aquatx-srna/issues).

## Authors

* **Kristen Brown** - Colorado State University - [biokcb](https://github.com/biokcb)

See also the list of [contributors](https://github.com/MontgomeryLab/aquatx-srna/contributors) who participated in this project.

## License

This project is licensed under the GPLv3 license (along with HTSeq, bowtie). License - see the [LICENSE.md](LICENSE.md) file for details
