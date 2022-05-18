# Getting Started

This folder (`START_HERE`) contains a working minimal configuration and a sample input dataset derived from C. elegans chromosome 1 of Wormbase WS279. We've assembled this configuration to make it easy to start using tinyRNA, and to provide a basis for your own project configuration.

Complete reference data for C. elegans WS279 can be found here: https://www.montgomerylab.org/resources.html

## Installation

See the [README](../README.md#installation) for installation instructions and tips.

## This folder

Here's what you'll find:
- **fastq_files**: contains library fastq files (10k line arbitrary subset)
- **reference_data**: contains reference annotation and genome files for chromosome 1
- **features.csv**: spreadsheet of selection rules for counting features
- **paths.yml**: configuration file for defining pipeline file inputs
- **run_config.yml**: configuration file for defining pipeline preferences for each step and for the overall pipeline run
- **samples.csv**: spreadsheet for defining library files and their associated group and replicate numbers

## First run
The configuration is tied together with `run_config.yml`, so this is what you will pass to the pipeline. Since we already have a working configuration let's run an end-to-end analysis on our sample data using the command:
```
cd START_HERE
tiny run --config run_config.yml
```
Did you receive "command not found"? Make sure that you activate the tinyrna environment before using it.
```
conda activate tinyrna
```
And when you're done, you can close your terminal or use `conda deactivate` to return to a normal shell.

### Terminal output
The output you see on your terminal is from `cwltool`, which coordinates the execution of the workflow CWL. The terminal output from individual steps is redirected to a logfile for later reference.

### File outputs
When the analysis is complete you'll notice a new folder has appeared whose name contains the date and time of the run. Inside you'll find subdirectories containing the file and terminal outputs for each step, and the processed Run Config file for auto-documentation of the run.

### Bowtie indexes
Bowtie indexes were built during this run because paths.yml didn't define an `ebwt` prefix. Now, you'll see the `ebwt` points to the freshly built indexes in your run directory. This means that indexes won't be rebuilt during any subsequent runs that use this `paths.yml` file. If you need to rebuild your indexes:
1. Change the value of ebwt to `ebwt: ''` in paths.yml
2. Ensure that your Run Config file contains `run_bowtie_build: True`

## Running Your Data
Expected runtime: ~15-30 minutes
1. Edit your GFF or GTF file so that it meets the requirements outlined in [the README](../README.md#requirements-for-user-provided-input-files)
2. Move your GFF and genome sequence files into the reference_data directory.
3. Edit features.csv and samples.csv file for your datasets and selection parameters.
4. Edit paths.yml as follows:
   - line 36: `ebwt: ''` (no value)
   - line 40: `reference_genome_files: your-fasta-formatted-DNA-sequence-file`
5. Run the pipeline with the command: `tiny run --config run_config.yml`