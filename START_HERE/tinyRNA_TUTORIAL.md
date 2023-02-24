# Getting Started With the tinyRNA Workflow

This folder (`START_HERE`) contains a working minimal configuration and a generated/simulated sample dataset. We've assembled this configuration to make it easy to start using tinyRNA, and to provide a basis for your own project configuration.

## Installation

See the [README](../README.md#tinyrna-installation) for installation instructions and tips.

## This folder

Here's what you'll find:
- **fastq_files**: contains generated sample FASTQ files
- **reference_data**: contains a reference genome file with random DNA sequences, and a reference annotation file with simulated features selected from the genome
- **features.csv**: spreadsheet of selection rules for counting features
- **paths.yml**: configuration file for defining the pipeline's main file inputs
- **run_config.yml**: configuration file for defining preferences for each pipeline step and the overall pipeline run
- **samples.csv**: spreadsheet for defining the group name, replicate number, etc. for each input FASTQ file

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
Bowtie indexes were built during this run because `paths.yml` didn't define an `ebwt` prefix. Now, you'll see the `ebwt` points to the freshly built indexes in your run directory. This means that indexes won't be rebuilt during any subsequent runs that use this `paths.yml` file. If you need to rebuild your indexes, simply delete the value to the right of `ebwt` in paths.yml

## Running Your Data
Expected runtime: ~10-60 minutes (expect longer runtimes if a bowtie index must be built)
1. Edit your GFF or GTF file so that it meets the requirements outlined in [the README](../README.md#requirements-for-user-provided-input-files)
2. Move your GFF and genome sequence files into the reference_data directory.
3. Edit features.csv and samples.csv file for your datasets and selection parameters.
4. Edit paths.yml as follows:
   - line 20: change the value after `path:` to point to your GFF or GTF file
   - line 46: delete the value after `ebwt:`
   - line 51: change the value after `- ` to point to your fasta formatted DNA sequence file
5. Run the pipeline with the command: `tiny run --config run_config.yml`