# Getting Started

This folder (`START_HERE`) contains a working configuration and a sample input dataset derived from C. elegans chromosome 1 of Wormbase WS279. We've assembled this configuration to make it easy to start using AQuATx, and to provide a basis from which you may model your project's configuration.

## Installation

Installation is pretty easy. You'll install AQuATx in an isolated environment so that it (and its dependencies) don't interfere with your daily. We've provided an environment file to make things simple. First, you'll need to install `conda` for managing the environment.

#### 1. Install conda
To install `conda` with the least time and disk space commitment, you can download and follow instructions for [Miniconda3](https://docs.conda.io/en/latest/miniconda.html), which contains only conda and its dependencies.

#### 2. Install R & DESeq2
Rather than installing R via conda, we recommend you install it yourself first from [CRAN](https://www.r-project.org/), then install DESeq2 following [instructions](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) in Bioconductor. 

#### 3a. XQuartz (mac only)
DESeq2 is quiet about it, but if you're on a mac you'll also need to download and install [XQuartz](https://www.xquartz.org/).

#### 3b. Install AQuATx

To install AQuATx and its remaining dependencies:
```
# Clone the repository into a local directory
git clone https://github.com/MontgomeryLab/aquatx-srna.git
cd aquatx-srna

# Install the aquatx-srna environment and dependencies
conda env create -f environment.yml
```

## This folder

Here's what you'll find:
- **bowtie_indexes**: contains indexes for chromosome 1
- **reference_data**: contains reference annotation and genome files for chromosome 1
- **sample_data**: contains library fastq files (only 10k lines each for fast runs)
- **features.csv**: spreadsheet of selection rules for counting features
- **paths.yml**: configuration file for defining pipeline file inputs
- **run_config.yml**: configuration file for defining pipeline preferences for each step and for the overall pipeline run
- **samples.csv**: spreadsheet for defining library files and their associated group and replicate numbers

## First run
All of the above files are mapped out in `run_config.yml` and the tree of configuration files that it references. Since we already have a working configuration let's run an end-to-end analysis on our sample data using the command:
```
aquatx run --config run_config.yml
```
Did you receive "command not found"? Make sure that you activate the aquatx environment before using it.
```
conda activate aquatx-srna
```
And when you're done, you can close your terminal or use `conda deactivate` to return to a normal shell.

### Terminal output
The output you see on your terminal is first from `cwltool`, which coordinates the execution of the workflow CWL, as well as the individual steps being executed.

### File outputs
When the analysis is complete you'll notice a new folder has appeared whose name contains the date and time of the run. This is for record keeping and to prevent outputs from being overwritten between runs. Inside you'll find folders which contain the file outputs of each pipeline step, feature counts, pipeline and counting statistics, and the processed Run Config file which determined the behavior of this run. 

You'll notice that there isn't a `bowtie-build` folder. The Run Config contains `run_bowtie_build: True`, so why was this step skipped? Let's investigate with a second run.

## Second Run
Let's pretend we don't already have prebuilt bowtie indexes and we'd like the pipeline to build them for us. We will need to change two configuration files in order to do so.
1. Open `run_config.yml` in your preferred text editor. Verify that `run_bowtie_build` is True.
2. Open `paths.yml`. The current value of `ebwt` should be a prefix for the indexes in the folder `bowtie_indexes`. Change this value to '' (that's two single-quotes, but two double quotes is fine too)
3. We already have our reference genome file listed in this configuration, but you'll want to make sure this is set and correct when you begin working with your own datasets.
4. Save and close `paths.yml`

Rerun the pipeline with the same command we used for the first run. You'll notice that now there is a `bowtie-build` folder in the most recent output folder. If you open `paths.yml` again, you'll notice that the value of `ebwt` is now a prefix for the indexes we just built. This means that subsequent runs will use these new indexes rather than rebuilding them each time. If you'd like, run the pipeline a third time to verify that indexes were not rebuilt and that bowtie ran successfully using the second run's indexes.

# Configuration Files
![Configuration File Tree](../images/config-files.png)
### Run Config
When setting up for running the pipeline with your own data, you'll want to start in your Run Config file. Make sure that the value for `paths_config` points to your Paths file. All other settings in this file are optional. You may change settings for each pipeline step here. Not all available options for each step are listed; we've included only the settings relevant to sRNA research, and default values reflecting our own "best practices."

>**Tip:**
> You may use relative or absolute paths in these configuration files. If you use relative paths, they will be evaluated relative to the file that contains them. This allows you to store these files in separate locations and more flexibly organize your project.

### Paths
Next, open the Paths file. This is where you will provide inputs for `bowtie-build` if you do not already have bowtie indexes for your project. Also verify that `samples_csv` and `features_csv` point to their corresponding files.

### Features and Samples Sheets
Finally, the Samples Sheet and Features Sheet will need to be filled out with information about your library files and the selection rules for the features you wish to count. These files should be edited with a spreadsheet editor such as Microsoft Excel or LibreOffice Calc.

# Feature Selection
sRNA research is often confounded by the fact that features of interest often overlap with superabundant features which may not be relevant. In these situations it is unclear which feature a sequence originated from, and when tallying feature counts this can lead to the inflation or deflation of the counts of features of interest. We provide a **Features Sheet** for defining rules to dictate which kinds of features should receive counts in these cases of overlap.

Rules can be broken down into two parts: selection parameters for feature attributes, and selection parameters for read attributes.

>**Important**: candidate features do not receive counts if they do not pass selection process described below

## Feature Attribute Parameters
| Attribute Key | Attribute Value | Hierarchy |
| --- | --- | --- |

The first round of selection is performed using this portion of the rule. Each rule must be assigned a hierarchy value, which is applies only to this round of selection. A lower hierarchy value indicates higher selection preference and rules may share hierarchy values.

Each feature associated with a read-locus (alignment) pair will be examined to see if its attributes contain the key-value pairs defined in your rules list. If multiple features match, then elimination will be performed using the hierarchy values of the rules they matched. The feature-rule pair(s) with the lowest hierarchy value will be selected for a second round of elimination. Internally, these matches are represented as `(hierarchy, rule, feature)` tuples which we call _hits_.

A feature may match multiple rules. When this happens, a _hit_ is produced for each matched rule and normal hierarchical elimination is performed. This means that the product of the first round of elimination is not just a list of features, but rather feature-rule pairs in the form of _hits_.

>**Tip**: These parameters are case sensitive. The capitalization in your rule must match the capitalization in your GFF files

## Read Attribute Parameters
| Strand | 5' End Nucleotide | Length | Match |
| --- | --- | --- | --- |

The second round of selection switches to attributes of the read to which features are being assigned. Recall that the first round of selection yields a list of _hits_. Now, only the rules contained in this list of _hits_ are used for selection. Contrast this with the first round in which the key-value pairs of all rules were considered.

### Strand
Valid values are `sense`, `antisense`, and `both`.

### Match
Valid values are `Partial` and `Full`. This parameter is referring to the required amount of overlap between a read alignment and a feature in order for that feature to be a candidate for selection. If a rule specifies a `Partial` Match value, then features which overlap a read alignment by at least one base will be considered for selection. If a rule specifies a `Full` Match value, then only features whose enpoints are fully contained by or equal to the read alignment will be considered for selection.

### 5' End Nucleotide and Length
| Parameter | Single | List | Range | Wildcard |
| --- | :---: | :---: | :---: | :---: |
| 5' end nt | X | X |  | X |
| Length | X | X | X | X |

Examples:
- **Single**: `G` or `22`
- **List**: `C,G,U` or `25, 26` (spaces do not matter)
- **Range**: `20-25`
- **Wildcard**: `all`
- **Mixed**: `19, 21-23, 25-30`

>**Tip:** these parameters are **not** case sensitive.

>**Tip:** you may specify U and T bases in your rules. Uracil bases will be converted to thymine when your Features Sheet is loaded.

### Misc
| Name Attribute | Feature Source |
| --- | --- |

You may specify a **Name Attribute** which will be used for the `Feature Name` column of the Feature Counts output table. The intention of this column is to provide a human-friendly name for each feature. For example, if one of your rules specifies a **Name Attribute** of `sequence_name` and gene1's `sequence_name` attribute is "abc123", then gene1's `Feature Name` column in the Feature Counts table will read "abc123".

The **Feature Source** field of a rule is tied only to the **Name Attribute**; rules are _not_ partitioned on a GFF file basis, and features parsed from these GFF files are similarly not partitioned as they all go into the same lookup table regardless of source. However, the **Name Attribute** is used to build the alias table for `Feature Name` when parsing its corresponding **Feature Source**. Additionally, each GFF file is parsed only once regardless of the number of times it occurs in the Features Sheet.

### The Nitty Gritty Details
You may encounter the following cases when you have more than one unique GFF file listed in your **Feature Source**s:
- Attribute value lists are supported. Values which contain commas are treated as lists of values when parsing GFF files
- If a feature is defined in one GFF file, then again in another but with differing attributes, then those attribute values will be appended
- A feature may have multiple **Name Attribute**s associated with it
- If a feature is defined in one GFF file, then again but under a different **Name Attribute**, then both aliases are retained and treated as a list. All aliases will be present in the `Feature Name` column of the Feature Counts output table. They will be comma separated.