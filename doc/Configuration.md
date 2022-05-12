- [Overview](#run-config)
    - [Paths File](#paths-file)
    - [Samples Sheet](#samples-sheet)
    - [Features Sheet](#features-sheet)
    - [Plot Stylesheet](#plot-sheet)

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

>**Tip**: Each of the following files will allow you to map out the paths to all of the input files required for analysis. You can use either relative or absolute paths to do so. **Relative paths will be evaluated relative to the file in which they are defined.** This gives you the freedom to organize your project as you see fit.

#### Run Config

The overall behavior of the pipeline and its steps is determined by the Run Config file (`run_config.yml`). This YAML file can be edited using a simple text editor. Within it you must specify the location of your Paths file (`paths.yml`). All other settings are optional. When the pipeline starts up, tinyRNA will process the Run Config based on the contents of it and your other configuration files, and the processed copy will be saved to your run directory. The processed configuration is what ultimately determines the behavior of the workflow. This provides auto-documentation for all of your pipeline runs.

#### Paths File

The locations of pipeline file inputs are defined in the Paths file (`paths.yml`). This YAML file includes paths to your Samples and Features Sheets, in addition to your bowtie index prefix (optional) and the final run directory name. The final run directory will contain all pipeline outputs. The directory name is prepended with the `run_name` and current date and time to keep outputs separate.

#### Samples Sheet

To make it simple to specify your fastq files and their locations, along with associated sample names and replicate numbers, we provide a CSV file (`samples.csv`) which can be edited in a spreadsheet editor such as Microsoft Excel or LibreOffice Calc.

#### Features Sheet

Small RNAs can often be classified by sequence characteristics, such as length, strandedness, and 5' nucleotide. We provide a Features Sheet (`features.csv`) in which you can define selection rules to more accurately capture counts for the small RNAs of interest.  Rules apply to features parsed from **all** Feature Sources, with the exception of "Alias by..." which only applies to the Feature Source on the same row. These rules are used to select among overlapping features at each alignment locus. Selection first takes place against feature attributes (GFF column 9), and is directed by defining the attribute you want to be considered (Select by...) and the acceptable values for that attribute (with value...). Rules that match features at this stage will be used in a second stage of selection which includes elimination by hierarchy and interval overlap characteristics. Remaining candidates pass to the third and final stage of selection which examines characteristics of the alignment itself: strand relative to the feature of interest, 5' end nucleotide, and length. See [Counter](doc/Counter.md) for more information.

>**Tip**: Don't worry about having duplicate Feature Source entries. Each GFF file is parsed only once.

#### Plot Stylesheet

Matplotlib styles can be optionally overridden using a matplotlib stylesheet. The stylesheet provided by `tiny get-template` contains the defaults utilized by Plotter. Please keep in mind that Plotter overrides these defaults for a few specific elements of certain plots. Feel free to reach out if there is a plot style you wish to override but you find that you are unable to.