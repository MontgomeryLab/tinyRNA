## Pipeline Commands
The following commands deal with pipeline operations for carrying out end-to-end analyses:

```shell
# End-to-end analysis
tiny run --config run_config.yml

# Resume prior analyses
tiny recount --config processed_run_config.yml
tiny replot --config processed_run_config.yml

# Retrieving config files
tiny get-template
tiny setup-cwl
```

### Running End-to-End
The `tiny run` command performs a comprehensive analysis of your [input files](../README.md#requirements-for-user-provided-input-files) according to the preferences defined in your [configuration files](Configuration.md).

### Resuming a Prior Analysis
The Counter and Plotter steps offer a wide variety of options for refining your analysis. You might find that repeat analyses are required while tuning these options to your goals. However, the earlier pipeline steps (fastp, Collapser, and bowtie) handle the largest volume of data and are resource intensive, so you can save time by reusing their outputs for subsequent analyses. One could do so by running the later steps individually (e.g. `tiny-count`, `tiny-deseq.r`, and `tiny-plot`), but assembling their commandline inputs by hand is labor-intensive and prone to spelling mistakes.

The commands `tiny recount` and `tiny replot` seek to solve this problem. As discussed in the [Run Config documentation](Configuration.md#the-processed-run-config), the Run Directory for each end-to-end analysis will contain a processed Run Config, and this is the file that determines the behavior of a resume run. 

You can modify the behavior of a resume run by changing settings in:
- The **processed** Run Config
- The **original** Features Sheet that was used for the end-to-end run (as indicated by the `features_csv` key in the **processed** Run Config)

#### The Steps
1. Make and save the desired changes in the files above
2. In your terminal, `cd` to the Run Directory of the end-to-end run you wish to resume
3. Run the desired resume command

#### A Note on File Inputs
File inputs are sourced from the **original** output subdirectories of prior steps in the target Run Directory. For `tiny replot`, this means that files from previous executions of `tiny recount` will **not** be used as inputs; only the original end-to-end outputs are used.

#### Where to Find Outputs from Resume Runs
Output subdirectories for resume runs can be found alongside the originals, and will have a timestamp appended to their name to differentiate them.
