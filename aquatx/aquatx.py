#!/usr/bin/env python
"""
The main entry point for AQuATx for small RNA data analysis. This tool provides an end-to-end
workflow for analyzing small RNA sequencing data from raw fastq files. This entry point also provides
options for only returning workflows that can be used separately or to retrieve template files.

When installed this script should be run as:

aquatx <subcommand> <args>
"""
import os
import sys
import subprocess
from shutil import copyfile
from pkg_resources import resource_filename

def main():
    """
    The main routine that determines what type of run to do.

    Options:

        run: Run the end-to-end analysis based on a config file.
        get-template: Get the input sheets & template config files.
        setup-cwl: Get the CWL workflow for a run
    """

    # Get the package data
    aquatx_cwl_path = resource_filename('aquatx', 'cwl/')
    aquatx_extras_path = resource_filename('aquatx', 'extras/')

    if sys.argv[1] == "run":
        print("Running the end-to-end analysis...")

        # First get the configuration file set up
        conf = subprocess.run(['aquatx-config'] + sys.argv[2:], stdout=subprocess.PIPE)
        config_file = conf.stdout.decode('utf-8')

        # Run with cwltool
        subprocess.run('cwltool ' + aquatx_cwl_path + 'workflows/aquatx_wf.cwl '
                       + config_file,
                       shell=True, stdout=subprocess.PIPE)

    elif sys.argv[1] == "get-template":
        print("Copying template input files to current directory...")

        copyfile(aquatx_extras_path + 'run_config_template.yml', './run_config_template.yml')
        copyfile(aquatx_extras_path + 'sample_sheet_template.csv', './sample_sheet_template.csv')
        copyfile(aquatx_extras_path + 'reference_sheet_template.csv', './reference_sheet_template.csv')

    elif sys.argv[1] == "setup-cwl":
        print("Creating cwl workflow...")

        if sys.argv[3] not in ('None', 'none'):
            print("The configuration file: ")
            subprocess.run(['aquatx-config'] + sys.argv[2:], stdout=subprocess.PIPE)

        print("\nThe workflow and files are under: cwl/tools/ and cwl/workflows/")
        os.makedirs('cwl/tools/', exist_ok=True)
        os.makedirs('cwl/workflows', exist_ok=True)
        copyfile(aquatx_cwl_path + 'tools/bowtie.cwl', 'cwl/tools/bowtie.cwl')
        copyfile(aquatx_cwl_path + 'tools/aquatx-collapse.cwl', 'cwl/tools/aquatx-collapse.cwl')
        copyfile(aquatx_cwl_path + 'tools/aquatx-count.cwl', 'cwl/tools/aquatx-count.cwl')
        copyfile(aquatx_cwl_path + 'tools/aquatx-deseq.cwl', 'cwl/tools/aquatx-deseq.cwl')
        copyfile(aquatx_cwl_path + 'tools/aquatx-merge.cwl', 'cwl/tools/aquatx-merge.cwl')
        copyfile(aquatx_cwl_path + 'tools/bowtie-build.cwl', 'cwl/tools/bowtie-build.cwl')
        copyfile(aquatx_cwl_path + 'tools/fastp.cwl', 'cwl/tools/fastp.cwl')
        copyfile(aquatx_cwl_path + 'workflows/aquatx_wf.cwl', 'cwl/workflows/aquatx_wf.cwl')

    elif sys.argv[1] == "setup-nextflow":
        print("Creating nextflow workflow...")
        print("This command is currently not implemented.")

    else:
        print("Command not recognized...")

if __name__ == '__main__':
    main()
