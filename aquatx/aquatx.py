#!/usr/bin/env python
"""
The main entry point for AQuATx for small RNA data analysis. This tool provides an end-to-end
workflow for analyzing small RNA sequencing data from raw fastq files. This entry point also provides
options for only returning workflows that can be used separately or to retrieve template files.

When installed this script should be run as:

aquatx <subcommand> --config <config-file>
"""
import srna.configuration_setup
import subprocess
import shutil

from pkg_resources import resource_filename
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser(
        description=
        "This tool provides an end-to-end workflow for analyzing small RNA sequencing data" +
        "from raw fastq files. This entry point also provides options for only returning" +
        "workflows that can be used separately or to retrieve template files."
    )

    # Parser for subcommands run, setup-cwl, setup-nextflow, get-template
    subparsers = parser.add_subparsers(required=True, dest='command')

    # get-template has no additional arguments
    subparsers.add_parser("get-template")

    # Commands that require a configuration file parameter
    commands_with_configfile = ["run", "setup-cwl", "setup-nextflow"]
    for command in commands_with_configfile:
        subparsers.add_parser(command).add_argument(
            '--config', metavar='configFile', required=True,
        )

    return parser.parse_args()


def run(aquatx_cwl_path, config_file):
    print("Running the end-to-end analysis...")

    # First get the configuration file set up for this run
    workflow_conf_file = srna.configuration_setup.setup_config(config_file)

    # Run with cwltool
    subprocess.run(f"cwltool {aquatx_cwl_path}workflows/aquatx_wf.cwl {workflow_conf_file}", shell=True)


def get_template(aquatx_extras_path):
    print("Copying template input files to current directory...")
    template_files = [
        'run_config_template.yml',
        'sample_sheet_template.csv',
        'reference_sheet_template.csv'
    ]

    for template in template_files:
        shutil.copyfile(f"{aquatx_extras_path}{template}", '.')


def setup_cwl(aquatx_cwl_path, config_file):
    print("Creating cwl workflow...")

    if config_file not in ('None', 'none'):
        # Set up the config file
        processed_config_location = srna.configuration_setup.setup_config(config_file)
        print("The processed configuration file is located at: " + processed_config_location)

    print("The workflow and files are under: cwl/tools/ and cwl/workflows/")
    shutil.copytree(aquatx_cwl_path, '.')


def main():
    """
    The main routine that determines what type of run to do.

    Options:
        run: Run the end-to-end analysis based on a config file.
        get-template: Get the input sheets & template config files.
        setup-cwl: Get the CWL workflow for a run
    """

    # Parse command line arguments
    args = get_args()

    # Get the package data
    aquatx_cwl_path = resource_filename('aquatx', 'cwl/')
    aquatx_extras_path = resource_filename('aquatx', 'extras/')

    # Execute appropriate command based on command line input
    if args.command == "run":
        run(aquatx_cwl_path, args.config)
    elif args.command == "get-template":
        get_template(aquatx_extras_path)
    elif args.command == "setup-cwl":
        setup_cwl(aquatx_cwl_path, args.config)
    elif args.command == "setup-nextflow":
        print("Creating nextflow workflow...")
        print("This command is currently not implemented.")


if __name__ == '__main__':
    main()
