#!/usr/bin/env python
"""
The main entry point for AQuATx for small RNA data analysis. This tool provides an end-to-end
workflow for analyzing small RNA sequencing data from raw fastq files. This entry point also provides
options for only returning workflows that can be used separately or to retrieve template files.

When installed this script should be run as:

aquatx <subcommand> --config <config-file>
"""
import aquatx.srna.configuration_setup
import subprocess
import shutil
import os

from pkg_resources import resource_filename
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser(
        description=
        "This tool provides an end-to-end workflow for analyzing small RNA sequencing data" +
        "from raw fastq files. This entry point also provides options for only returning" +
        "workflows that can be used separately or to retrieve template files."
    )

    # Parser for subcommands: (run, setup-cwl, get-template, etc.)
    subparsers = parser.add_subparsers(required=True, dest='command')

    # Subcommand get-template has no additional arguments
    subparsers.add_parser("get-template")

    # Subcommands that require a configuration file argument
    commands_with_configfile = ["run", "setup-cwl", "setup-nextflow"]
    for command in commands_with_configfile:
        subparsers.add_parser(command).add_argument(
            '--config', metavar='configFile', required=True,
        )

    return parser.parse_args()


def run(aquatx_cwl_path, config_file):
    print("Running the end-to-end analysis...")

    # First get the configuration file set up for this run
    workflow_conf_file = aquatx.srna.configuration_setup.setup_config(config_file)

    # Run with cwltool
    subprocess.run(f"cwltool {aquatx_cwl_path}workflows/aquatx_wf.cwl {workflow_conf_file}", shell=True)


def get_template(aquatx_extras_path):
    print("Copying template input files to current directory...")
    template_files = [
        'run_config_template.yml',
        'sample_sheet_template.csv',
        'reference_sheet_template.csv'
    ]

    # Copy template files to the current working directory
    for template in template_files:
        shutil.copyfile(f"{aquatx_extras_path}/{template}", f"{os.getcwd()}/{template}")


def setup_cwl(aquatx_cwl_path, config_file):
    print("Creating cwl workflow...")

    # If the word "None" or "none" is supplied, simply copy workflow files. No config file processing.
    if config_file not in ('None', 'none'):
        # Set up the config file
        processed_config_location = aquatx.srna.configuration_setup.setup_config(config_file)
        print("The processed configuration file is located at: " + processed_config_location)

    # Copy the entire cwl directory to the current working directory
    shutil.copytree(aquatx_cwl_path, os.getcwd() + "/cwl/")
    print("The workflow and files are under: cwl/tools/ and cwl/workflows/")


def setup_nextflow(config_file):
    print("Creating nextflow workflow...")
    print("This command is currently not implemented.")


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
    command_map = {
        "run": lambda: run(aquatx_cwl_path, args.config),
        "setup-cwl": lambda: setup_cwl(aquatx_cwl_path, args.config),
        "get-template": lambda: get_template(aquatx_extras_path),
        "setup-nextflow": lambda: setup_nextflow(args.config)
    }

    command_map[args.command]()


if __name__ == '__main__':
    main()
