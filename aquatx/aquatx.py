#!/usr/bin/env python
"""The main entry point for AQuATx for small RNA data analysis.

This tool provides an end-to-end workflow for analyzing small RNA sequencing
data from raw fastq files. This entry point also provides options for only
returning template files and workflows that can be used separately.

Subcommands:
    - get-template
    - setup-cwl
    - run

When installed, run and setup-cwl should be invoked with:
    aquatx <subcommand> --config <config-file>

A configuration file should be supplied for the run subcommand (required)
and for the setup-cwl subcommand (optional; alternatively you may use the
word "None" or "none" to obtain only the workflow files). This config file
will be processed and rewritten to reflect the workflow inputs defined under
the keys ` samples_csv` and `reference_sheet_file` (use get-template for
more info). Config files that share the same name as the template config file
will be renamed.
"""

import cwltool.factory
import subprocess
import shutil
import os

from pkg_resources import resource_filename
from aquatx.srna.Configuration import Configuration
from argparse import ArgumentParser


def get_args():
    """Parses command line input"""

    parser = ArgumentParser(description=__doc__)

    # Parser for subcommands: (run, setup-cwl, get-template, etc.)
    subparsers = parser.add_subparsers(required=True, dest='command')
    subcommands_with_configfile = {
        "run": "Processes the provided config file and executes the workflow it specifies.",
        "setup-cwl": 'Processes the provided config file and copies workflow files to the current directory',
        "setup-nextflow": "This subcommand is not yet implemented"
    }

    # Subcommands that require a configuration file argument
    for command, desc in subcommands_with_configfile.items():
        subparsers.add_parser(command).add_argument(
            '--config', metavar='configFile', required=True, help=desc
        )

    # Subcommand get-template has no additional arguments
    subparsers.add_parser("get-template",
                          help="Copies run config, sample, and reference templates to current directory")

    return parser.parse_args()


def run(aquatx_cwl_path: str, config_file: str) -> None:
    """Processes the provided config file and executes the workflow it defines

    The provided configuration file will be processed and rewritten to reflect the content
    of the sample and reference csv files. The location of these files is defined under the
    config file's ` samples_csv` and `reference_sheet_file` keys. Config files named
    "run_config_template.yml" will be left unmodified, and the processed config will
    instead be written under a file whose name reflects the current date and time.

    Args:
        aquatx_cwl_path: The path to the project's CWL workflow file directory
        config_file: The configuration file for this run.
    """

    print("Running the end-to-end analysis...")

    # First get the configuration file set up for this run
    config_object = Configuration(config_file)
    run_directory = config_object.create_run_directory()
    cwl_conf_file = config_object.write_processed_config()

    # Run with cwltool
    debug = False
    subprocess.run(f"cwltool --outdir {run_directory} --copy-outputs --timestamps "
                   f"{'--leave-tmpdir --debug --js-console ' if debug else ''}"
                   f"{aquatx_cwl_path}/workflows/aquatx_wf.cwl {cwl_conf_file}", shell=True)

    # runtime_context = cwltool.factory.RuntimeContext()
    # runtime_context.outdir = os.path.join('.', config.get('run_directory'))
    # runtime_context.on_error = "continue"
    #
    # loading_context = cwltool.factory.LoadingContext()
    # loading_context.jobdefaults = config.config
    #
    # cwl = cwltool.factory.Factory(runtime_context=runtime_context, loading_context=loading_context)
    # cwl.make(f"{aquatx_cwl_path}/workflows/aquatx_wf.cwl")


def get_template(aquatx_extras_path: str) -> None:
    """Retrieves the template run configuration file, and the sample/reference csv templates

    Args:
        aquatx_extras_path: The path to the project's extras directory. This directory
            contains templates for the run configuration, sample inputs, and reference
            inputs.
    """

    print("Copying template input files to current directory...")

    # Copy template files to the current working directory
    for template in ['run_config_template.yml', 'samples.csv', 'features.csv']:
        shutil.copyfile(f"{aquatx_extras_path}/{template}", f"{os.getcwd()}/{template}")


def setup_cwl(aquatx_cwl_path: str, config_file: str) -> None:
    """Retrieves the project's workflow files, and if provided, processes the run config file

    Args:
        aquatx_cwl_path: The path to the project's CWL workflow file directory
        config_file: The configuration file to be processed (or None/none to skip processing)

    """

    print("Creating cwl workflow...")

    # If the word "None" or "none" is supplied, simply copy workflow files. No config file processing.
    if config_file not in ('None', 'none'):
        # Set up the config file
        processed_config_location = Configuration(config_file).write_processed_config()
        print("The processed configuration file is located at: " + processed_config_location)

    # Copy the entire cwl directory to the current working directory
    shutil.copytree(aquatx_cwl_path, os.getcwd() + "/cwl/")
    print("The workflow and files are under: cwl/tools/ and cwl/workflows/")


def setup_nextflow(config_file: str) -> None:
    """This function is not yet implemented

    Args:
        config_file: The YML run configuration file to be converted for use with Nextflow

    """

    print("Creating nextflow workflow...")
    print("This command is currently not implemented.")


def main():
    """The main routine that determines what type of run to do.

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
