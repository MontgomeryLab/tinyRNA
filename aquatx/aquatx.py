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
word "None" or "none" to obtain only the workflow files). The config file
will be processed to generate pipeline settings from your input files.
"""

import cwltool.executors
import cwltool.factory
import cwltool.secrets
import subprocess
import shutil
import os

from cwltool.context import LoadingContext
from pkg_resources import resource_filename
from argparse import ArgumentParser

from aquatx.srna.Configuration import Configuration
from cwltool.main import setup_loadingContext


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

    The provided configuration file will be processed to reflect the contents of all
    related config files, and saved to the run directory. Both single and parallel
    library processing is supported.

    Args:
        aquatx_cwl_path: The path to the project's CWL workflow file directory
        config_file: The configuration file for this run

    """

    print("Running the end-to-end analysis...")

    # First get the configuration file set up for this run
    config_object = Configuration(config_file)
    run_directory = config_object.create_run_directory()
    cwl_conf_file = config_object.write_processed_config()

    debug = False
    if config_object['run_native']:  # experimental
        # Execute the CWL runner via native Python
        run_native(config_object, aquatx_cwl_path, run_directory,
                   debug=debug, parallel=config_object['run_parallel'])
    else:
        if config_object['run_parallel']:
            print("WARNING: parallel execution with cwltool is an experimental feature")

        # Use the cwltool CWL runner via command line
        cwl_runner = f"cwltool --outdir {run_directory} --copy-outputs --timestamps " \
                     f"{'--leave-tmpdir --debug --js-console ' if debug else ''}" \
                     f"{'--parallel ' if config_object['run_parallel'] else ''}" \
                     f"{aquatx_cwl_path}/workflows/aquatx_wf.cwl {cwl_conf_file}"

        subprocess.run(cwl_runner, shell=True)

        # This would be the invocation for Toil, but it only supports a CWL version that breaks
        #   when output names contain special characters. Keeping this here for future releases
        #
        # cwl_runner = f"toil-cwl-runner --outdir {run_directory} --stats " \
        #              f"{f'--debug-Worker --writeLogs {run_directory}' if debug else ''}" \
        #              f"{aquatx_cwl_path}/workflows/aquatx_wf.cwl {cwl_conf_file}"


def run_native(config_object, cwl_path, run_directory, debug=False, parallel=False):
    """Executes the workflow using native Python rather than subprocess "command line"

    Args:
        config_object: the processed configuration produced by Configuration.py
        cwl_path: the path to the project's CWL workflow file directory
        run_directory: the destination folder for workflow outputs
        debug: instruct the CWL runner to provide additional debug info
        parallel: process libraries in parallel where possible

    Returns: None

    """

    def furnish_if_file_record(file_dict):
        if isinstance(file_dict, dict) and file_dict.get('class', None) == 'File':
            file_dict['basename'] = os.path.basename(file_dict['path'])
            file_dict['location'] = file_dict['path']
            file_dict['contents'] = None

    for _, config_param in config_object.config.items():
        if isinstance(config_param, list):
            for config_dict in config_param:
                furnish_if_file_record(config_dict)
        else:
            furnish_if_file_record(config_param)

    runtime_context = cwltool.factory.RuntimeContext({
        'secret_store': cwltool.secrets.SecretStore(),
        'default_stdout': subprocess.PIPE,
        'default_stderr': subprocess.PIPE,
        'outdir': run_directory,
        'on_error': "continue",
        'debug': debug
    })

    loading_context = setup_loadingContext(LoadingContext(), runtime_context, {})

    cwl = cwltool.factory.Factory(
        runtime_context=runtime_context,
        loading_context=loading_context,
        executor=cwltool.executors.MultithreadedJobExecutor()   # Run jobs in parallel
        if parallel else cwltool.executors.SingleJobExecutor()  # Run one library at a time
    )

    pipeline = cwl.make(f"{cwl_path}/workflows/aquatx_wf.cwl")
    pipeline(**config_object.config)


def get_template(aquatx_extras_path: str) -> None:
    """Copies all configuration file templates to the current working directory

    Args:
        aquatx_extras_path: The path to the project's extras directory. This directory
            contains templates for the run configuration, sample inputs, feature selection
            rules, and paths for all the above.

    """

    print("Copying template input files to current directory...")

    # Copy template files to the current working directory
    for template in ['run_config_template.yml', 'samples.csv', 'features.csv', 'paths.yml']:
        shutil.copyfile(f"{aquatx_extras_path}/{template}", f"{os.getcwd()}/{template}")


def setup_cwl(aquatx_cwl_path: str, config_file: str) -> None:
    """Retrieves the project's workflow files, and if provided, processes the run config file

    Args:
        aquatx_cwl_path: The path to the project's CWL workflow file directory
        config_file: The configuration file to be processed (or None/none to skip processing)

    """

    # If the word "None" or "none" is supplied, simply copy workflow files. No config file processing.
    if config_file not in ('None', 'none'):
        # Set up the config file
        print("Processing configuration file...")
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
