#!/usr/bin/env python
"""The main entry point for tinyRNA for small RNA data analysis.

This tool provides an end-to-end workflow for analyzing small RNA sequencing
data from raw fastq files. This entry point also provides options for only
returning template files and workflows that can be used separately.

Subcommands:
    - get-template
    - setup-cwl
    - recount
    - replot
    - run

When installed, run, recount and setup-cwl should be invoked with:
    tiny <subcommand> --config <config-file>

The Run Config file should be supplied for the run subcommand (required)
and for the setup-cwl subcommand (optional; alternatively you may use the
word "None" or "none" to obtain only the workflow files). The config file
will be processed to generate pipeline settings from your input files.

The recount and replot commands must be invoked from within the Run Directory
of the run you wish to resume. The config file you supply must be the processed
Run Configuration file within the Run Directory.
"""

import cwltool.factory
import cwltool.secrets
import coloredlogs
import subprocess
import functools
import logging
import shutil
import sys
import os

from cwltool.context import LoadingContext, RuntimeContext
from cwltool.executors import SingleJobExecutor, MultithreadedJobExecutor
from cwltool.utils import DEFAULT_TMP_PREFIX
from pkg_resources import resource_filename
from argparse import ArgumentParser

from tiny.rna.Configuration import Configuration, ConfigBase
from tiny.rna.resume import ResumeCounterConfig, ResumePlotterConfig
from tiny.rna.util import report_execution_time


def get_args():
    """Parses command line input"""

    parser = ArgumentParser(description=__doc__)

    # Parser for subcommands: (run, recount, replot, setup-cwl, get-template, etc.)
    subparsers = parser.add_subparsers(required=True, dest='command')
    subcommands_with_configfile = {
        "run": "Processes the provided config file and executes the workflow it specifies.",
        "replot": "Resume pipeline at the Plotter step using the PROCESSED run config provided",
        "recount": "Resume pipeline at the Counter step using the PROCESSED run config provided",
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


@report_execution_time("Pipeline runtime")
def run(tinyrna_cwl_path: str, config_file: str) -> None:
    """Processes the provided config file and executes the workflow it defines

    The provided configuration file will be processed to reflect the contents of all
    related config files, and saved to the run directory. Both single and parallel
    library processing is supported.

    Args:
        tinyrna_cwl_path: The path to the project's CWL workflow file directory
        config_file: The configuration file for this run

    Returns: None

    """

    print("Running the end-to-end analysis...")

    # First get the configuration file set up for this run
    config_object = Configuration(config_file)
    run_directory = config_object.create_run_directory()
    cwl_conf_file = config_object.save_run_profile()

    workflow = f"{tinyrna_cwl_path}/workflows/tinyrna_wf.cwl"
    parallel = config_object['run_parallel']
    loudness = config_object['verbosity']

    if config_object['run_native']:  # experimental
        # Execute the CWL runner via native Python
        return_code = run_native(config_object, workflow, run_directory, verbosity=loudness)
    else:
        if config_object['run_parallel']:
            print("WARNING: parallel execution with cwltool is an experimental feature")

        # Use the cwltool CWL runner via command line
        return_code = run_cwltool_subprocess(
            cwl_conf_file, workflow,
            run_directory=run_directory,
            parallel=parallel, verbosity=loudness)

    # If the workflow completed without errors, we want to update
    # the Paths Sheet to point to the new bowtie index prefix
    if config_object['run_bowtie_build'] and return_code == 0:
        paths_sheet_filename = config_object.paths.inf
        config_object.paths.write_processed_config(paths_sheet_filename)


@report_execution_time("Pipeline resume runtime")
def resume(tinyrna_cwl_path: str, config_file: str, step: str) -> None:
    """Resumes pipeline execution at either the Counter or Plotter step

    The user must invoke this from the RUN DIRECTORY for which they wish to
    resume, and the config file they provide must be the PROCESSED run config
    within that directory. The previous pipeline outputs and processed config
    will be used to resume a truncated workflow which begins at the chosen step.

    Output directories will be timestamped to keep results separate between
    resume runs. A copy of each run's Features Sheet will be provided in the
    (timestamped) counts output directory.

    Args:
        tinyrna_cwl_path: The path to the project's CWL file directory
        config_file: The PROCESSED config file for the run to resume
        step: The stepname to serve as the new workflow entrypoint

    Returns: None

    """

    # Maps step to Configuration class
    entry_config = {
        "Counter": ResumeCounterConfig,
        "Plotter": ResumePlotterConfig
    }

    print(f"Resuming pipeline execution at the {step} step...")

    # Make appropriate config and workflow for this step; write modified workflow to disk
    config = entry_config[step](config_file, f"{tinyrna_cwl_path}/workflows/tinyrna_wf.cwl")
    resume_wf = f"{tinyrna_cwl_path}/workflows/tiny-resume.cwl"
    config.write_workflow(resume_wf)

    if config['run_native']:
        # We can pass our config object directly without writing to disk first
        run_native(config, resume_wf, verbosity=config['verbosity'])
    else:
        # Processed Run Config must be written to disk first
        resume_conf_file = "resume_" + os.path.basename(config_file)
        config.write_processed_config(resume_conf_file)
        run_cwltool_subprocess(resume_conf_file, resume_wf, verbosity=config['verbosity'])

    if os.path.isfile(resume_wf):
        # We don't want the generated workflow to be returned by a call to setup-cwl
        os.remove(resume_wf)


def run_cwltool_subprocess(config_file: str, workflow: str, run_directory=None, parallel=False, verbosity="normal") -> int:
    """Executes the workflow using a command line invocation of cwltool

    Args:
        config_file: the processed configuration file produced by Configuration.py
        workflow: the path to the workflow to be executed
        run_directory: the destination folder for workflow output subdirectories (default: CWD)
        parallel: process libraries in parallel where possible
        verbosity: controls the depth of information written to terminal by cwltool

    Returns: None

    """

    command = ['cwltool --copy-outputs --timestamps --relax-path-checks']
    if verbosity == 'debug': command.append('--debug --js-console')
    if verbosity == 'quiet': command.append('--quiet')
    if run_directory: command.append(f'--outdir {run_directory}')
    if parallel: command.append('--parallel')

    cwl_runner = ' '.join(command + [workflow, config_file])
    return subprocess.run(cwl_runner, shell=True).returncode


def run_native(config_object: 'ConfigBase', workflow: str, run_directory: str = '.', verbosity="normal") -> int:
    """Executes the workflow using native Python rather than subprocess "command line"

    Args:
        config_object: a constructed ConfigBase-derived object
        workflow: the path to the workflow to be executed
        run_directory: the destination folder for workflow output subdirectories (default: CWD)
        parallel: process libraries in parallel where possible
        verbosity: controls the depth of information written to terminal by cwltool

    Returns: None

    """

    def furnish_if_file_record(file_dict):
        if isinstance(file_dict, dict) and file_dict.get('class', None) == 'File':
            file_dict['basename'] = os.path.basename(file_dict['path'])
            file_dict['location'] = file_dict['path']
            file_dict['contents'] = None

    # Upgrade file entries in Run Config with extra descriptors cwltool expects
    for _, config_param in config_object.config.items():
        if isinstance(config_param, list):
            for config_dict in config_param:
                furnish_if_file_record(config_dict)
        else:
            furnish_if_file_record(config_param)

    # Set overall config for cwltool
    runtime_context = RuntimeContext({
        'secret_store': cwltool.secrets.SecretStore(),
        'outdir': run_directory,
        'move_outputs': "copy",
        'on_error': "continue",
        'js_console': verbosity == "debug",
        'debug': verbosity == "debug"
    })

    # Set proper temp directory for Mac users
    if sys.platform == "darwin":
        default_mac_path = "/private/tmp/docker_tmp"
        if runtime_context.tmp_outdir_prefix == DEFAULT_TMP_PREFIX:
            runtime_context.tmp_outdir_prefix = default_mac_path
        if runtime_context.tmpdir_prefix == DEFAULT_TMP_PREFIX:
            runtime_context.tmpdir_prefix = default_mac_path

    # Enable rich terminal output (timestamp, color, formatting)
    logger = logging.getLogger("cwltool")
    logger.handlers.clear()  # executors.py loads a default handler; outputs are printed twice if we don't clear it
    level = 'DEBUG' if verbosity == 'debug' else 'WARN' if verbosity == "quiet" else "INFO"
    coloredlogs.install(logger=logger, stream=sys.stderr, fmt="[%(asctime)s] %(levelname)s %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S", level=level, isatty=True)

    # Create a wrapper for the executors so that we may pass our logger to them (unsupported by Factory)
    parallel: MultithreadedJobExecutor = functools.partial(MultithreadedJobExecutor(), logger=logger)
    serial: SingleJobExecutor = functools.partial(SingleJobExecutor(), logger=logger)

    # Instantiate Factory with our run preferences
    cwl = cwltool.factory.Factory(
        runtime_context=runtime_context,
        loading_context=LoadingContext({'relax_path_checks': True}),
        executor=parallel if parallel else serial
    )

    try:
        # Load the workflow document and execute
        pipeline = cwl.make(workflow)
        pipeline(**config_object.config)
    except cwltool.factory.WorkflowStatus:
        # For now, return non-zero if workflow did not complete
        return 1

    return 0


def get_template(tinyrna_extras_path: str) -> None:
    """Copies all configuration file templates to the current working directory

    Args:
        tinyrna_extras_path: The path to the project's extras directory. This directory
            contains templates for the run configuration, sample inputs, feature selection
            rules, the project's matplotlib stylesheet, and paths for all the above.

    Returns: None

    """

    print("Copying template input files to current directory...")

    template_files = ['run_config_template.yml', 'samples.csv', 'features.csv',
                      'paths.yml', 'tinyrna-light.mplstyle']

    # Copy template files to the current working directory
    for template in template_files:
        shutil.copyfile(f"{tinyrna_extras_path}/{template}", f"{os.getcwd()}/{template}")


def setup_cwl(tinyrna_cwl_path: str, config_file: str) -> None:
    """Retrieves the project's workflow files, and if provided, processes the run config file

    Args:
        tinyrna_cwl_path: The path to the project's CWL workflow file directory
        config_file: The configuration file to be processed (or None/none to skip processing)

    Returns: None

    """

    # If the word "None" or "none" is supplied, simply copy workflow files. No config file processing.
    if config_file not in ('None', 'none'):
        # Set up the config file
        print("Processing configuration file...")
        outfile_name = "processed_" + os.path.basename(config_file)
        Configuration(config_file).write_processed_config(filename=outfile_name)
        print("The processed configuration file is located at: " + outfile_name)

    # Copy the entire cwl directory to the current working directory
    shutil.copytree(tinyrna_cwl_path, os.getcwd() + "/cwl/")
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
        recount: Resume pipeline execution at the Counter step
        replot: Resume pipeline execution at the Plotter step
        get-template: Get the input sheets & template config files.
        setup-cwl: Get the CWL workflow for a run
    """

    # Parse command line arguments
    args = get_args()

    # Get the package data
    tinyrna_cwl_path = resource_filename('tiny', 'cwl')
    tinyrna_extras_path = resource_filename('tiny', 'extras')

    # Execute appropriate command based on command line input
    command_map = {
        "run": lambda: run(tinyrna_cwl_path, args.config),
        "replot": lambda: resume(tinyrna_cwl_path, args.config, "Plotter"),
        "recount": lambda: resume(tinyrna_cwl_path, args.config, "Counter"),
        "setup-cwl": lambda: setup_cwl(tinyrna_cwl_path, args.config),
        "get-template": lambda: get_template(tinyrna_extras_path),
        "setup-nextflow": lambda: setup_nextflow(args.config)
    }

    command_map[args.command]()


if __name__ == '__main__':
    main()
