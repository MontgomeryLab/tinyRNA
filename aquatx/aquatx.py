#!/usr/bin/env python
"""The main entry point for AQuATx for small RNA data analysis.

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
    aquatx <subcommand> --config <config-file>

The Run Config file should be supplied for the run subcommand (required)
and for the setup-cwl subcommand (optional; alternatively you may use the
word "None" or "none" to obtain only the workflow files). The config file
will be processed to generate pipeline settings from your input files.

The recount and replot commands must be invoked from within the Run Directory
of the run you wish to resume. The config file you supply must be the processed
Run Configuration file within the Run Directory.
"""

import cwltool.executors
import cwltool.factory
import cwltool.secrets
import subprocess
import shutil
import sys
import os

from cwltool.context import LoadingContext
from cwltool.utils import DEFAULT_TMP_PREFIX
from pkg_resources import resource_filename
from argparse import ArgumentParser

from aquatx.srna.Configuration import Configuration, ConfigBase
from aquatx.srna.resume import ResumeCounterConfig, ResumePlotterConfig


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


def run(aquatx_cwl_path: str, config_file: str) -> None:
    """Processes the provided config file and executes the workflow it defines

    The provided configuration file will be processed to reflect the contents of all
    related config files, and saved to the run directory. Both single and parallel
    library processing is supported.

    Args:
        aquatx_cwl_path: The path to the project's CWL workflow file directory
        config_file: The configuration file for this run

    Returns: None

    """

    print("Running the end-to-end analysis...")

    # First get the configuration file set up for this run
    config_object = Configuration(config_file)
    run_directory = config_object.create_run_directory()
    cwl_conf_file = config_object.save_run_profile()

    workflow = f"{aquatx_cwl_path}/workflows/aquatx_wf.cwl"
    parallel = config_object['run_parallel']
    debug = False

    if config_object['run_native']:  # experimental
        # Execute the CWL runner via native Python
        run_native(config_object, workflow, run_directory, debug=debug, parallel=parallel)
    else:
        if config_object['run_parallel']:
            print("WARNING: parallel execution with cwltool is an experimental feature")

        # Use the cwltool CWL runner via command line
        run_cwltool_subprocess(cwl_conf_file, workflow, run_directory=run_directory, parallel=parallel, debug=debug)


def resume(aquatx_cwl_path: str, config_file: str, step: str) -> None:
    """Resumes pipeline execution at either the Counter or Plotter step

    The user must invoke this from the RUN DIRECTORY for which they wish to
    resume, and the config file they provide must be the PROCESSED run config
    within that directory. The previous pipeline outputs and processed config
    will be used to resume a truncated workflow which begins at the chosen step.

    Output directories will be timestamped to keep results separate between
    resume runs. A copy of each run's Features Sheet will be provided in the
    (timestamped) counts output directory.

    Args:
        aquatx_cwl_path: The path to the project's CWL file directory
        config_file: The PROCESSED config file for the run to resume
        step: The stepname to serve as the new workflow entrypoint

    Returns: None

    """

    entry_config = {
        "Counter": ResumeCounterConfig,
        "Plotter": ResumePlotterConfig
    }

    print(f"Resuming pipeline execution at the {step} step...")
    config = entry_config[step](config_file, f"{aquatx_cwl_path}/workflows/aquatx_wf.cwl")
    resume_wf = f"{aquatx_cwl_path}/workflows/aquatx-resume.cwl"
    config.write_workflow(resume_wf)

    debug = False
    if config['run_native']:
        # Don't need to write processed config, pass in directly
        run_native(config, resume_wf, run_directory=".", debug=debug)
    else:
        resume_conf_file = "resume_" + os.path.basename(config_file)
        config.write_processed_config(resume_conf_file)
        run_cwltool_subprocess(resume_conf_file, resume_wf, debug=debug)

    if os.path.isfile(resume_wf):
        # We don't want the generated workflow to be returned by a call to setup-cwl
        os.remove(resume_wf)


def run_cwltool_subprocess(config_file: str, workflow: str, run_directory=None, parallel=False, debug=False) -> None:
    """Executes the workflow using a command line invocation of cwltool

    Args:
        config_file: the processed configuration file produced by Configuration.py
        workflow: the path to the workflow to be executed
        run_directory: the destination folder for workflow output subdirectories (default: CWD)
        parallel: process libraries in parallel where possible
        debug: instruct the CWL runner to provide additional debug info

    Returns: None

    """

    cwl_runner = "cwltool --copy-outputs --timestamps --relax-path-checks " \
                 f"{'--leave-tmpdir --debug --js-console ' if debug else ''}" \
                 f"{'--outdir ' + run_directory + ' ' if run_directory else ''}" \
                 f"{'--parallel ' if parallel else ''}" \
                 f"{workflow} {config_file}"

    subprocess.run(cwl_runner, shell=True)


def run_native(config_object: 'ConfigBase', workflow: str, run_directory: str = None, parallel=False, debug=False) -> None:
    """Executes the workflow using native Python rather than subprocess "command line"

    Args:
        config_object: a constructed ConfigBase-derived object
        workflow: the path to the workflow to be executed
        run_directory: the destination folder for workflow output subdirectories (default: CWD)
        parallel: process libraries in parallel where possible
        debug: instruct the CWL runner to provide additional debug info

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
        'move_outputs': "copy",
        'on_error': "continue",
        'debug': debug
    })

    if sys.platform == "darwin":
        default_mac_path = "/private/tmp/docker_tmp"
        if runtime_context.tmp_outdir_prefix == DEFAULT_TMP_PREFIX:
            runtime_context.tmp_outdir_prefix = default_mac_path
        if runtime_context.tmpdir_prefix == DEFAULT_TMP_PREFIX:
            runtime_context.tmpdir_prefix = default_mac_path

    cwl = cwltool.factory.Factory(
        runtime_context=runtime_context,
        loading_context=LoadingContext({'relax_path_checks': True}),
        executor=cwltool.executors.MultithreadedJobExecutor()   # Run jobs in parallel
        if parallel else cwltool.executors.SingleJobExecutor()  # Run one library at a time
    )

    pipeline = cwl.make(workflow)
    pipeline(**config_object.config)


def get_template(aquatx_extras_path: str) -> None:
    """Copies all configuration file templates to the current working directory

    Args:
        aquatx_extras_path: The path to the project's extras directory. This directory
            contains templates for the run configuration, sample inputs, feature selection
            rules, the project's matplotlib stylesheet, and paths for all the above.

    Returns: None

    """

    print("Copying template input files to current directory...")

    template_files = ['run_config_template.yml', 'samples.csv', 'features.csv',
                      'paths.yml', 'aquatx-srna-light.mplstyle']

    # Copy template files to the current working directory
    for template in template_files:
        shutil.copyfile(f"{aquatx_extras_path}/{template}", f"{os.getcwd()}/{template}")


def setup_cwl(aquatx_cwl_path: str, config_file: str) -> None:
    """Retrieves the project's workflow files, and if provided, processes the run config file

    Args:
        aquatx_cwl_path: The path to the project's CWL workflow file directory
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
        recount: Resume pipeline execution at the Counter step
        replot: Resume pipeline execution at the Plotter step
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
        "replot": lambda: resume(aquatx_cwl_path, args.config, "Plotter"),
        "recount": lambda: resume(aquatx_cwl_path, args.config, "Counter"),
        "setup-cwl": lambda: setup_cwl(aquatx_cwl_path, args.config),
        "get-template": lambda: get_template(aquatx_extras_path),
        "setup-nextflow": lambda: setup_nextflow(args.config)
    }

    command_map[args.command]()


if __name__ == '__main__':
    main()
