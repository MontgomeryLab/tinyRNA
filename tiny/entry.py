#!/usr/bin/env python
"""tinyRNA provides an all-in-one workflow for precision analysis of sRNA-seq data.
At the core of tinyRNA is a highly flexible counting utility, tiny-count, that allows
for hierarchical assignment of reads to features based on their attributes.
"""

import cwltool.factory
import cwltool.secrets
import coloredlogs
import subprocess
import functools
import argparse
import logging
import shutil
import sys
import os

from cwltool.context import LoadingContext, RuntimeContext
from cwltool.executors import SingleJobExecutor, MultithreadedJobExecutor
from cwltool.utils import DEFAULT_TMP_PREFIX
from pkg_resources import resource_filename

from tiny.rna.configuration import Configuration, ConfigBase, get_templates
from tiny.rna.resume import ResumeCounterConfig, ResumePlotterConfig
from tiny.rna.util import report_execution_time, SmartFormatter, add_transparent_help

# Global variables
cwl_path = resource_filename('tiny', 'cwl')
templates_path = resource_filename('tiny', 'templates')
template_files = [os.path.basename(f) for f in os.listdir(templates_path)]
optional_files = [f for f in template_files if f.endswith('.mplstyle')]


def get_args():
    """Parses command line input"""

    primary_parser = argparse.ArgumentParser(
        description=__doc__,         # Borrow docstring found at the beginning of the file
        add_help=False,              # We later add --help argument "transparently" so it isn't listed in help string
        conflict_handler='resolve',  # Otherwise argparse gets upset about subparsers also having "transparent" help
    )

    add_transparent_help(primary_parser)

    subcommands = primary_parser.add_subparsers(
        title="Available subcommands",
        description="Run tiny SUBCOMMAND -h for subcommand details",
        required=True,
        dest='command',
    )

    def add_subcommand(command, brief, detail, config: str = None, **kwargs):
        parser = subcommands.add_parser(
            command,
            help=brief,
            description=detail,
            conflict_handler='resolve',
            add_help=False,
            **kwargs
        )
        add_transparent_help(parser)
        if config is not None:
            args_section = parser.add_argument_group("Required arguments")
            args_section.add_argument(
                '--config',
                help=f"The file path to your {config.replace('_', ' ')}",
                metavar=config.upper(),
                required=True
            )

    def priority_sort(x):
        priority = ('.yml', '.csv', '.mplstyle')
        ext = os.path.splitext(x)[1]
        return priority.index(ext) if ext in priority else sys.maxsize

    ext_sorted_templates = sorted(template_files, key=priority_sort)
    add_subcommand(
        "get-templates",
        brief="Copy configuration files to the current directory",
        detail=f"This command copies configuration file templates to the current directory. "
               "These files include: \n\t" + '\n\t'.join(ext_sorted_templates) + "\n\n"
               "The following files are not required for pipeline execution:\n  "
               + ', '.join(optional_files),
        formatter_class=SmartFormatter
    )

    add_subcommand(
        "run",
        brief="Run an end-to-end analysis",
        detail="This command coordinates a comprehensive end-to-end analysis of your input files "
               "according to the preferences defined in your configuration files.\n\n"
               "The primary configuration file, or Run Config, contains user preferences "
               "for each step of the workflow, and a reference to a second configuration "
               "file: the Paths File. Within the Paths File you must specify the location "
               "of your file inputs, including two final configuration files: the Samples "
               "Sheet and Features Sheet.",
        config='run_config',
        formatter_class=SmartFormatter
    )

    resume_helpstring_template = (
        "This command resumes the workflow at the {step_name} step using outputs "
        "from a prior end-to-end analysis. It must be executed within the run "
        "directory whose prior outputs you wish to reuse. Additionally, you must "
        "provide the processed Run Config which is located in that directory.")

    add_subcommand(
        "recount",
        brief="Resume an analysis at the tiny-count step",
        detail=resume_helpstring_template.format(step_name='tiny-count'),
        config='processed_run_config'
    )

    add_subcommand(
        "replot",
        brief="Resume an analysis at the tiny-plot step",
        detail=resume_helpstring_template.format(step_name='tiny-plot'),
        config='processed_run_config'
    )

    add_subcommand(
        "setup-cwl",
        brief="Produce a processed Run Config and workflow files",
        detail="This command should only be used if you intend to use an alternative "
               "workflow runner (advanced). It will produce a processed copy of your "
               "Run Config, and it will copy the workflow CWL files to the current "
               'directory. You can also pass the word "none" for your Run Config to '
               "skip processing and only copy the workflow files.",
        config='run_config'
    )

    return primary_parser.parse_args()


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
    config_object = Configuration(config_file, validate_inputs=True)
    run_directory = config_object.create_run_directory()
    config_object.save_run_profile()

    workflow = f"{tinyrna_cwl_path}/workflows/tinyrna_wf.cwl"

    if config_object['run_native']:
        # Execute the CWL runner via native Python
        return_code = run_cwltool_native(config_object, workflow, run_directory)
    else:
        # Use the cwltool CWL runner via command line
        return_code = run_cwltool_subprocess(config_object, workflow, run_directory)

    config_object.execute_post_run_tasks(return_code)


@report_execution_time("Pipeline resume runtime")
def resume(tinyrna_cwl_path: str, config_file: str, step: str) -> None:
    """Resumes pipeline execution at either the tiny-count or tiny-plot step

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
        "tiny-count": ResumeCounterConfig,
        "tiny-plot": ResumePlotterConfig
    }

    print(f"Resuming pipeline execution at the {step} step...")

    # Make appropriate config and workflow for this step; write modified workflow to disk
    config = entry_config[step](config_file, f"{tinyrna_cwl_path}/workflows/tinyrna_wf.cwl")
    resume_wf = f"{tinyrna_cwl_path}/workflows/tiny-resume.cwl"
    config.write_workflow(resume_wf)

    if config['run_native']:
        # We can pass our config object directly without writing to disk first
        run_cwltool_native(config, resume_wf)
    else:
        # Processed Run Config must be written to disk first
        resume_conf_file = config.get_outfile_path()
        config.write_processed_config(resume_conf_file)
        run_cwltool_subprocess(config, resume_wf)

    if os.path.isfile(resume_wf):
        # We don't want the generated workflow to be returned by a call to setup-cwl
        os.remove(resume_wf)


def run_cwltool_subprocess(config_object: 'ConfigBase', workflow: str, run_directory='.') -> int:
    """Executes the workflow using a command line invocation of cwltool

    Args:
        config_object: a constructed ConfigBase-derived object
        workflow: the path to the workflow to be executed
        run_directory: the destination folder for workflow output subdirectories (default: CWD)

    Returns: None

    """

    processed_configfile = config_object.get_outfile_path()
    tmpdir = config_object['tmp_directory']
    parallel = config_object['run_parallel']
    verbosity = config_object['verbosity']

    command = ['cwltool --timestamps --relax-path-checks --on-error continue']
    command.append(f'--log-dir {run_directory}/{config_object["dir_name_logs"]}')
    command.append(f'--outdir {run_directory}')

    if tmpdir is not None: command.append(f'--tmpdir-prefix {tmpdir}')
    if verbosity == 'debug': command.append('--debug --js-console --leave-tmpdir')
    if verbosity == 'quiet': command.append('--quiet')
    if parallel: command.append('--parallel')

    cwl_runner = ' '.join(command + [workflow, processed_configfile])
    return subprocess.run(cwl_runner, shell=True).returncode


def run_cwltool_native(config_object: 'ConfigBase', workflow: str, run_directory: str = '.') -> int:
    """Executes the workflow using native Python rather than subprocess "command line"

    Args:
        config_object: a constructed ConfigBase-derived object
        workflow: the path to the workflow to be executed
        run_directory: the destination folder for workflow output subdirectories (default: CWD)

    Returns: None

    """

    tmpdir = config_object['tmp_directory']
    parallel = config_object['run_parallel']
    verbosity = config_object['verbosity']

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
        'log_dir': f"{run_directory}/{config_object['dir_name_logs']}",
        'secret_store': cwltool.secrets.SecretStore(),
        'outdir': run_directory,
        'on_error': "continue",
        'js_console': verbosity == "debug",
        'debug': verbosity == "debug"
    })

    # Set temp directory for intermediate outputs
    if tmpdir is not None:
        runtime_context.tmp_outdir_prefix = tmpdir
    elif sys.platform == "darwin":
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
    parallel_exec: MultithreadedJobExecutor = functools.partial(MultithreadedJobExecutor(), logger=logger)
    serial_exec: SingleJobExecutor = functools.partial(SingleJobExecutor(), logger=logger)

    # Instantiate Factory with our run preferences
    cwl = cwltool.factory.Factory(
        runtime_context=runtime_context,
        loading_context=LoadingContext({'relax_path_checks': True}),
        executor=parallel_exec if parallel else serial_exec
    )

    try:
        # Load the workflow document and execute
        pipeline = cwl.make(workflow)
        pipeline(**config_object.config)
    except cwltool.factory.WorkflowStatus:
        # For now, return non-zero if workflow did not complete
        return 1

    return 0


def setup_cwl(tinyrna_cwl_path: str, config_file: str) -> None:
    """Retrieves the project's workflow files, and if provided, processes the run config file

    Args:
        tinyrna_cwl_path: The path to the project's CWL workflow file directory
        config_file: The configuration file to be processed (or none to skip processing)

    Returns: None

    """

    if config_file.upper() == "NONE":
        # Simply copy workflow files. No config file processing.
        print("Processing configuration file...")
        outfile_name = "processed_" + os.path.basename(config_file)
        Configuration(config_file).write_processed_config(filename=outfile_name)
        print("The processed configuration file is located at: " + outfile_name)

    # Copy the entire cwl directory to the current working directory
    shutil.copytree(tinyrna_cwl_path, os.getcwd() + "/cwl/")
    print("The workflow and files are under: cwl/tools/ and cwl/workflows/")


def main():
    """The main routine that determines what type of run to do.

    Options:
        run: Run an end-to-end analysis based on a config file.
        recount: Resume pipeline execution at the tiny-count step
        replot: Resume pipeline execution at the tiny-plot step
        get-templates: Get the input sheets & template config files.
        setup-cwl: Get the CWL workflow for a run
    """

    # Parse command line arguments
    args = get_args()

    # Execute appropriate command based on command line input
    command_map = {
        "run": lambda: run(cwl_path, args.config),
        "replot": lambda: resume(cwl_path, args.config, "tiny-plot"),
        "recount": lambda: resume(cwl_path, args.config, "tiny-count"),
        "setup-cwl": lambda: setup_cwl(cwl_path, args.config),
        "get-templates": lambda: get_templates("tiny")
    }

    command_map[args.command]()


if __name__ == '__main__':
    main()
