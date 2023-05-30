"""This script helps manage conda environment specifications for the project.

A separate development environment is created to hold its dependencies
(Python modules and command-line utilities). The script can be executed in any environment
and will utilize the dependencies in the dev environment via `conda run` and modification
of the system PATH variable."""
import platform
import subprocess
import itertools
import argparse
import json
import sys
import os

## ADDITIONAL IMPORTS from dev environment after calling redirect_imports() ##

from tempfile import TemporaryDirectory, NamedTemporaryFile
from typing import Dict


osx_lockfile = "conda-osx-64.lock"
linux_lockfile = "conda-linux-64.lock"

dev_env_name = "tinyrna_dev_v1.1"
dev_py_vers = "3.9"
dev_deps = [
    "conda=22.9",  # Only for compatibility with conda-lock=1.4.0 (see https://github.com/conda/conda-lock/issues/408)
    "python=" + dev_py_vers,
    "mamba=1.4.2",
    "conda-lock=1.4.0",
    "jake=3.0.0",
    "conda-pack=0.7.0"
]


def get_args():
    parser = argparse.ArgumentParser()
    subcommands = parser.add_subparsers(
        title="Available subcommands",
        description=f"Use `python3 {os.path.basename(__file__)} <subcommand> -h` for subcommand details.",
        required=True,
        dest="command"
    )

    new_env = subcommands.add_parser("new-env", help="Create an environment.yml file for a new Python version.")
    new_env.add_argument(
        "--py",
        help='The target Python version formatted as "3.XX".',
        required=True
    )
    new_env.add_argument(
        "--template",
        help="The environment file to serve as a template for the new file. Dependencies should be minimally pinned.",
        required=True
    )

    rebuild_locks = subcommands.add_parser("rebuild-locks", help="Rebuild osx-64 and linux-64 lockfiles from env-file.")
    rebuild_locks.add_argument(
        "--env-file",
        help="The environment[...].yml to use for rebuilding lockfiles.",
        required=True
    )

    make_archive = subcommands.add_parser("make-archive", help="Create an environment archive for offline installation.",
                                          description="NOTE: archive is OS-specific. Install using conda-pack.")
    make_archive.add_argument(
        "--env-name",
        help="The name of the environment to archive.",
        required=True
    )

    return parser.parse_args()


def get_conda_envs() -> Dict[str, str]:
    """Returns a dict of {env_name: env_path} for all conda environments."""

    proc = subprocess.run(
        ["conda", "env", "list"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8"
    )

    # Check for errors
    try:
        proc.check_returncode()
        raw_env_list = proc.stdout
    except subprocess.CalledProcessError:
        print(proc.stderr, file=sys.stderr)
        sys.exit("Failed to get conda environment list.")

    envs = [line.replace('*', ' ').split(' ', 1)
            for line in raw_env_list.splitlines()
            if len(line) and line[0] not in ("#", " ")]

    return {env_name.strip(): env_path.strip()
            for env_name, env_path in envs}


def get_dev_env_path() -> str:
    """Returns the path to the dev environment, creating it if necessary."""

    print("Checking for dev environment ...")

    envs = get_conda_envs()
    if dev_env_name in envs:
        print("Dev environment found.")
        return envs[dev_env_name]

    with open("dev_env_create.log", 'w') as log:
        print(f"Creating dev environment {dev_env_name} ...")
        proc = subprocess.run(
            ["conda", "create",
             "--name", dev_env_name,
             "--channel", "conda-forge",
             "--yes",  *dev_deps],
            stdout=log, stderr=log,
            encoding="utf-8"
        )

    try:
        proc.check_returncode()
        print("Dev environment created.")
    except subprocess.CalledProcessError:
        sys.exit("Failed to create dev environment. See dev_env_create.log for details.")

    return get_conda_envs()[dev_env_name]


class EnvSolver:

    def __init__(self, template_filename):
        self.template_filename = template_filename

        import ruamel.yaml  # redirect_inputs() must be called first
        self.yaml = ruamel.yaml.YAML()

    def solve(self, py_version: str):
        py_bounded = self.get_bounded_py_spec(py_version)
        template_yml = self.get_template(py_bounded)
        solved_versions = self.solve_specs(template_yml)
        solved_template = self.get_updated_template(template_yml, solved_versions)

        # Write the new environment file
        output_file = f"environment_{py_version}.yml"
        with open(output_file, 'w') as f:
            self.yaml.dump(solved_template, f)
            print(f"Created environment file: {output_file}")

    def get_template(self, py_bounded: str) -> dict:
        """Returns a parsed template for the bounded Python version"""

        with open(self.template_filename, 'r') as f:
            versioned = f.read().replace("{PYTHON_VERSION}", py_bounded)
            return self.yaml.load(versioned)

    def solve_specs(self, template_yml: dict) -> dict:
        """Returns the solved dependency tree as a dictionary of {package_name: package_version}"""

        with NamedTemporaryFile('w+', suffix=".yml") as tmp_env_file, \
                TemporaryDirectory() as tmp_prefix:

            # Write Python-versioned template to a temporary file
            self.yaml.dump(template_yml, tmp_env_file)
            tmp_env_file.flush()

            # Get the complete environment spec by solving the template with mamba
            solved_versions = mamba_solve(tmp_env_file, tmp_prefix)

        return solved_versions

    def get_updated_template(self, template_yml: dict, solutions: dict) -> dict:
        """Updates each unpinned dependency in the template with its solved version"""

        for i, spec in enumerate(template_yml["dependencies"]):
            if set(spec) & {"=", "<", ">"}:  # Skip pinned specs
                continue

            spec_no_channel = spec.split("::")[-1]
            version = solutions[spec_no_channel]
            template_yml["dependencies"][i] = f"{spec}={version}"

        return template_yml

    def get_bounded_py_spec(self, py_version: str) -> str:
        """Place an upper bound on the specified version to prevent the solver from overshooting"""
        v_split = py_version.split('.')
        return f"{py_version}, <{v_split[0]}.{int(v_split[1]) + 1}"


def mamba_solve(env_file, env_prefix) -> Dict[str, str]:
    """Calls Mamba subprocess to solve the env_file specs and returns the result as a dict
    of {package_name: package_version}.

    This function is a reduced version of conda_lock.conda_solver.solve_specs_for_arch()."""

    # Call mamba from the dev environment regardless of the current environment
    env_ctx =   ["conda", "run", "-n", dev_env_name]
    solve_cmd = ["mamba", "env", "create",
                 "--file", env_file.name,
                 "--prefix", env_prefix,
                 "--dry-run",
                 "--json"]

    print("Solving target environment ...")
    proc = subprocess.run(
        env_ctx + solve_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8"
    )

    # Check for errors
    try:
        proc.check_returncode()
    except subprocess.CalledProcessError:
        try:
            err_json = json.loads(proc.stdout)
            message = err_json.get("message", err_json)
        except json.JSONDecodeError as e:
            print(f"Failed to parse json, {e}", file=sys.stderr)
            message = proc.stdout
        sys.exit(message)

    # Parse the output
    proc_stdout = proc.stdout
    json_resp = json.loads(proc_stdout[proc_stdout.index("{"): proc_stdout.rindex("}") + 1])
    actions = json_resp["actions"]

    # Return a dict of {package_name: package_version}
    return {
        pkg["name"]: pkg["version"]
        for pkg in itertools.chain(actions['FETCH'], actions['LINK'])
    }

def rebuild_lockfiles(env_file):
    """Rebuilds the osx-64 and linux-64 lockfiles from the provided environment file."""

    supported_platforms = ["osx-64", "linux-64"]
    env_ctx = ["conda", "run", "-n", dev_env_name]  # Run command from dev environment regardless of current env

    for platform in supported_platforms:
        lock_cmd = ["conda-lock",
                    "--file", env_file,
                    "--kind", "explicit",
                    "--platform", platform]

        print(f"Building {platform} lockfile ...")
        proc = subprocess.run(
            env_ctx + lock_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            encoding="utf-8"
        )

        try:
            proc.check_returncode()
        except subprocess.CalledProcessError:
            log_file = f"{platform}_lockfile.log"
            with open(log_file, 'w') as f:
                f.write(proc.stdout)
            sys.exit(f"Error occurred while rebuilding the {platform} lockfile. See {log_file}.")

        print(f"{platform} lockfile complete.")


def scan_lockfiles():
    """Scans for vulnerable packages in both lockfiles using the Sonatype OSS Index"""

    env_ctx =  ["conda", "run", "--no-capture-output", "-n", dev_env_name]  # Run command from dev environment
    scan_cmd = ["jake", "ddt", "-t", "CONDA"]

    for lock in (osx_lockfile, linux_lockfile):
        with open(lock, 'rb') as f:
            packages = f.read()

        subprocess.run(
            env_ctx + scan_cmd,
            input=packages
        )

        # Indicate which file the report was for
        print('\n'.join(["", "=" * 50, f"{lock.upper()}: SCAN COMPLETE", "=" * 50]))


def make_env_archive(env_name):
    """Just a wrapper for conda-pack."""

    reported_os = platform.system()
    os_name = {'Darwin': "osx"}.get(reported_os, reported_os).lower()
    archive_name = f"{env_name}_{os_name}.tar.gz"

    env_ctx = ["conda", "run", "--no-capture-output", "-n", dev_env_name]  # Run command from dev environment regardless of current env
    command = ["conda-pack",
               "-n", env_name,
               "-o", archive_name]

    subprocess.run(env_ctx + command)


def redirect_imports(denv):
    """Adds the dev environment's Python site-packages location to $PATH so that
    dependencies for this script can be imported from it. This allows the script
    to be executed from any environment. It also allows us to create the dev
    environment and utilize its dependencies all in the same run."""

    sys.path.insert(0, os.path.join(denv, "lib", f"python{dev_py_vers}", "site-packages"))


def main():
    args = get_args()
    denv = get_dev_env_path()
    redirect_imports(denv)

    if args.command == "new-env":
        EnvSolver(args.template).solve(args.py)
    if args.command == "rebuild-locks":
        rebuild_lockfiles(args.env_file)
        scan_lockfiles()
    if args.command == "make-archive":
        make_env_archive(args.env_name)

    print("Done.")


if __name__ == '__main__':
    main()