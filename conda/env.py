import subprocess
import itertools
import argparse
import json
import sys
import os

from tempfile import TemporaryDirectory, NamedTemporaryFile
from typing import Dict

dev_env_name = "tinyrna_dev_v1.1"
dev_py_vers = "3.9"
dev_deps = [
    "conda=22.9",  # Only for compatibility with conda-lock=1.4.0 (see https://github.com/conda/conda-lock/issues/408)
    "python=" + dev_py_vers,
    "conda-forge::mamba=1.4.2",
    "conda-forge::conda-lock=1.4.0",
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

    update_lock = subcommands.add_parser("rebuild-locks", help="Rebuild osx-64 and linux-64 lockfiles from env-file.")
    update_lock.add_argument(
        "--env-file",
        help="The environment[...].yml to use for rebuilding lockfiles.",
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


def new_env_file(env_template, py_version, dev_env_path):
    """Creates a new environment file for the given Python version by solving the template with mamba and pinning"""

    # The dev environment's /bin first needs to be prepended to PATH in main()
    import ruamel.yaml

    yaml = ruamel.yaml.YAML()
    with open(env_template) as template, \
            TemporaryDirectory() as tmp_prefix, \
            NamedTemporaryFile('w+', suffix=".yml") as tmp_env_file:

        # Update Python version in the template contents and write it to a temporary file
        versioned_template = yaml.load(template.read().replace("{PYTHON_VERSION}", py_version))
        yaml.dump(versioned_template, tmp_env_file)
        tmp_env_file.flush()

        # Get the complete environment spec by solving the template with mamba
        solved_versions = mamba_solve(tmp_env_file, tmp_prefix)

    # Pin each unpinned dependency in the template with its solved version
    for i, spec in enumerate(versioned_template["dependencies"]):
        if set(spec) & {"=", "<", ">"}:  # Skip pinned template specs
            continue

        spec_no_channel = spec.split("::")[-1]
        version = solved_versions[spec_no_channel]
        versioned_template["dependencies"][i] = f"{spec}={version}"

    # Write the new environment file
    output_file = f"environment_{py_version}.yml"
    with open(output_file, 'w') as f:
        yaml.dump(versioned_template, f)
        print(f"Created environment file: {output_file}")


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


def main():
    args = get_args()
    denv = get_dev_env_path()

    # Import modules from dev env from here forward
    sys.path.insert(0, os.path.join(denv, "lib", f"python{dev_py_vers}", "site-packages"))

    if args.command == "new-env":
        new_env_file(args.template, args.py, denv)
    if args.command == "rebuild-locks":
        rebuild_lockfiles(args.env_file)

    print("Done.")


if __name__ == '__main__':
    main()