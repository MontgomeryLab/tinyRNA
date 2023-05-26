# For Users
These files are utilized during tinyRNA's setup. During normal usage you shouldn't need to access or modify them. To install or update the tinyRNA environment and codebase, navigate to the project's root directory and run `./setup.sh`

# For Developers
Environment lockfiles allow us to skip Conda's time-intensive preparation steps, e.g. "Solving environment" and "Collecting package metadata", when the user installs/updates the tinyRNA environment.

Dependencies are strictly pinned so that the user's environment exactly matches the environment where testing has been performed.

## Files in this Directory
| File                       | Description                                                                                                                                                                                                           |
|----------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `environment_3.[...].yml`  | Environment files for different Python versions. These files are platform-agnostic to simplify testing. Just one environment file is needed for creating both the Linux and macOS lockfiles.                          |
| `environment_template.yml` | Template environment file used for creating new environment files for a specific Python version.                                                                                                                      |
| `conda-[...]-64.lock`      | Platform-specific lockfiles that contain a fully solved dependency tree specification for the latest Python version we support. These are the files that the setup script uses when creating the tinyRNA environment. | 
| `env.py`                   | A script for creating and updating these files.                                                                                                                                                                       |

## The Development Environment
A separate conda environment is created to hold the dependencies required for creating and updating lockfiles. Keeping development and tinyRNA dependencies separate helps to avoid dependency conflicts.

## Updating the Default Python Version
An environment's Python version determines the environment's dependency tree. Creating specifications for a new Python version is a multi-step process, but a script has been provided to help automate it. **The script can be executed from any environment.**

The idea is to start with a template environment file containing all of our dependencies where nearly all dependencies are unpinned. Normally only the Python version is pinned in the template, but right now Pysam is also pinned because Mamba doesn't currently solve for the version we need by default. We then use Mamba to solve the unpinned dependency versions which are used to create a final, strictly pinned environment file for the new Python version. The new file is then used to create platform-specific lockfiles.

### 1. Create an Environment File for the new Python Version
```shell
python3 env.py new-env --py 3.XX --template environment_template.yml
```

Replace `3.XX` with the target Python version. This command will produce a new environment file named `environment_3.XX.yml` that contains each dependency's solved version.

### 2. Create Platform-Specific Lockfiles
```shell
python3 env.py update-locks --env-file environment_3.XX.yml
```

Replace `environment_3.XX.yml` with the name of the environment file you created in step 1. This command will rebuild the dependency tree in both lockfiles (linux and osx) according to the new environment file.



