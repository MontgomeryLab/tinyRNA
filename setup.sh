#!/usr/bin/env bash

# USAGE: if you would like to install under a different conda environment name,
# you can pass the preferred name as the first argument to the script:
#   ./setup.sh preferred_name

env_name=${1:-tinyrna}
miniconda_version="23.3.1-0"
cwd="$(dirname "$0")"

# This is the default Python version that will be used by Miniconda (if Miniconda requires installation).
# Note that this isn't the same as the tinyRNA environment's Python version.
# The tinyRNA environment's Python version is instead specified in the platform lockfile.
miniconda_python_version="310"

function success() {
  check="✓"
  green_on="\033[1;32m"
  green_off="\033[0m"
  printf "${green_on}${check} %s${green_off}\n" "$*"
}

function status() {
  blue_on="\033[1;34m"
  blue_off="\033[0m"
  printf "${blue_on}%s${blue_off}\n" "$*"
}

function fail() {
  nope="⃠"
  red_on="\033[1;31m"
  red_off="\033[0m"
  printf "${red_on}${nope} %s${red_off}\n" "$*"
}

function stop() {
  kill -TERM -$$
}
# Ensures that when user presses Ctrl+C, the script stops
# rather than stopping current task and proceeding to the next
trap 'stop' SIGINT

function verify_conda_checksum() {
  local installer_file="$1"

  # Download the Miniconda installer hash list and parse out only the sha256 sums
  local hash_list=$(
    curl -s https://raw.githubusercontent.com/conda/conda-docs/master/docs/source/miniconda_hashes.rst 2> /dev/null \
    | grep -o '[0-9a-f]\{64\}')

  # Check if the downloaded installer has a matching hash
  if grep -q "$(shasum -a 256 "$installer_file" | cut -f 1 -d ' ')" <<< "$hash_list"; then
    success "Miniconda installer checksum verified"
  else
    fail "Miniconda checksum verification failed. Something fishy might be happening."
    rm "$installer_file"
    exit 1
  fi
}

function setup_environment() {
  # Setup tinyRNA environment using our generated lock file
  status "Setting up $env_name environment (this may take a while)..."
  conda create --file $platform_lockfile --name $env_name 2>&1 | tee "env_install.log"
  if ! tr -d \\n < env_install.log | grep -q "Executing transaction: ...working... done"; then
    fail "$env_name environment setup failed"
    echo "Console output has been saved to env_install.log."
    exit 1
  else
    success "$env_name environment setup complete"
  fi
}

function setup_macOS_command_line_tools() {
  # Install Xcode command line tools if necessary
  if ! xcode-select --print-path > /dev/null 2>&1; then
    status "Installing Xcode command line tools. Follow prompts in new window..."
    if xcode-select --install; then
      success "Command line tools setup complete"
    else
      fail "Command line tools installation failed"
      exit 1
    fi
  else
    success "Xcode command line tools are already installed"
  fi
}

################################################################################
# Main Routine
################################################################################

if [[ $CONDA_DEFAULT_ENV == "$env_name" ]]; then
    fail "You must deactivate the $env_name environment before running this script."
    exit 1
fi

# Check if os is mac or linux
if [[ "$OSTYPE" == "darwin"* ]]; then
  success "macOS detected"
  arch=$(uname -m)  # Support Apple Silicon
  shell=$(basename "$(dscl . -read ~/ UserShell | cut -f 2 -d " ")")
  miniconda_installer="Miniconda3-py${miniconda_python_version}_${miniconda_version}-MacOSX-${arch}.sh"
  platform_lockfile="${cwd}/conda/conda-osx-64.lock"
  setup_macOS_command_line_tools
elif [[ "$OSTYPE" == "linux-gnu" ]]; then
  success "Linux detected"
  shell="$(basename "$SHELL")"
  miniconda_installer="Miniconda3-py${miniconda_python_version}_${miniconda_version}-Linux-x86_64.sh"
  platform_lockfile="${cwd}/conda/conda-linux-64.lock"
else
  fail "Unsupported OS"
  exit 1
fi

# Check if Conda is installed
if grep -q "conda init" ~/."$shell"rc ~/."$shell"_profile 2> /dev/null; then
  success "Conda is already installed for $shell"
  eval "$(conda shell."$shell" hook)"
  miniconda_installed=0
else
  status "Downloading Miniconda..."
  curl -O -# https://repo.anaconda.com/miniconda/$miniconda_installer
  if [ -f $miniconda_installer ]; then
    success "Miniconda downloaded"
    verify_conda_checksum $miniconda_installer
    status "Running interactive Miniconda installer..."
    # Use bash since the installer appears to no longer work with zsh
    if ! bash $miniconda_installer; then
      fail "Miniconda installation failed"
      exit 1
    fi
  else
    fail "Miniconda failed to download"
    exit 1
  fi

  # Finalize installation
  source ~/."$shell"rc
  eval "$(conda shell."$shell" hook)"
  conda config --set auto_activate_base false

  success "Miniconda installed"
  miniconda_installed=1
  rm $miniconda_installer
fi

# Check if the conda environment $env_name exists
if conda env list | grep -q "^${env_name}\s"; then
  echo
  echo "The Conda environment $env_name already exists."
  echo "It must be removed and recreated."
  echo
  read -p "Would you like to proceed? [y/n]: " -n 1 -r

  if [[ $REPLY =~ ^y$ ]]; then
    echo
    echo
    status "Removing $env_name environment..."
    conda env remove -n "$env_name" -y > /dev/null 2>&1
    success "Environment removed"
    setup_environment
  elif [[ $REPLY =~ ^n$ ]]; then
    echo
    echo
    fail "Exiting..."
    exit 1
  else
    echo
    fail "Invalid option: $REPLY"
    exit 1
  fi
else
  # Environment doesn't already exist. Create it.
  setup_environment
fi

# Activate environment and set environment variable config for Linux stability
conda activate $env_name
conda env config vars set PYTHONNOUSERSITE=1 > /dev/null  # FYI: cannot be set by lockfile

# Install the tinyRNA codebase
status "Installing tinyRNA codebase via pip..."
if ! pip install "$cwd" > "pip_install.log" 2>&1; then
  fail "Failed to install tinyRNA codebase"
  echo "Check the pip_install.log file for more information."
  exit 1
fi
success "tinyRNA codebase installed"

success "Setup complete"
if [[ $miniconda_installed -eq 1 ]]; then
  status "First, run this one-time command to finalize the Miniconda installation:"
  echo
  echo "  source ~/${shell}.rc"
  echo
fi
status "To activate the environment, run:"
echo
echo "  conda activate $env_name"
echo
