#!/usr/bin/env bash

[[ $# -eq 1 ]] && env_name=$1 || env_name="tinyrna"
bioc_version="3.14"
tested_bioc_versions="3.1[2-4]"

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
  conda create --file $platform_lock_file --name $env_name 2>&1 | tee "env_install.log"
  if ! grep -q "Executing transaction: ...working... done" env_install.log; then
    fail "$env_name environment setup failed"
    echo "Console output has been saved to env_install.log."
    exit 1
  else
    success "$env_name environment setup complete"
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
  shell=$(basename "$(dscl . -read ~/ UserShell | cut -f 2 -d " ")")
  miniconda_installer="Miniconda3-latest-MacOSX-x86_64.sh"
  platform_lock_file="./conda/conda-r-osx-64.lock"
elif [[ "$OSTYPE" == "linux-gnu" ]]; then
  success "Linux detected"
  shell="$(basename "$SHELL")"
  miniconda_installer="Miniconda3-latest-Linux-x86_64.sh"
  platform_lock_file="./conda/conda-r-linux-64.lock"
else
  fail "Unsupported OS"
  exit 1
fi

# Check if Conda is installed
if grep -q "conda init" ~/."$shell"rc ~/."$shell"_profile 2> /dev/null; then
  success "Conda is already installed for $shell"
  eval "$(conda shell."$shell" hook)"
else
  status "Downloading Miniconda..."
  curl -O -# https://repo.anaconda.com/miniconda/$miniconda_installer
  if [ -f $miniconda_installer ]; then
    success "Miniconda downloaded"
    verify_conda_checksum $miniconda_installer
    status "Running interactive Miniconda installer..."
    if ! $shell $miniconda_installer; then
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
  rm $miniconda_installer
fi

# By default, assume that host environment does not contain an appropriate DESeq2 version
install_R_deseq2=true

# Check if R is installed
if [ -x "$(command -v R)" ]; then
  status "Checking host R environment..."
  # Check if DESeq2 is installed
  if Rscript -e "library(DESeq2); print(TRUE)" 2>&1 | tail -n 1 | grep -q TRUE; then
    # Get installed Bioconductor version
    host_bioc_vers=$(Rscript -e "library(BiocManager); BiocManager::version()" 2>&1 | tail -n 1 | grep -Eo '[0-9]+\.[0-9]+')
    # Check to see if host_bioc_version is in our tested range
    if [[ $host_bioc_vers =~ $tested_bioc_versions ]]; then
      success "DESeq2 is already installed in the host environment"
      install_R_deseq2=false
    else
      echo
      echo "tinyRNA has been tested with DESeq2 in Bioconductor release $tested_bioc_versions." \
      "The installer found v$host_bioc_vers on your system. tinyRNA can use your copy or we can" \
      "install a tested version of DESeq2 and R in the isolated tinyRNA environment." | fold -s
      echo
      echo "BEWARE: installation of DESeq2 will take over 20 minutes."
      echo
      read -p "Would you like tinyRNA to use your copy of DESeq2? [y/n]: " -n 1 -r
      echo
      if [[ $REPLY =~ ^[Yy]$ ]]; then
        success "The host's DESeq2 installation will be used"
        install_R_deseq2=false
      elif [[ $REPLY =~ ^[^YyNn]$ ]]; then
        fail "Invalid option: $REPLY"
        exit 1
      fi
    fi # End of Bioconductor version check
  fi # End DESeq2 check
fi # End of R check

if [[ $install_R_deseq2 == false ]]; then
  # Switch to using non-R lock file
  platform_lock_file="${platform_lock_file//-r/}"
fi

# Check if the conda environment $env_name exists
if conda env list | grep -q "$env_name"; then
  echo
  echo "The Conda environment $env_name already exists."
  echo "    1) Update environment"
  echo "    2) Remove and recreate environment"
  echo
  read -p "Select an option [1/2]: " -n 1 -r

  if [[ $REPLY =~ ^[1]$ ]]; then
    echo
    echo
    status "Updating $env_name environment..."
    conda update --file $platform_lock_file --name "$env_name" 2>&1 | tee "env_update.log"
    if ! grep -q "Executing transaction: ...working... done" env_update.log; then
      fail "Failed to update the environment"
      echo "Check the env_update.log file for more information."
      exit 1
    fi
    success "$env_name environment updated"
  elif [[ $REPLY =~ ^[2]$ ]]; then
    echo
    echo
    status "Removing $env_name environment..."
    conda env remove -n "$env_name" -y > /dev/null 2>&1
    success "Environment removed"
    setup_environment
  else
    echo
    fail "Invalid option: $REPLY"
    exit 1
  fi
else
  # Environment doesn't already exist. Create it.
  setup_environment
fi

# Activate tinyRNA environment
conda activate $env_name

# Install pip dependencies and our codebase
status "Installing pip dependencies..."
if ! pip --use-feature=in-tree-build install htseq==0.13.5 . > "pip_install.log" 2>&1; then
  fail "Failed to install pip dependencies"
  echo "Check the pip_install.log file for more information."
  exit 1
fi
success "pip dependencies installed"

if [[ $install_R_deseq2 == true ]]; then
  # Install DESeq2 from Bioconductor
  status "Installing DESeq2 from Bioconductor (this may take over 20 minutes)..."
  status 'To check status run "tail -f deseq2_install.log" from another terminal'
  Rscript -e "install.packages(\"BiocManager\", version=\"$bioc_version\", repos=\"https://cloud.r-project.org\")" > "deseq2_install.log" 2>&1
  Rscript -e "BiocManager::install(\"DESeq2\", version=\"$bioc_version\")" >> "deseq2_install.log" 2>&1

  # Check if DESeq2 installation was successful
  if grep -q "DONE (DESeq2)" "deseq2_install.log"; then
    success "DESeq2 installation was successful"
  else
    fail "DESeq2 installation failed"
    echo "See deseq2_install.log for more information"
    exit 1
  fi
fi

success "Setup complete"
status "To activate the environment, run:"
echo
echo "  conda activate $env_name"
echo
