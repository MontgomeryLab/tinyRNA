#!/usr/bin/env bash

env_name="tinyrna"
bioc_version="3.14"
deseq2_version="1.34.0"

# Get the current shell
shell="$(basename "$SHELL")"

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
# rather than proceeding to the next step
trap 'stop' SIGINT

# check if os is mac or linux
if [[ "$OSTYPE" == "darwin"* ]]; then
  success "macOS detected"
  miniconda_installer="Miniconda3-latest-MacOSX-x86_64.sh"
  platform_lock_file="conda-osx-64.lock"
elif [[ "$OSTYPE" == "linux-gnu" ]]; then
  success "Linux detected"
  miniconda_installer="Miniconda3-latest-Linux-x86_64.sh"
  platform_lock_file="conda-linux-64.lock"
else
  fail "Unsupported OS"
  exit 1
fi

function verify_conda_checksum() {
  installer_file="$1"

  # Download the Miniconda installer hash list and parse out only the sha256 sums
  hash_list=$(
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
  if ! conda create --file $platform_lock_file --name $env_name  > "env_install.log" 2>&1; then
    fail "$env_name environment setup failed"
    echo "Check the env_install.log file for more information."
    exit 1
  else
    success "$env_name environment setup complete"
  fi
}

# Check if Conda is installed
if grep -q "conda init" ~/."${shell}"rc; then
  success "Conda is already installed for ${shell}"
else
  status "Downloading Miniconda..."
  curl -O -# https://repo.anaconda.com/miniconda/$miniconda_installer
  if [ -f $miniconda_installer ]; then
    success "Miniconda downloaded"
    verify_conda_checksum $miniconda_installer
    status "Running interactive Miniconda installer..."
    if ! $SHELL $miniconda_installer; then
      fail "Miniconda installation failed"
      exit 1
    fi
  else
    fail "Miniconda failed to download"
    exit 1
  fi

  # Cleanup
  success "Miniconda installed"
  rm $miniconda_installer
fi

# Enable use of Conda through this script
eval "$(conda shell."$shell" hook)"

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
    if ! conda update --file $platform_lock_file --name "$env_name" > "env_update.log" 2>&1; then
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

# Check if Deseq2 is up to date (NOTE: $Version is an R directive, not a shell variable, hence single quotes)
if ! Rscript -e 'packageDescription("DESeq2")$Version' 2>&1 | grep -q $deseq2_version; then
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
else
  success "DESeq2 $deseq2_version is already installed"
fi

success "Setup complete"
status "To activate the environment, run:"
echo
echo "  conda activate $env_name"
echo