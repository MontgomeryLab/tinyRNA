#!/usr/bin/env bash

# USAGE: if you would like to specify a different environment file,
# you can pass the preferred filename as the first argument to the script.
# Otherwise, environment.yml is used.
#   ./update-lockfiles.sh preferred_file

[[ $# -eq 1 ]] && env_file=$1 || env_file="environment.yml"

dev_env_name="tinyrna_dev_v1.0"
dev_deps="conda-forge::conda-lock python=3.9"

#==== FUNCTION DEFINITIONS ================================================================

function check_lockfile_update_status() {
  if [[ $1 -eq 4 ]]; then
    # Supposedly reported with --check-input-hash
    # However this appears to not apply to subdependencies
    success "No updates for lockfile"
  elif [[ $1 -eq 0 ]]; then
    success "Lockfile has been updated"
  else
    fail "Error occurred while updating lockfile. Please see the associated log file."
    exit 1
  fi
}

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

#==== MAIN ROUTINE ========================================================================

# Check if os is mac or linux
if [[ "$OSTYPE" == "darwin"* ]]; then
  success "macOS detected"
  shell=$(basename "$(dscl . -read ~/ UserShell | cut -f 2 -d " ")")
elif [[ "$OSTYPE" == "linux-gnu" ]]; then
  success "Linux detected"
  shell="$(basename "$SHELL")"
else
  fail "Unsupported OS"
  exit 1
fi

# Enable the use of Conda commands in this script
eval "$(conda shell."$shell" hook)"

# Create the dev environment if needed
if ! conda env list | grep -q ${dev_env_name}; then
  status "Creating ${dev_env_name} environment..."
  conda create --name "${dev_env_name}" --yes "${dev_deps}" 2>&1 | tee "dev_env_create.log"
fi

# Activate dev environment
conda activate "${dev_env_name}"

# Update Linux lockfile
status "Updating linux-64 lockfile..."
conda-lock lock --kind explicit --platform linux-64 --file "$env_file" > "linux_lockfile.log" 2>&1

# Check result for Linux lockfile
check_lockfile_update_status $?

# Update macOS lockfile
status "Updating osx-64 lockfile..."
conda-lock lock --kind explicit --platform osx-64 --file "$env_file" > "osx_lockfile.log" 2>&1

# Check result for OSX lockfile
check_lockfile_update_status $?