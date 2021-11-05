#!/usr/bin/env bash

# Get the current shell
shell="$(basename "$SHELL")"

function success() {
  check="✓"
  green_on="\033[1;32m"
  green_off="\033[0m"
  printf "${green_on}${check} %s${green_off}\n" "$*"
}

function fail() {
  nope="⃠"
  red_on="\033[1;31m"
  red_off="\033[0m"
  printf "${red_on}${nope} %s${red_off}\n" "$*"
}

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

function verify_xquartz_checksum() {
  installer_file="$1"

  # XQuartz doesn't list their checksums so we'll have to check against a hardcoded hash
  echo "0eb477f3f8f3795738df9a6a8a15e3858fe21fe9d7103b8d13861bb3aa509e3b  $installer_file" | shasum -a 256 -c -
  if [ $? -eq 0 ]; then
    success "XQuartz installer checksum verified"
  else
    fail "XQuartz checksum verification failed. Something fishy might be happening."
    rm "$installer_file"
    exit 1
  fi
}

function verify_r_checksum() {
  installer_file="$1"

  # Get SHA-1 checksum from R website
  hash_list=$(curl -s https://cran.r-project.org/bin/macosx/ | grep -o '[0-9a-f]\{40\}')

  # See if the downloaded installer has a matching hash
  if shasum -a 1 "$installer_file" | grep -q "$hash_list"; then
    success "R installer checksum verified"
  else
    fail "R checksum verification failed. Something fishy might be happening."
    rm "$installer_file"
    exit 1
  fi
}

function prompt_sudo() {
  if [ "$shell" == "zsh" ]; then
    sudo -n true 2> /dev/null || (echo "Please enter your password to install Miniconda and XQuartz" && sudo -v)
  else
    sudo -v || exit 1
  fi
}

# check if os is mac or linux
if [[ "$OSTYPE" == "darwin"* ]]; then
  # Mac OSX
  success "OSX detected"

  # Check if user has root privileges
  if [[ $EUID -ne 0 ]]; then
    echo
    echo "XQuartz and R require root privileges for installation."
    echo "Would you like to provide your password? This is optional."
    echo "If you answer no, you'll have to run their installers manually."
    echo

    if [ "$shell" == "zsh" ]; then
      read -r -k 1 "Do you want to provide your sudo login? [y/N] "
    else
      read -r -n 1 -p "Do you want to provide your sudo login? [Y/n] "
    fi

    if [[ $REPLY =~ ^[Yy] ]]; then
      sudo -v
      download_only=false
    else
      download_only=true
    fi

  # Check if Conda is installed
  if grep -q "conda init" ~/.${shell}rc; then
    success "Conda already installed for ${shell}"
  else
    echo "Downloading Miniconda..."
    curl -O -# https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    if [ -f "Miniconda3-latest-MacOSX-x86_64.sh" ]; then
      success "Miniconda downloaded"
      verify_conda_checksum "Miniconda3-latest-MacOSX-x86_64.sh"
      $SHELL Miniconda3-latest-MacOSX-x86_64.sh
      if [ $? -eq 1 ]; then
        fail "Miniconda installation failed"
        exit 1
      fi
    else
      fail "Miniconda failed to download"
      exit 1
    fi

    # Cleanup
    success "Miniconda installed"
    rm Miniconda3-latest-MacOSX-x86_64.sh
    rm miniconda_hashes.rst
    rm Miniconda_install.log
  fi

  # Check if XQuartz is installed
  if [ -x "$(command -v xquartz)" ]; then
    success "XQuartz is already installed"
  else
    echo "Downloading XQuartz..."
    curl -O -# https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.1/XQuartz-2.8.1.dmg
    if [ ! -f "XQuartz-2.8.1.dmg" ]; then
      fail "XQuartz failed to download"
      exit 1
    fi

    success "XQuartz downloaded"
    verify_xquartz_checksum "XQuartz-2.8.1.dmg"

    if ! "$download_only"; then
        sudo hdiutil attach XQuartz-2.8.1.dmg
        sudo installer -verbose -pkg /Volumes/XQuartz-2.8.1/XQuartz.pkg -target / 2>&1 "XQuartz_install.log"

        # installer return code is 1 if it fails
        if [ $? -eq 1 ]; then
          fail "XQuartz installation failed"
          echo "See XQuartz_install.log for more information"
          exit 1
        fi

        # Cleanup
        success "XQuartz installation was successful"
        sudo hdiutil detach /Volumes/XQuartz-2.8.1
        rm XQuartz-2.8.1.dmg
        rm XQuartz_install.log
      fi
    fi
  fi

  # Check if R is installed
  if [ -x "$(command -v R)" ]; then
    success "R is already installed"
  else
    echo "Downloading R..."
    curl -O -# https://cran.r-project.org/bin/macosx/base/R-4.1.2.pkg
    if [ ! -f "R-4.1.2.pkg" ]; then
      fail "R failed to download"
      exit 1
    fi

    success "R downloaded"
    verify_r_checksum "R-4.1.2.pkg"

    if ! "$download_only"; then
      sudo installer -verbose -pkg R-4.1.2.pkg -target / 2>&1 "R_install.log"

      # installer return code is 1 if it fails
      if [ $? -eq 1 ]; then
        fail "R installation failed"
        echo "See R_install.log for more information"
        exit 1
      fi

      # Cleanup
      success "R installation was successful"
      rm R-4.1.2.pkg
      rm R_install.log
    fi
  fi

  # Enable use of Conda through this script
  eval "$(conda shell."$shell" hook)"

  # Install pip dependencies and our setup.py
  pip install htseq==0.13.5
  pip --use-feature=in-tree-build install .

  # Setup tinyRNA environment using our generated lock file
  conda create -n tinyRNA --file conda-osx-64.lock


elif [[ "$OSTYPE" == "linux-gnu" ]]; then
  # linux
  echo "Linux"
else
  echo "Your operating system was not recognized"

fi

