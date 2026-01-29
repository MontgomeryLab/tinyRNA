#!/usr/bin/env bash

# USAGE: if you would like to install under a different conda environment name,
# you can pass the preferred name as the first argument to the script:
#   ./setup.sh preferred_name

env_name=${1:-tinyrna}
miniconda_version="25.1.1-2"
cwd="$(cd "$(dirname "$0")" && pwd -P)"
export ts=$(date +%Y-%m-%d_%H-%M-%S) && readonly ts
# Ensure common Conda locations are discoverable before detection/bootstrapping
export PATH="/opt/homebrew/bin:/usr/local/bin:${HOME}/.local/bin:${HOME}/bin:${HOME}/miniconda3/bin:${HOME}/anaconda3/bin:${HOME}/miniforge3/bin:${HOME}/mambaforge/bin:${HOME}/opt/miniconda3/bin:${HOME}/opt/anaconda3/bin:${PATH}"
# This is the default Python version that will be used by Miniconda (if installation of Miniconda is required).
# Note that this isn't the same as the tinyRNA environment's Python version.
# The tinyRNA environment's Python version is instead specified in the platform lockfile.
miniconda_python_version="310"


######------------------------------ HELPER FUNCTIONS -------------------------------######


# Ensure `conda activate` works reliably inside this non-interactive script.
function conda_bootstrap() {
  local conda_exe conda_root conda_sh

  # micromamba doesn't use conda.sh like conda/mamba; use the hook
  if [[ "${CONDA:-conda}" == "micromamba" ]]; then
    conda_exe="${MICROMAMBA_EXE:-$(type -P micromamba 2>/dev/null)}"
    [[ -n "$conda_exe" && -x "$conda_exe" ]] || return 1
    eval "$("$conda_exe" shell hook -s bash)" || return 1
    return 0
  fi

  # Prefer conda for locating base
  conda_exe="$(type -P conda 2>/dev/null)" || conda_exe="$(type -P "${CONDA:-conda}" 2>/dev/null)" || return 1
  conda_root="$("$conda_exe" info --base 2>/dev/null)" || return 1
  conda_sh="${conda_root}/etc/profile.d/conda.sh"

  if [[ -f "$conda_sh" ]]; then
    # shellcheck disable=SC1090
    source "$conda_sh"
    return 0
  fi

  return 1
}

# Configure build environment for pip-installed native extensions on macOS.
# Forces use of the system Apple clang toolchain and an older deployment target
# to avoid compilation failures (e.g. _Float16 errors) caused by Conda compiler
# wrappers and newer macOS SDK headers. This affects build-time behavior only
# and does not modify the Conda environment or lockfile contents.
function configure_build_env() {
  status "Configuring build environment"

  if [[ "$OSTYPE" == "darwin"* ]]; then
    # Force an x86_64-safe macOS SDK surface.
    # Prevents _Float16 errors when building pip extensions on newer macOS SDKs.
    export MACOSX_DEPLOYMENT_TARGET=10.15
    export SDKROOT="$(xcrun --sdk macosx --show-sdk-path)"

    export CC="/usr/bin/clang"
    export CXX="/usr/bin/clang++"

    export CFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}"
    export CXXFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}"
    export LDFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}"

    status "macOS build target set to ${MACOSX_DEPLOYMENT_TARGET}"
  fi
}

function success() {
  local check="✓"
  local green_on="\033[1;32m"
  local green_off="\033[0m"
  printf "${green_on}${check} %s${green_off}\n" "$*"
}

function status() {
  local blue_on="\033[1;34m"
  local blue_off="\033[0m"
  printf "${blue_on}%s${blue_off}\n" "$*"
}

function warn() {
  local exclaim="⚠"
  local yellow_on="\033[1;33m"
  local yellow_off="\033[0m"
  printf "${yellow_on}${exclaim} %s${yellow_off}\n" "$*"
}

function fail() {
  local nope="⃠"
  local red_on="\033[1;31m"
  local red_off="\033[0m"
  printf "${red_on}${nope} %s${red_off}\n" "$*"
}

function stop() {
  kill -TERM -$$
}
# Ensures that when user presses Ctrl+C, the script stops
# rather than stopping current task and proceeding to the next
trap 'stop' SIGINT

function ensure_conda_on_path() {
  local target_user target_home
  target_user="${SUDO_USER:-$USER}"

  if command -v getent >/dev/null 2>&1; then
    target_home="$(getent passwd "$target_user" | cut -d: -f6)"
  else
    target_home="$HOME"
  fi
  export CONDA_TARGET_HOME="$target_home"

  # If conda is already visible, do nothing
  if command -v conda >/dev/null 2>&1 || command -v mamba >/dev/null 2>&1 || command -v micromamba >/dev/null 2>&1; then
    return 0
  fi

  export PATH="/opt/homebrew/bin:/usr/local/bin:${target_home}/.local/bin:${target_home}/bin:${target_home}/miniconda3/bin:${target_home}/anaconda3/bin:${target_home}/miniforge3/bin:${target_home}/mambaforge/bin:/opt/conda/bin:${PATH}"
}

function get_host_conda_command() {
  if [[ -n "${MICROMAMBA_EXE:-}" && -x "${MICROMAMBA_EXE:-}" ]]; then
    echo "micromamba"
  elif command -v conda > /dev/null 2>&1; then
    echo "conda"
  elif command -v mamba > /dev/null 2>&1; then
    echo "mamba"
  elif command -v micromamba > /dev/null 2>&1; then
    echo "micromamba"
  else
    return 1
  fi
}

function get_shell_rcfile() {
  local shell="$1"
  local os="$2"
  case $shell in
    bash)
      # see https://github.com/conda/conda/pull/11849
      if [[ $os == "macOS" ]]; then
        echo "${HOME}/.bash_profile"
      else
        echo "${HOME}/.bashrc"
      fi;;
    zsh)  # |xonsh|tcsh)
      echo "${HOME}/.${shell}rc";;
    # fish)
    #   echo "${HOME}/.config/fish/config.fish";;
    # nu)
    #   echo "${HOME}/.config/nushell/config.nu";;
    *)
      return 1  # shell isn't supported
  esac
}

function get_shell_hook() {
  local shell_current="$1"
  if [[ "$CONDA" == "conda" ]]; then
    "$CONDA" shell."$shell_current" hook
  elif [[ "$CONDA" == "mamba" || "$CONDA" == "micromamba" ]]; then
    "$CONDA" shell hook -s "$shell_current"
  fi
}

function get_init_block_regex() {
  if [[ "$CONDA" == "conda" ]]; then
    echo '/^# >>> conda initialize >>>/,/^# <<< conda initialize <<</p'
  elif [[ "$CONDA" == "mamba" || "$CONDA" == "micromamba" ]]; then
    echo '/^# >>> mamba initialize >>>/,/^# <<< mamba initialize <<</p'
  fi
}

function download_and_install_miniconda() {
  local miniconda_installer="$1"
  status "Downloading Miniconda..."
  curl -O -# "https://repo.anaconda.com/miniconda/${miniconda_installer}"
  if [ -f "$miniconda_installer" ]; then
    success "Miniconda downloaded"
    if ! verify_miniconda_checksum "$miniconda_installer"; then
      fail "Miniconda checksum verification failed"
      stop
    fi
    status "Running interactive Miniconda installer..."
    # Use bash since the installer appears to no longer work with zsh
    if ! bash "$miniconda_installer"; then
      fail "Miniconda installation failed"
      stop
    fi
  else
    fail "Miniconda download failed"
    stop
  fi
  
  # Ensure freshly-installed Miniconda is discoverable in this script.
  export PATH="${CONDA_TARGET_HOME:-$HOME}/miniconda3/bin:${CONDA_TARGET_HOME:-$HOME}/opt/miniconda3/bin:${PATH}"
  export CONDA=conda
  conda_bootstrap || { fail "Failed to initialize conda (conda.sh not found)"; stop; }
  local new_conda
  new_conda="$(type -P conda 2>/dev/null)" || true
  if [[ -n "$new_conda" ]]; then
    export PATH="$(cd "$(dirname "$new_conda")" && pwd -P):${PATH}"
  fi
}

function verify_miniconda_checksum() {
  local installer_file; local repo_index; local installer_hash; local expected_hash;

  installer_file="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    installer_hash=$(set -o pipefail && sha256sum "$installer_file" | awk '{print $1}') || return 1
  else
    installer_hash=$(set -o pipefail && shasum -a 256 "$installer_file" | awk '{print $1}') || return 1
  fi

  # Get HTML table of all Miniconda versions and their info, including checksums
  if ! repo_index=$(curl -s https://repo.anaconda.com/miniconda/ 2> /dev/null); then
    fail "Failed to download the list of Miniconda installer checksums"
    return 1
  fi

  # Parse installer's expected hash from the table
  expected_hash=$(awk -v target="$installer_file" \
    'BEGIN { FS = "</?td\.*>"; RS = "</?tr>" }
    NF==8 && index($2, target) { print $7; exit; }' \
    <<< "$repo_index")

  if [[ "$installer_hash" == "$expected_hash" ]]; then
    success "Miniconda installer checksum verified"
  elif grep -Fq "$installer_hash" <<< "$repo_index"; then
    # Fallback incase table HTML changes in the future
    success "Miniconda installer checksum verified (fallback)"
  else
    fail "SHA256 checksum for $installer_file"
    fail "Expected: $expected_hash"
    fail "Actual:   $installer_hash"
    rm "$installer_file"
    return 1
  fi
}

function env_prefix() {
  local name="$1"
  local base
  base="$("$CONDA" info --base 2>/dev/null)" || return 1
  echo "${base}/envs/${name}"
}

function env_exists() {
  local name="$1"
  local prefix
  prefix="$(env_prefix "$name" 2>/dev/null)" || prefix=""

  # 1) Fast path: if prefix dir exists, treat as existing (even if conda doesn't list it)
  if [[ -n "$prefix" && -d "$prefix" ]]; then
    return 0
  fi

  # 2) Otherwise fall back to parsing env list (works for conda/mamba/micromamba)
  "$CONDA" env list 2>/dev/null | awk -v target="$name" '
    NF==0 {next}
    $1=="#" {next}
    $1=="Name" {next}
    { env=$1; gsub(/\*/, "", env); if (env==target) found=1 }
    END { exit(found?0:1) }
  '
}

function prune_env_record() {
  local prefix="$1"
  local envs_file="${HOME}/.conda/environments.txt"

  [[ -f "$envs_file" ]] || return 0

  # Remove exact matching line (portable for macOS sed)
  /usr/bin/grep -Fvx "$prefix" "$envs_file" > "${envs_file}.tmp" 2>/dev/null || true
  mv "${envs_file}.tmp" "$envs_file" 2>/dev/null || true
}

function remove_environment() {
  local env_name="$1"
  local logfile="env_remove_${ts}.log"
  local prefix
  prefix="$(env_prefix "$env_name" 2>/dev/null)" || prefix=""

  status "Removing $env_name environment..."

  # Try by name first
  if "$CONDA" env remove -n "$env_name" -y > "$logfile" 2>&1; then
    success "$env_name environment removed"
    return 0
  fi

  # If that fails, try by prefix (helps when conda env list doesn't know about it)
  if [[ -n "$prefix" && -d "$prefix" ]]; then
    if "$CONDA" env remove -p "$prefix" -y >> "$logfile" 2>&1; then
      success "$env_name environment removed (by prefix)"
      return 0
    fi

    # Last resort: directory exists but conda can't remove it -> remove folder
    warn "Conda couldn't remove $env_name, but prefix exists at: $prefix"
    warn "Removing directory directly..."
    rm -rf "$prefix" >> "$logfile" 2>&1 || { fail "Failed to remove $prefix (see ${logfile})"; stop; }
    prune_env_record "$prefix"
    success "$env_name environment directory removed"
    return 0
  fi

  fail "Failed to remove environment (see ${logfile})"
  stop
}

function setup_environment() {
  local env_name="$1"
  local platform_lockfile="$2"
  local logfile="env_install_${ts}.log"

  status "Setting up $env_name environment (this may take a while)..."

  if [[ "$CONDA" == "micromamba" ]]; then
    if [[ -n "${conda_subdir:-}" ]]; then
      if ! CONDA_SUBDIR="$conda_subdir" "$CONDA" create -n "$env_name" --file "$platform_lockfile" -y > "$logfile" 2>&1; then
        CONDA_SUBDIR="$conda_subdir" "$CONDA" create -n "$env_name" -f "$platform_lockfile" -y >> "$logfile" 2>&1 || stop
      fi
    else
      if ! "$CONDA" create -n "$env_name" --file "$platform_lockfile" -y > "$logfile" 2>&1; then
        "$CONDA" create -n "$env_name" -f "$platform_lockfile" -y >> "$logfile" 2>&1 || stop
      fi
    fi
  else
    if [[ -n "${conda_subdir:-}" ]]; then
      CONDA_SUBDIR="$conda_subdir" "$CONDA" create --file "$platform_lockfile" --name "$env_name" -y > "$logfile" 2>&1 || stop
    else
      "$CONDA" create --file "$platform_lockfile" --name "$env_name" -y > "$logfile" 2>&1 || stop
    fi
  fi

  if ! env_exists "$env_name"; then
    fail "$env_name environment setup failed (see ${logfile})"
    stop
  else
    success "$env_name environment setup complete"
  fi
}


function setup_macOS_command_line_tools() {
  # Install Xcode command line tools if necessary
  if xcode-select --print-path > /dev/null 2>&1; then
    success "Xcode command line tools are already installed"
    return 0
  fi

  status "Installing Xcode command line tools. Follow prompts in new window..."
  # This only REQUESTS the install and often returns immediately.
  xcode-select --install >/dev/null 2>&1 || true

  # Wait until xcode-select --print-path succeeds (install completed)
  local i
  for i in {1..120}; do  # up to ~10 minutes (120 * 5s)
    if xcode-select --print-path > /dev/null 2>&1; then
      success "Command line tools setup complete"
      return 0
    fi
    sleep 5
  done

  fail "Command line tools are not installed yet."
  fail "Finish the installation in the pop-up, then re-run this script."
  stop
}


######--------------------------------- PRECHECKS -----------------------------------######


if [[ $CONDA_DEFAULT_ENV == "$env_name" ]]; then
    fail "You must deactivate the $env_name environment before running this script"
    exit 1
fi


######--------------------------------- HOST INFO -----------------------------------######


if [[ "$OSTYPE" == "darwin"* ]]; then
  platform="macOS"
  arch=$(uname -m)
  conda_subdir=""

  if [[ "$arch" == "arm64" ]]; then
    # On Apple Silicon, install Rosetta if it isn't already.
    warn "Apple Silicon detected: creating an osx-64 (x86_64/Rosetta) environment to match the lockfile"
    conda_subdir="osx-64"
    if ! /usr/bin/pgrep -q oahd 2>/dev/null && ! /usr/sbin/pkgutil --pkg-info com.apple.pkg.RosettaUpdateAuto >/dev/null 2>&1; then
      status "Rosetta is required for x86_64 tools; requesting admin permission if needed..."
      /usr/bin/sudo -n /usr/sbin/softwareupdate --install-rosetta --agree-to-license >/dev/null 2>&1 \
        || /usr/bin/sudo /usr/sbin/softwareupdate --install-rosetta --agree-to-license >/dev/null 2>&1 \
        || { warn "Rosetta install was not completed (may require admin approval)."; stop; }
    fi
  fi

  shell_preferred=$(dscl . -read "/Users/$USER" UserShell 2>/dev/null | awk '{print $2}' | xargs basename)
  miniconda_installer="Miniconda3-py${miniconda_python_version}_${miniconda_version}-MacOSX-${arch}.sh"
  platform_lockfile="${cwd}/conda/conda-osx-64.lock"
  setup_macOS_command_line_tools

elif [[ "$OSTYPE" == "linux-gnu" ]]; then
  platform="linux"
  shell_preferred="$(basename "$SHELL")"
  miniconda_installer="Miniconda3-py${miniconda_python_version}_${miniconda_version}-Linux-x86_64.sh"
  platform_lockfile="${cwd}/conda/conda-linux-64.lock"

else
  fail "Unsupported OS"
  exit 1
fi

success "$platform detected"

######-------------------------------- SHELL INFO -----------------------------------######

shell_current=$(ps -o comm= $PPID | cut -f 1 -d " ")
shell_current=${shell_current#-}  # remove the leading dash that login shells have
shell_current=$(basename "$shell_current") # Convert /usr/local/bin/bash -> bash

if [[ "$shell_current" != "$shell_preferred" ]]; then
  warn "The current shell is $shell_current but your default is $shell_preferred"
fi

if ! shellrc=$(get_shell_rcfile "$shell_current" "$platform"); then
  fail "The shell \"$shell_current\" is not supported"
  exit 1
fi


######--------------------------- MINICONDA INSTALLATION ----------------------------######

ensure_conda_on_path

if CONDA=$(get_host_conda_command); then
  export CONDA && readonly CONDA
  success "$CONDA is already installed for $shell_current"

  if ! conda_bootstrap; then
    fail "Failed to initialize conda (conda.sh not found in expected locations)"
    stop
  fi

  miniconda_installed=0
else
  warn "Couldn't find an existing Conda/Mamba installation"
  download_and_install_miniconda "$miniconda_installer"
  export CONDA="conda" && readonly CONDA

  if ! conda_bootstrap; then
    fail "Failed to initialize conda (conda.sh not found in expected locations)"
    stop
  fi

  "$CONDA" config --set auto_activate_base false

  success "Miniconda installed"
  miniconda_installed=1
  rm "$miniconda_installer"
fi


######----------------------------- CREATE ENVIRONMENT ------------------------------######


if env_exists "$env_name"; then
  echo
  echo "The Conda environment \"$env_name\" already exists"
  echo "It must be removed and recreated"
  echo
  read -p "Would you like to proceed? [y/n]: " -n 1 -r
  echo

  if [[ $REPLY =~ ^y$ ]]; then
    remove_environment "$env_name"
  elif [[ $REPLY =~ ^n$ ]]; then
    fail "Exiting..."
    exit 1
  else
    fail "Invalid option: $REPLY"
    exit 1
  fi
fi

# Apply macOS build-time compiler overrides *after* activating the Conda environment.
# This ensures pip builds use the system Apple clang instead of Conda's compiler
# wrappers, which can fail with newer macOS SDKs.
setup_environment "$env_name" "$platform_lockfile"

# NOTE: For micromamba, activation works because conda_bootstrap
# runs `micromamba shell hook -s bash`
if [[ "$CONDA" == "micromamba" ]]; then
  micromamba activate "$env_name"
else
  "$CONDA" activate "$env_name"
fi

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  fail "Environment activation failed"
  stop
fi

configure_build_env

# Avoid writing conda internal metadata (conda-meta/state). If you want
# PYTHONNOUSERSITE, set it in the user's shell or via activation hooks.
# export PYTHONNOUSERSITE=1


######---------------------------- tinyRNA INSTALLATION -----------------------------######


status "Installing tinyRNA codebase via pip..."
logfile="pip_install_${ts}.log"

if ! python -m pip install "$cwd" > "$logfile" 2>&1; then
  fail "Failed to install tinyRNA codebase (see ${logfile})"
  exit 1
fi
success "tinyRNA codebase installed"


######---------------------------------- FINALIZE -----------------------------------######


success "Setup complete"
if [[ $miniconda_installed -eq 1 ]]; then
  echo
  echo "First, run this one-time command to finalize the Miniconda installation:"
  echo
  echo "  source $shellrc"
fi
echo
echo "To activate the environment, run:"
echo
echo "  conda activate $env_name"
echo
