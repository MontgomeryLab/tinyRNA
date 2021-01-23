#!/bin/bash
eval "$(conda shell.bash hook)"

# This script creates a conda environment for every version of Python
# that is listed in versions array. It then installs aquatx in each
# environment and reports if there was an error during installation.

versions=( 3.7 3.8 3.9 )

for v in "${versions[@]}"; do
  echo "Testing v$v"
  # Create/overwrite environment
  conda create -n aquatx-srna-$v python=$v --yes > /dev/null || exit

  # Install aquatx. Read conda output for error strings in real-time
  conda env update -f ../environment.yml -n aquatx-srna-$v 2>&1 | while read line; do
    if [[ $line == "Found conflicts!"* || $line == "UnsatisfiableError"* ]]; then
      echo "Failed env update."
      # Send SIGINT (CTRL+C) to conda for graceful exit from env update
      # This avoids corrupting test environments which in some cases would break this script
      pid=$(ps -e -o pid,command | grep 'environment.yml -n aquatx-srna' | awk 'NR==1{print $1}')
      kill -INT $pid
    fi
  done

  # Activate the environment to verify intended Python version
  conda activate aquatx-srna-$v
  if [[ $(python -V 2>&1) != "Python $v"* ]]; then echo "Failed version check."; fi
  conda deactivate
done