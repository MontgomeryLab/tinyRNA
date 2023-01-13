
"""Addresses the error 'RWrapper terminated, exit code: 134 ... rpath ...'

This issue is specific to macOS. The R Language plugin for Pycharm crashes when the R has been installed
via Conda. The issue is described here: https://youtrack.jetbrains.com/issue/R-1271

The issue is still present in Pycharm 2022.3.1.
This patch will need to be applied every time the plugin is updated."""

import platform
import sys
import os
import re

from glob import glob

if platform.system() != 'Darwin':
    sys.exit("This patch is only for macOS")

# Get PyCharm directory most recently modified (assumed to be the latest version)
home_dir = os.path.expanduser("~")
pycharm_dirs = glob(f"{home_dir}/Library/Application Support/JetBrains/PyCharm*")
latest_pycharm = sorted(pycharm_dirs, key=os.path.getmtime)[-1]

# The replacement R function
patched_fn = """getLDLibraryPath <- function() {
  conda_path <- Sys.getenv("CONDA_EXE", unset = NA)
  if (!is.na(conda_path)) conda_path = dirname(dirname(conda_path))
  r_path <- Sys.getenv("R_HOME", unset = NA)
  if (!is.na(r_path)) r_path = dirname(dirname(dirname(dirname(r_path))))
  
  if (get_os() == "osx" && r_path != conda_path) Sys.getenv("DYLD_FALLBACK_LIBRARY_PATH")
  else if (get_os() == "linux") Sys.getenv("LD_LIBRARY_PATH")
  else ""
}"""

# Apply the patch
target_file = f"{latest_pycharm}/plugins/r-plugin/R/GetEnvVars.R"
with open(target_file, "r+") as f:
    file_contents = f.read()
    patched = re.sub(r"getLDLibraryPath <- function\(\) {.*}", patched_fn, file_contents, flags=re.DOTALL)

    # Overwrite the file
    f.seek(0)
    f.write(patched)
    f.truncate()