#!/usr/bin/env python
import os
import sys
import pysam
import setuptools

from setuptools.command.install import install
from Cython.Build import cythonize

# Package metadata
NAME = 'tinyrna'
DESCRIPTION = 'Precision analysis of small RNA high-throughput sequencing data'
URL = 'https://github.com/MontgomeryLab/tinyrna/'
EMAIL = 'ajtate@colostate.edu'
AUTHOR = 'Kristen Brown, Alex Tate'
PLATFORM = 'Unix'
REQUIRES_PYTHON = '>=3.9.0'
VERSION = '1.4.0'
REQUIRED = []  # Required packages are installed via Conda's environment.yml


class PrereqAndExec(install):
    """These checks are performed prior to installation to ensure that this setup routine
    is only executed within a Conda environment. Users should perform installation via
    setup.sh, not setup.py."""

    def run(self):
        if not self.in_conda_env():
            sys.exit("CRITICAL ERROR: you appear to be installing %s outside of a conda environment.\n" % (NAME,) +
                     "Instead, please run: ./setup.sh" )
        else:
            install.run(self)

    def in_conda_env(self):
        return all([os.getenv(conda_var) for conda_var in ["CONDA_PREFIX", "CONDA_DEFAULT_ENV"]])


def get_cython_extension_defs():
    """Returns a list of Extension objects corresponding to the contents of pyx_files.
    Extensions indicated as optional in pyx_files will NOT cause the installation to
    error out if there are build issues, and therefore must be used as optional imports."""

    cxx_extension_args = {
        'extra_compile_args': ['-std=c++11', '-O3'],
        'extra_link_args': [],
        'language': 'c++'
    }

    if sys.platform == "darwin" and os.getenv('SDKROOT') is None:
        # If building with the environment variables set by Conda, this isn't necessary
        # However build shortcuts in some IDEs will strip these vars
        sdk_root = get_macos_sdk_path()
        cxx_extension_args['extra_compile_args'] += ['-isystem', os.path.join(sdk_root, '/usr/include')]
        cxx_extension_args['extra_link_args'] += ['-L' + os.path.join(sdk_root, '/usr/lib')]

    StepVector = setuptools.Extension(
        "tiny.rna.counter.stepvector._stepvector",
        sources=["tiny/rna/counter/stepvector/_stepvector.pyx"],
        optional=True,
        **cxx_extension_args
    )

    StepVector_test = setuptools.Extension(
        "tests.cython_tests.stepvector.test_cython",
        sources=["tests/cython_tests/stepvector/test_cython.pyx"],
        optional=True,
        **cxx_extension_args
    )

    AlignmentIter = setuptools.Extension(
        "tiny.rna.counter.parsing.alignments",
        sources=["tiny/rna/counter/parsing/alignments.pyx"],
        extra_link_args=pysam.get_libraries(),
        include_dirs=pysam.get_include(),
        define_macros=pysam.get_defines()
    )

    return [StepVector, StepVector_test, AlignmentIter]


def get_macos_sdk_path():
    """Determines the SDK path for compiler dependencies in a manner that will be
    compatible with conda-build pipelines. The following code was copied from
    https://github.com/python-pillow/Pillow/blob/main/setup.py """

    try:
        import subprocess
        sdk_path = subprocess.check_output(["xcrun", "--show-sdk-path"]).strip().decode('latin1')
    except Exception:
        sdk_path = None

    xcode_sdk_path = "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk"
    commandlinetools_sdk_path = "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk"

    if not sdk_path or sdk_path == xcode_sdk_path:
        if os.path.exists(commandlinetools_sdk_path):
            sdk_path = commandlinetools_sdk_path
        else:
            raise RuntimeError("The macOS SDK path could not be found. This is required for Cython. "
                               "Please run xcode-select --install in your terminal.")

    return sdk_path


if os.getenv('tiny_count_only') == '1':
    print("Only tiny-count will be installed because the environment variable \n"
          '"tiny_count_only" is set to 1', file=sys.stderr)

    console_scripts = ['tiny-count = tiny.rna.counter.counter:main']
    scripts = []
else:
    console_scripts = [
        'tiny = tiny.entry:main',
        'tiny-config = tiny.rna.configuration:Configuration.main',
        'tiny-collapse = tiny.rna.collapser:main',
        'tiny-count = tiny.rna.counter.counter:main',
        'tiny-plot = tiny.rna.plotter:main'
    ]
    scripts = ['tiny/rna/tiny-deseq.r']

setuptools.setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=EMAIL,
    description=DESCRIPTION,
    cmdclass={'install': PrereqAndExec},
    include_package_data=True,
    packages=['tiny'],
    zip_safe=False,
    entry_points={'console_scripts': console_scripts},
    ext_modules=cythonize(
        get_cython_extension_defs(),
        compiler_directives={'language_level': '3'},
        include_path=pysam.get_include(),
        gdb_debug=False
    ),
    scripts=scripts,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: Unix',
    ],
)
