#!/usr/bin/env python
import os
import sys
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
VERSION = '1.0'
REQUIRED = []  # Required packages are installed via Conda's environment.yml


class PrereqAndExec(install):
    def run(self):
        if not self.in_conda_env():
            sys.exit("CRITICAL ERROR: you appear to be installing %s outside of a conda environment.\n" % (NAME,) +
                     "Instead, please run: ./setup.sh" )
        else:
            install.run(self)

    def in_conda_env(self):
        return all([os.getenv(conda_var) for conda_var in ["CONDA_PREFIX", "CONDA_DEFAULT_ENV"]])


pyx_files = [
    'tiny/rna/counter/stepvector/_stepvector.pyx',
    'tests/cython_tests/stepvector/test_cython.pyx'
]
cxx_extension_args = {
    'extra_compile_args': ['-stdlib=libc++', '-std=c++11', '-O3'],
    'extra_link_args': [],
    'language': 'c++'
}
if sys.platform == "darwin":
    cxx_extension_args['extra_compile_args'] += ['-isystem', '/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include']
    cxx_extension_args['extra_link_args'] += ['-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib']
cython_extensions = [
    setuptools.Extension(
        pyx_filename.replace('/', '.').rstrip('.pyx'),
        sources=[pyx_filename],
        **cxx_extension_args
    )
    for pyx_filename in pyx_files
]

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
    entry_points={
        'console_scripts': [
            'tiny = tiny.entry:main',
            'tiny-config = tiny.rna.configuration:Configuration.main',
            'tiny-collapse = tiny.rna.collapser:main',
            'tiny-count = tiny.rna.counter.counter:main',
            'tiny-plot = tiny.rna.plotter:main'
        ]
    },
    ext_modules=cythonize(
        cython_extensions,
        compiler_directives={'language_level': '3'},
        gdb_debug=False
    ),
    scripts=['tiny/rna/tiny-deseq.r'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: Unix',
    ],
)
