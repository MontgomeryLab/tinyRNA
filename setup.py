#!/usr/bin/env python
import os
import sys

from setuptools import setup, Extension
from setuptools.command.install import install

from Cython.Build import cythonize

# Package metadata
NAME = 'tinyrna'
DESCRIPTION = 'Comprehensive analysis of small RNA high-throughput sequencing data'
URL = 'https://github.com/MontgomeryLab/tinyrna/'
EMAIL = 'ajtate@colostate.edu'
AUTHOR = 'Kristen Brown, Alex Tate'
PLATFORM = 'Unix'
REQUIRES_PYTHON = '>=3.7.0'
VERSION = '0.1'

# Required packages are installed via Conda's environment.yml
# See PreFlight below...
REQUIRED = []


class PreFlight(install):
    def run(self):
        if not all([os.getenv(conda_var) for conda_var in ["CONDA_PREFIX", "CONDA_DEFAULT_ENV"]]):
            sys.exit("CRITICAL ERROR: you appear to be installing %s outside of a conda environment.\n"
                     "Instead, please run: ./setup.sh" % (NAME,))
        else:
            install.run(self)


bitpack_dir = "tiny/rna/bitpack/"
short_seq_common_compile_args = [
    '-stdlib=libc++',
    '-std=c++11',
    "-O3",
    '-march=native']

extensions = [
    Extension("tiny.rna.bitpack.short_seq",
              sources=[bitpack_dir + 'short_seq.pyx'],
              extra_compile_args=short_seq_common_compile_args,
              language='c++'),
    Extension("tiny.rna.bitpack.short_seq_128",
              sources=[bitpack_dir + 'short_seq_128.pyx'],
              extra_compile_args=short_seq_common_compile_args,
              language='c++'),
    Extension("tiny.rna.bitpack.short_seq_64",
              sources=[bitpack_dir + 'short_seq_64.pyx'],
              extra_compile_args=short_seq_common_compile_args,
              language='c++'),
    Extension("tiny.rna.bitpack.short_seq_util",
              sources=[bitpack_dir + 'short_seq_util.pyx'],
              extra_compile_args=short_seq_common_compile_args,
              language='c++'),
    Extension("tiny.rna.bitpack.fast_read",
              sources=[bitpack_dir + 'fast_read.pyx'],
              extra_compile_args=short_seq_common_compile_args,
              language='c++'),
    Extension("tiny.rna.bitpack.umi",
              sources=[bitpack_dir + 'umi.pyx'],
              extra_compile_args=short_seq_common_compile_args,
              language='c++'),
]


setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=EMAIL,
    description=DESCRIPTION,
    cmdclass={'install': PreFlight},
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
    ext_modules=cythonize(extensions, compiler_directives={'language_level': '3'}),
    scripts=['tiny/rna/tiny-deseq.r'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: Unix',
    ],
)
