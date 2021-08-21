#!/usr/bin/env python
import os
import sys
import setuptools

from setuptools.command.install import install

# Package metadata
NAME = 'tinyrna'
DESCRIPTION = 'Automated Quantitative Analysis of Transcript Expression'
URL = 'https://github.com/MontgomeryLab/tinyrna/'
EMAIL = 'ajtate@colostate.edu'
AUTHOR = 'Kristen Brown, Alex Tate'
REQUIRES_PYTHON = '>=3.7.0'
VERSION = '0.1'

# Required packages
REQUIRED = []


class PreFlight(install):
    def run(self):
        if not all([os.getenv(conda_var) for conda_var in ["CONDA_PREFIX", "CONDA_DEFAULT_ENV"]]):
            sys.exit("CRITICAL ERROR: you appear to be installing %s outside of a conda environment.\n"
                     "Instead, please run: conda env create -f environment.yml" % (NAME,))
        else:
            install.run(self)


setuptools.setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=EMAIL,
    description=DESCRIPTION,
    cmdclass={'install': PreFlight},
    packages=setuptools.find_packages(exclude=['tests/*']),
    include_package_data=True,
    package_data={'tinyrna': ['cwl/tools/*.cwl',
                             'cwl/workflows/*.cwl', 
                             'extras/*']},
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'tinyrna = tinyrna.entry:main',
            'tinyrna-config = tinyrna.srna.Configuration:Configuration.main',
            'tinyrna-collapse = tinyrna.srna.collapser:main',
            'tinyrna-count = tinyrna.srna.counter.counter:main',
            'tinyrna-plot = tinyrna.srna.plotter:main'
        ]
    },
    scripts=['tinyrna/srna/tinyrna-deseq.r'],
    python_requires=REQUIRES_PYTHON,
    install_requires=REQUIRED,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
)
