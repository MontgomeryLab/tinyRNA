#!/usr/bin/env python

import setuptools 

# Package metadata
NAME = 'aquatx'
DESCRIPTION = 'Automated Quantitative Analysis of Transcript Expression'
URL = 'https://github.com/MontgomeryLab/aquatx-srna/'
EMAIL = 'kristen.brown@colostate.edu'
AUTHOR = 'Kristen Brown'
REQUIRES_PYTHON = '>=3.7.0'
VERSION = '0.1'

# Required packages
REQUIRED = [
    'cwltool',
    'htseq',
    'numpy',
    'pandas',
    'matplotlib',
]

setuptools.setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=EMAIL,
    description=DESCRIPTION,
    packages=setuptools.find_packages(exclude=['tests/*']),
    include_package_data=True,
    package_data={'aquatx': ['cwl/tools/*.cwl', 
                             'cwl/workflows/*.cwl', 
                             'extras/*.csv', 
                             'extras/*.yml']},
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'aquatx = aquatx.aquatx:main',
            'aquatx-config = aquatx.srna.Configuration:Configuration.main',
            'aquatx-collapse = aquatx.srna.collapser.collapser:main',
            'aquatx-count = aquatx.srna.counter:main'
        ]
    },
    scripts=['aquatx/srna/aquatx-deseq.r'],
    python_requires=REQUIRES_PYTHON,
    install_requires=REQUIRED,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
)
