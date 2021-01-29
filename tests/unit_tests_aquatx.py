#!/usr/bin/env python
""" Testing module for aquatx.py """

import unit_test_helpers as helpers
import aquatx.aquatx as aquatx
import contextlib
import subprocess
import unittest
import psutil
import shutil
import json
import io
import os

from pkg_resources import resource_filename

"""

======== PRE-INSTALL TESTS ========
These tests interface directly with the package Python code
located within the source directory.

"""
class test_aquatx_preinstall(unittest.TestCase):
    def setUp(self):
        self.stdout_capture = io.StringIO()
        self.aquatx_cwl_path = '../aquatx/cwl'
        self.aquatx_extras_path = '../aquatx/extras'
        self.config_file = './testdata/run_config_template.yml'
        self.expected_cwl_dir_tree = {
            'cwl': {
                'files': [],
                'tools': {
                    'files': [
                        'aquatx-deseq.cwl', 'bowtie.cwl',
                        'aquatx-collapse.cwl', 'bowtie-build.cwl',
                        'aquatx-count.cwl', 'aquatx-merge.cwl', 'fastp.cwl'
                    ]
                },
                'workflows': {
                    'files': ['aquatx_wf.cwl']
                }
            }
        }


    """
    Testing that aquatx.get_template copies the correct files 
    to the current directory.
    """
    def test_get_template(self):
        template_files = [
            'run_config_template.yml', 'sample_sheet_template.csv',
            'reference_sheet_template.csv'
        ]

        try:
            before_count = len(os.listdir('.'))
            aquatx.get_template(self.aquatx_extras_path)
            self.assertEqual(
                len(os.listdir('.')) - before_count, 3,
                "Abnormal number of template files. Expected 3.")

            for file in template_files:
                self.assertTrue(os.path.isfile(file),
                                "An expected template file wasn't copied.")
                os.remove(file)
        finally:
            for file in template_files: os.remove(file)


    """
    Testing that aquatx.setup_cwl with --config None/none 
    copies workflow files without mentioning the config file
    """
    def test_setup_cwl_noconfig(self):
        no_config = ['None', 'none']
        for config in no_config:
            try:
                with contextlib.redirect_stdout(self.stdout_capture):
                    aquatx.setup_cwl(self.aquatx_cwl_path, config)

                self.assertNotIn(
                    "configuration", self.stdout_capture.getvalue(),
                    "Setup mentioned configfile when None was provided")
                self.assertEqual(helpers.get_dir_tree('./cwl'),
                                  self.expected_cwl_dir_tree,
                                 "The expected local cwl directory tree was not found")
            finally:
                if os.path.isdir('./cwl'): shutil.rmtree('./cwl')


    """
    Testing that aquatx.setup_cwl WITH config file properly
    processes the input configfile, then copies workflow files.
    """
    def test_setup_cwl_withconfig(self):
        try:
            with contextlib.redirect_stdout(self.stdout_capture):
                aquatx.setup_cwl(self.aquatx_cwl_path, self.config_file)

            self.assertEqual(helpers.get_dir_tree('./cwl'),
                              self.expected_cwl_dir_tree)
            self.assertIn("The processed configuration file is located at: ",
                          self.stdout_capture.getvalue(),
                          "Setup failed to mention the location of the processed config file.")
        finally:
            if os.path.isdir('./cwl'): shutil.rmtree('./cwl')


    """
    Test that aquatx.run spawns a cwltool subprocess
    """
    def test_run(self):
        try:
            with contextlib.redirect_stdout(self.stdout_capture):
                aquatx.run(self.aquatx_cwl_path, self.config_file)
                self_pid = psutil.Process(os.getpid())
                self.assertTrue('cwltool' in self_pid.children)
        finally:
            print(self.stdout_capture.getvalue())

"""

======== POST-INSTALL TESTS ========
These tests interface with the installed aquatx package via 
the command line

"""
class test_aquatx_postinstall(unittest.TestCase):
    def setUp(self):
        self.MIN_PYTHON = 3.7
        os.system(f"conda env create -n aquatx-test python={self.MIN_PYTHON} --yes > /dev/null")
        os.system(f"conda env update -n aquatx-text -f ../environment.yml > /dev/null")

    #def test_get_template(self):


if __name__ == '__main__':
    unittest.main()
