#!/usr/bin/env python
""" Testing module for aquatx.py """

import aquatx.aquatx as aquatx
import contextlib
import unittest
import filecmp
import shutil
import io
import os

from pkg_resources import resource_filename

global min_python_version
min_python_version = 3.7

class test_aquatx_preinstall(unittest.TestCase):

    def setUp(self):
        self.f = io.StringIO()
        self.aquatx_cwl_path = resource_filename('aquatx', 'cwl/')
        self.aquatx_extras_path = resource_filename('aquatx', 'extras/')

    """
    Testing that aquatx.get_template copies the correct files to the current directory.
    """
    def test_get_template(self):
        template_files = ['run_config_template.yml', 'sample_sheet_template.csv', 'reference_sheet_template.csv']
        before_count = len(os.listdir('.'))
        aquatx.get_template(self.aquatx_extras_path)
        self.assertEquals(len(os.listdir('.')) - before_count, 3, "Abnormal number of template files. Expected 3.")

        for file in template_files:
            self.assertTrue(os.path.isfile(file), "An expected template file wasn't copied.")
            os.remove(file)

    """
    Testing that aquatx.setup_cwl with --config None/none copies workflow files without mentioning the config file
    """
    def test_setup_cwl_noconfig(self):
        no_config = ['None', 'none']
        for config in no_config:
            with contextlib.redirect_stdout(self.f):
                aquatx.setup_cwl(self.aquatx_cwl_path, config)
                filecmp.dircmp(self.aquatx_extras_path, './cwl')

            self.assertTrue(os.path.isdir('cwl/tools') and os.path.isdir('cwl/workflows'), "Expected CWL outputs don't exist.")
            self.assertNotIn("configuration", self.f.getvalue(), "Setup mentioned configfile when None was provided")
            self.assertEquals([f for f in os.walk('./cwl/')], [f for f in os.walk(self.aquatx_cwl_path)])
            shutil.rmtree('./cwl')

    """
    Testing that aquatx.setup
    """
    def test_setup_cwl_withconfig(self):
        config_file = './testdata/run_config_template.yml'
        aquatx.setup_cwl(self.aquatx_cwl_path, config_file)
        self.assertTrue(os.path.isdir('cwl/tools') and os.path.isdir('cwl/workflows'))


class test_aquatx_postinstall(unittest.TestCase):

    def setUp(self):
        os.system(f"conda env create -n aquatx-test python={min_python_version} --yes > /dev/null")
        os.system(f"conda env update -n aquatx-text -f ../environment.yml > /dev/null")

    #def test_get_template(self):




if __name__ == '__main__':
    unittest.main()
