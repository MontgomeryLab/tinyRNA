#!/usr/bin/env python

import unit_test_helpers as helpers
import aquatx.aquatx as aquatx
import contextlib
import subprocess
import threading
import unittest
import psutil
import shutil
import time
import io
import os

"""

Contains tests for aquatx.py, both from direct source-level calls
as well as post-install testing of invocation by terminal. Each
test covers both environments.

"""
class test_aquatx_preinstall(unittest.TestCase):
    def setUp(self):
        # For pre-install tests
        self.aquatx_cwl_path = '../aquatx/cwl'
        self.aquatx_extras_path = '../aquatx/extras'

        # For post-install tests
        os.system("pip install ../ > /dev/null")

        # For both pre and post install
        self.config_file = './testdata/run_config_template.yml'
        self.stdout_capture = io.StringIO()
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
    Testing that get-template copies the correct files 
    to the current directory.
    """
    def test_get_template(self):
        functions = [
            "aquatx.get_template(self.aquatx_extras_path)", # The pre-install invocation
            'os.system("aquatx get-template")'              # The post-install command
        ]
        template_files = [
            'run_config_template.yml', 'sample_sheet_template.csv',
            'reference_sheet_template.csv'
        ]

        for fn in functions:
            try:
                before_count = len(os.listdir('.'))
                eval(fn, globals(), locals())

                # Check that exactly 3 files were produced by the command
                self.assertEqual(
                    len(os.listdir('.')) - before_count, 3,
                    f"Abnormal number of template files. Expected 3. Function: {fn}")

                # Check that each expected file was produced
                for file in template_files:
                    self.assertTrue(os.path.isfile(file),
                                    f"An expected template file wasn't copied: {file}, function: {fn}")
                    os.remove(file)
            finally:
                # Remove the local template files if necessary, even if an exception was thrown above
                for file in template_files:
                    if os.path.isfile(file): os.remove(file)


    """
    Testing that setup-cwl with a None/none config file 
    copies workflow files without mentioning a config file
    """
    def test_setup_cwl_noconfig(self):
        functions = [
            "aquatx.setup_cwl(self.aquatx_cwl_path, config)",   # The pre-install invocation
            'os.system(f"aquatx setup-cwl --config {config}")'  # The post-install command
        ]
        no_config = ['None', 'none']
        for fn in functions:
            for config in no_config:
                try:
                    # Execute the given function and capture its stdout stream
                    with contextlib.redirect_stdout(self.stdout_capture):
                        eval(fn, globals(), locals())

                    # Check that the function did not mention the configuration file
                    self.assertNotIn(
                        "configuration", self.stdout_capture.getvalue(),
                        "Setup mentioned configfile when None was provided")

                    # Check (by name and directory structure) that the expected files/folders were produced
                    self.assertEqual(helpers.get_dir_tree('./cwl'),
                                      self.expected_cwl_dir_tree,
                                     "The expected local cwl directory tree was not found")
                finally:
                    # Remove the copied workflow files even if an exception was thrown above
                    if os.path.isdir('./cwl'): shutil.rmtree('./cwl')


    """
    Testing that setup-cwl WITH config file mentions the location of the 
    processed input configfile, then copies workflow files. Correctness
    of processed config file will be checked in the setup_config tests.
    """
    def test_setup_cwl_withconfig(self):
        functions = [
            "aquatx.setup_cwl(self.aquatx_cwl_path, self.config_file)",     # The pre-install invocation
            'os.system(f"aquatx setup-cwl --config {self.config_file}")'    # The post-install command
        ]
        for fn in functions:
            try:
                # Execute the given function and capture its stdout stream
                with contextlib.redirect_stdout(self.stdout_capture):
                    eval(fn, globals(), locals())

                # Check (by name and directory structure) that the expected files/folders were produced
                self.assertEqual(helpers.get_dir_tree('./cwl'),
                                  self.expected_cwl_dir_tree)

                # Check that the function mentioned the config file
                self.assertIn("The processed configuration file is located at: ",
                              self.stdout_capture.getvalue(),
                              "Setup failed to mention the location of the processed config file.")
            finally:
                # Remove the copied workflow files even if an exception was thrown above
                if os.path.isdir('./cwl'): shutil.rmtree('./cwl')


    """
    Test that run invocation produces a cwltool subprocess. Since subprocess.run() 
    is a blocking call, we need to call aquatx.run() in its own thread so we can 
    measure its behavior. Otherwise we would have to wait for the whole pipeline 
    to finish before being able to measure it here, and by that time the relevant 
    subprocesses would have already exited. The post-install command accomplishes
    the same thing but in another process rather than another thread.
    """
    def test_run(self):
        def pre_install_command():
            # Run as daemon, meaning the thread exits once test_run finishes
            run_thread = threading.Thread(target=aquatx.run, args=(self.aquatx_cwl_path, self.config_file), daemon=True)
            return run_thread.start()

        def post_install_command():
            # Simulate running the command from the terminal
            return subprocess.Popen(["aquatx", "run", "--config", self.config_file])

        functions = [pre_install_command, post_install_command]
        for fn in functions:
            fn()
            # Allow some time to pass for the cwltool subprocess to start
            time.sleep(1)
            # Get the names of all descendents of the current process
            subprocs = psutil.Process(os.getpid()).children(recursive=True)
            sub_names = [sub.name() for sub in subprocs]
            self.assertIn('cwltool', sub_names, "The cwltool subprocess does not appear to have started.")

            # Remove child processes so that we have accurate results when testing the next function
            for subproc in reversed(subprocs):
                subproc.terminate()

if __name__ == '__main__':
    unittest.main()
