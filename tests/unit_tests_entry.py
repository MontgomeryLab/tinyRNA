#!/usr/bin/env python

import os
import shutil
import sys
import time
import unittest
from glob import glob

import psutil

import tiny.entry as entry
import unit_test_helpers as helpers

"""

Contains tests for entry.py, both from direct source-level calls
as well as post-install testing of invocation by terminal. Each
test covers both environments.

"""


class test_entry(unittest.TestCase):
    @classmethod
    def setUpClass(self):

        # For pre-install tests
        self.cwl_path = '../tiny/cwl'
        self.templates_path = '../tiny/templates'

        # For both pre and post install
        self.config_file = f'{self.templates_path}/run_config_template.yml'
        self.expected_cwl_dir_tree = {
            'cwl': {
                'files': set(),
                'tools': {
                    'files': {
                        'tiny-deseq.cwl', 'tiny-plot.cwl', 'bowtie.cwl', 'tiny-collapse.cwl',
                        'bowtie-build.cwl', 'tiny-count.cwl', 'fastp.cwl', 'make-subdir.cwl'
                    }
                },
                'workflows': {
                    'files': {
                        'tinyrna_wf.cwl', 'preprocessing.cwl'
                    }
                }
            }
        }

    """
    Testing that get-templates copies the correct files to the current directory.
    """

    def test_get_templates(self):
        test_functions = [
            helpers.LambdaCapture(lambda: entry.get_templates(self.templates_path)),  # The pre-install invocation
            helpers.ShellCapture("tiny get-templates")                                # The post-install command
        ]
        template_files = ['run_config_template.yml', 'samples.csv', 'features.csv',
                          'paths.yml', 'tinyrna-light.mplstyle']

        def dir_entry_ct():
            return len(os.listdir('.'))

        for test_context in test_functions:
            try:
                # Count number of entries in current directory before test
                dir_before_count = dir_entry_ct()
                with test_context as test:
                    test()

                # Check that exactly 5 files were produced by the command
                self.assertEqual(dir_entry_ct() - dir_before_count, len(template_files),
                    f"Abnormal number of template files. Function: {test_context.__class__.__name__}")

                # Check that each expected file was produced
                for file in template_files:
                    self.assertTrue(os.path.isfile(file), f"An expected template file wasn't copied: {file}, function: "
                                                          f"{test_context.__class__.__name__}")
                    os.remove(file)
            finally:
                # Remove the local template files if necessary, even if an exception was thrown above
                for file in template_files:
                    if os.path.isfile(file): os.remove(file)

    """
    Testing that setup-cwl with a None/none config file copies workflow files without mentioning a config file
    """

    def test_setup_cwl_noconfig(self):
        no_config = ['None', 'none']
        for config in no_config:
            test_functions = [
                helpers.LambdaCapture(lambda: entry.setup_cwl(self.cwl_path, config)),
                helpers.ShellCapture(f"tiny setup-cwl --config {config}")
            ]

            for test_context in test_functions:
                try:
                    with test_context as test:
                        test()

                        # Check that the function did not mention the configuration file
                        self.assertNotIn("configuration", test.get_stdout(),
                            f"Setup mentioned configfile when {no_config} was provided, function: "
                            f"{test_context.__class__.__name__}")

                        # Check (by name and directory structure) that the expected files/folders were produced
                        self.assertEqual(helpers.get_dir_tree('./cwl'), self.expected_cwl_dir_tree,
                                         f"The expected local cwl directory tree was not found, function: "
                                         f"{test_context.__class__.__name__}")
                finally:
                    # Remove the copied workflow files even if an exception was thrown above
                    if os.path.isdir('./cwl'): shutil.rmtree('./cwl')

    """
    Testing that setup-cwl WITH config file mentions the location of the processed input configfile, 
    then copies workflow files. Correctness of processed config file will be checked in the setup_config tests.
    """

    def test_setup_cwl_withconfig(self):
        test_functions = [
            helpers.LambdaCapture(lambda: entry.setup_cwl(self.cwl_path, self.config_file)),
            helpers.ShellCapture(f"tiny setup-cwl --config {self.config_file}")
        ]
        for test_context in test_functions:
            # So that we may reference the filename in the finally block below
            config_file_location = ""

            try:
                # Execute the given function and capture its stdout stream
                with test_context as test:
                    test()
                    stdout_result = test.get_stdout()
                    # Get the name of the processed config file
                    config_file_message = stdout_result.splitlines()[1]

                # Check that the function mentioned the config file with a complete name
                config_file_basename = os.path.basename(self.config_file)
                self.assertEqual(config_file_message,
                                 r'The processed configuration file is located at: processed_' + config_file_basename,
                                 f"Setup failed to mention the location of the processed config file. Function: "
                                 f"{test_context.__class__.__name__}")

                # Check that the processed configuration file exists
                config_file_location = config_file_message.split(": ")[1]
                self.assertTrue(os.path.isfile(config_file_location),
                                f"The processed config file does not exist: {config_file_location}. Function: "
                                f"{test_context.__class__.__name__}")

                os.remove(config_file_location)

                # Check (by name and directory structure) that the expected files/folders were produced
                self.assertDictEqual(helpers.get_dir_tree('./cwl'), self.expected_cwl_dir_tree)
            finally:
                # Remove the copied workflow files even if an exception was thrown above
                if os.path.isdir('./cwl'): shutil.rmtree('./cwl')
                # Remove the output config file
                if os.path.isfile(config_file_location): os.remove(config_file_location)

    """
    Test that run invocation produces a cwltool subprocess. Since subprocess.run() 
    is a blocking call, we need to call tiny.run() in its own thread so we can 
    measure its behavior. Otherwise we would have to wait for the whole pipeline 
    to finish before being able to measure it here, and by that time the relevant 
    subprocesses would have already exited. The post-install command accomplishes
    the same thing but in another process rather than another thread.
    """

    def test_run(self):
        # Non-blocking test functions (invocations continue to run in background until test_context is left)
        test_functions = [
            helpers.LambdaCapture(lambda: entry.run(self.cwl_path, self.config_file), blocking=False),
            helpers.ShellCapture(f"tiny run --config {self.config_file}", blocking=False)
        ]

        def get_children():
            return psutil.Process(os.getpid()).children(recursive=True)

        for test_context in test_functions:
            try:
                with test_context as test:
                    test()

                    # Check for cwltool in child processes up to 5 times, waiting 1 second in between
                    for i in range(10):
                        time.sleep(1)
                        sub_names = {sub.name() for sub in get_children()}
                        if 'node' in sub_names:
                            break

                    self.assertIn('node', sub_names,
                                  f"The cwltool subprocess does not appear to have started. Function: "
                                  f"{test_context.__class__.__name__}")
            except:
                print("Captured stdout:\n" + f'"{test.get_stdout()}"')
                print('\n\n')
                print("Captured stderr:\n" + f'"{test.get_stderr()}"')
            finally:
                run_dirs = glob("./testdata/entry_test_*_run_directory")
                for dir in run_dirs:
                    shutil.rmtree(dir)

    """
    A very minimal test for the subprocess context manager that is used
    to execute post-install tiny commands via a shell.
    """

    def test_ShellCapture_helper(self):
        # Test blocking capture with stdout (though we can still check stderr with a blocking capture...)
        with helpers.ShellCapture('echo "Today is the day"') as fn:
            # Pre-execution test
            self.assertFalse(fn.is_complete())
            self.assertEqual(fn.get_stdout(), '')
            self.assertEqual(fn.get_stderr(), '')
            self.assertEqual(fn.get_exit_status(), None)

            fn()
            self.assertTrue(fn.is_complete())
            self.assertEqual(fn.get_stdout(), "Today is the day\n")
            self.assertEqual(fn.get_stderr(), '')
            self.assertEqual(fn.get_exit_status(), 0)

        # Test non-blocking capture with stderr (though we can still check stdout with a non-blocking capture...)
        with helpers.ShellCapture('echo "Today was not the day" > /dev/stderr && sleep 1', blocking=False) as fn:
            # Pre-execution test
            self.assertFalse(fn.is_complete())
            self.assertEqual(fn.get_stdout(), '')
            self.assertEqual(fn.get_stderr(), '')
            self.assertEqual(fn.get_exit_status(), None)

            fn()
            self.assertFalse(fn.is_complete())
            time.sleep(1)
            self.assertEqual(fn.get_stdout(), '')
            self.assertEqual(fn.get_stderr(), "Today was not the day\n")
            self.assertEqual(fn.get_exit_status(), 0)
            self.assertTrue(fn.is_complete())

    """
    A very minimal test for the function context manager that is used
    to execute pre-install invocations of tiny Python functions
    """

    def test_LambdaCapture_helper(self):
        # Test stdout capture
        with helpers.LambdaCapture(lambda: print("Today is the day")) as fn:
            # Pre-execution test
            self.assertFalse(fn.is_complete())
            self.assertEqual(fn.get_stdout(), '')
            self.assertEqual(fn.get_stderr(), '')

            fn()
            self.assertTrue(fn.is_complete())
            self.assertEqual(fn.get_stdout(), "Today is the day\n")
            self.assertEqual(fn.get_stderr(), '')

        # Test stderr capture
        with helpers.LambdaCapture(lambda: print("Today wasn't the day", file=sys.stderr)) as fn:
            # Pre-execution test
            self.assertFalse(fn.is_complete())
            self.assertEqual(fn.get_stdout(), '')
            self.assertEqual(fn.get_stderr(), '')

            fn()
            self.assertTrue(fn.is_complete())
            self.assertEqual(fn.get_stdout(), '')
            self.assertEqual(fn.get_stderr(), "Today wasn't the day\n")

        # Test non-blocking execution
        with helpers.LambdaCapture(lambda: time.sleep(1), blocking=False) as fn:
            fn()
            self.assertFalse(fn.is_complete())
            time.sleep(2)
            self.assertTrue(fn.is_complete())
