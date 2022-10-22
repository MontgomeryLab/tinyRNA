import contextlib
import io
import os
import unittest
from unittest.mock import patch, mock_open, call

from tiny.rna.configuration import Configuration, SamplesSheet
from unit_test_helpers import csv_factory


class BowtieIndexesTest(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.root_cfg_dir = os.path.abspath("../tiny/templates")
        self.run_config = self.root_cfg_dir + "/run_config_template.yml"
        self.paths = self.root_cfg_dir + "/paths.yml"

        self.default_prefix = os.path.join(
            self.root_cfg_dir,
            Configuration(self.run_config)['run_directory'],
            "bowtie-build/ram1"
        )
        self.maxDiff = 1522

    """============ Helper functions ============"""

    def config_with(self, prefs):
        config = Configuration(self.run_config)
        for key, val in prefs.items():
            config[key] = val
        return config

    def bt_idx_files_from_prefix(self, prefix):
        return [
            {'path': f"{prefix}.{subext}.ebwt", 'class': 'File'}
            for subext in ['1', '2', '3', '4', 'rev.1', 'rev.2']
        ]

    """================ Tests =================="""

    """Does get_ebwt_prefix() produce the expected prefix path?"""

    def test_get_ebwt_prefix(self):
        config = Configuration(self.run_config)
        actual_prefix = config.get_ebwt_prefix()
        expected_prefix = self.default_prefix

        self.assertEqual(actual_prefix, expected_prefix)

    """Does get_ebwt_prefix() throw an error if reference genome files aren't provided?"""

    def test_get_ebwt_prefix_no_genome(self):
        config = Configuration(self.run_config)
        config['reference_genome_files'] = None

        with self.assertRaises(ValueError):
            config.get_ebwt_prefix()

    """Does get_bt_index_files() output the paths of indexes that have already been built?"""

    def test_get_bt_index_files_prebuilt_indexes(self):
        config = self.config_with({'run_bowtie_build': False})
        prefix = config.paths['ebwt'] = os.path.abspath("./testdata/counter/validation/ebwt/ram1")
        expected = self.bt_idx_files_from_prefix(prefix)
        self.assertListEqual(config.get_bt_index_files(), expected)

    """Does get_bt_index_files() output the paths of the index files that are expected
    to be built from the reference genome?"""

    def test_get_bt_index_files_unbuilt_indexes_with_genome(self):
        config = self.config_with({'run_bowtie_build': True})
        prefix = config.paths['ebwt'] = "mock_prefix"
        expected = self.bt_idx_files_from_prefix(prefix)
        self.assertListEqual(config.get_bt_index_files(), expected)

    """Does get_bt_index_files() produce an error and quit when index files are
    missing and a reference genome has not been provided?"""

    def test_get_bt_index_files_missing_indexes_without_genome(self):
        config = self.config_with({'run_bowtie_build': False, 'reference_genome_files': None})
        prefix = config.paths['ebwt'] = "missing"
        errmsg = '\n'.join([
            "The following Bowtie index file couldn't be found:",
            "\t" + f"{prefix}.1.ebwt",
            "\nPlease either correct your ebwt prefix or add reference genomes in the Paths File."
        ])

        with self.assertRaisesRegex(SystemExit, errmsg):
            config.get_bt_index_files()

    """Does get_bt_index_files() produce an error without quitting when index files
    are missing but a reference genome was provided, and does it return the list of
    index files that will be built from the genome?"""

    def test_get_bt_index_files_missing_indexes_with_genome(self):
        config = self.config_with({'run_bowtie_build': False})
        bad_prefix = config.paths['ebwt'] = "missing"
        genome_prefix = self.default_prefix

        expected_files = self.bt_idx_files_from_prefix(genome_prefix)
        expected_error = '\n'.join([
            "The following Bowtie index file couldn't be found:",
            "\t" + f"{bad_prefix}.1.ebwt",
            "\nIndexes will be built from your reference genome files during this run.",
            ""
        ])

        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            actual = config.get_bt_index_files()

        self.assertEqual(stderr.getvalue(), expected_error)
        self.assertListEqual(actual, expected_files)

    """Does verify_bowtie_build_outputs() update the paths in ["bt_index_files"] and rewrite 
    these changes to the processed Run Config if long indexes were produced? Does it also
    write to the Paths File to update the new ebwt prefix?"""

    def test_verify_bowtie_build_outputs(self):
        ebwt_short = ["1.ebwt", "2.ebwt", "3.ebwt"]
        ebwt_long = ["1.ebwtl", "2.ebwtl", "3.ebwtl"]
        run_conf_ebwt = [Configuration.cwl_file(f, verify=False) for f in ebwt_short]
        expected_ebwt = [Configuration.cwl_file(f, verify=False) for f in ebwt_long]

        config = self.config_with({'bt_index_files': run_conf_ebwt})

        with patch('tiny.rna.configuration.open', mock_open()) as mo, \
                patch('tiny.rna.configuration.glob', return_value=ebwt_long) as g:
            config.verify_bowtie_build_outputs()

        expected_writes = [
            call(self.paths, 'w'),
            call(os.path.join(self.root_cfg_dir, config['run_directory'], os.path.basename(self.run_config)), 'w')
        ]

        self.assertListEqual(config['bt_index_files'], expected_ebwt)
        self.assertListEqual(mo.call_args_list, expected_writes)


class SamplesSheetTest(unittest.TestCase):

    """Does SamplesSheet catch multi-assignment of control condition?"""

    def test_validate_control_group(self):
        sheet = csv_factory("samples.csv", [
            {'File': '1.fastq', 'Group': 'G1', 'Replicate': '1', 'Control': True, 'Normalization': ''},  # Good
            {'File': '2.fastq', 'Group': 'G1', 'Replicate': '2', 'Control': True, 'Normalization': ''},  # Good
            {'File': '3.fastq', 'Group': 'G2', 'Replicate': '1', 'Control': True, 'Normalization': ''}   # Bad
        ])                                                                 # ^^^

        exp_contains = r".*(multiple control conditions).*"
        with self.assertRaisesRegex(AssertionError, exp_contains), \
                patch('tiny.rna.configuration.open', mock_open(read_data=sheet)), \
                patch('tiny.rna.configuration.os.path.isfile', return_value=True):
            SamplesSheet('mock_filename')

    """Does SamplesSheet catch duplicate entries for the same group and rep?"""
    def test_validate_group_rep(self):
        sheet = csv_factory("samples.csv", [
            {'File': '1.fastq', 'Group': 'G1', 'Replicate': '1', 'Control': True, 'Normalization': ''},  # Good
            {'File': '2.fastq', 'Group': 'G1', 'Replicate': '2', 'Control': True, 'Normalization': ''},  # Good
            {'File': '3.fastq', 'Group': 'G1', 'Replicate': '2', 'Control': True, 'Normalization': ''}   # Bad
        ])                             # ^^^                ^^^

        exp_contains = r".*(same group and replicate).*"
        with self.assertRaisesRegex(AssertionError, exp_contains), \
                patch('tiny.rna.configuration.open', mock_open(read_data=sheet)), \
                patch('tiny.rna.configuration.os.path.isfile', return_value=True):
            SamplesSheet('mock_filename')

    """Does SamplesSheet catch fastq files that don't exist, have a bad file extension, or are listed more than once?"""
    def test_validate_fastq_filepath(self):
        csv_rows = [
            {'File': '1.fastq', 'Group': 'G1', 'Replicate': '1', 'Control': True, 'Normalization': ''},  # Good
            {'File': '1.fastq', 'Group': 'G1', 'Replicate': '2', 'Control': True, 'Normalization': ''},  # Bad
            {'File': '2.fasta', 'Group': 'G2', 'Replicate': '1', 'Control': True, 'Normalization': ''}   # Bad
        ]           # ^^^^^^^
        sheet = csv_factory("samples.csv", csv_rows)

        # File doesn't exist
        exp_contains = r".*(was not found).*"
        with self.assertRaisesRegex(AssertionError, exp_contains), \
                patch('tiny.rna.configuration.open', mock_open(read_data=sheet)):
            SamplesSheet('mock_filename')

        # Duplicate filename
        exp_contains = r".*(listed more than once).*"
        with self.assertRaisesRegex(AssertionError, exp_contains), \
                patch('tiny.rna.configuration.open', mock_open(read_data=sheet)), \
                patch('tiny.rna.configuration.os.path.isfile', return_value=True):
            SamplesSheet('mock_filename')

        # Bad file extension
        exp_contains = r".*(\.fastq\(\.gz\) extension).*"
        csv_rows.pop(0)
        sheet = csv_factory("samples.csv", csv_rows)
        with self.assertRaisesRegex(AssertionError, exp_contains), \
                patch('tiny.rna.configuration.open', mock_open(read_data=sheet)), \
                patch('tiny.rna.configuration.os.path.isfile', return_value=True):
            SamplesSheet('mock_filename')


if __name__ == '__main__':
    unittest.main()
