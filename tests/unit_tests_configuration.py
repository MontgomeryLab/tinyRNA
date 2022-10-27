import contextlib
import io
import os
import random
import unittest
from unittest.mock import patch, mock_open, call

from tiny.rna.configuration import Configuration, SamplesSheet, PathsFile
from unit_test_helpers import csv_factory, paths_template_file, make_paths_file


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


class PathsFileTest(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.template_dir = os.path.dirname(paths_template_file)

    """Does PathsFile automatically resolve paths when queried?"""

    def test_getitem_single(self):
        config = make_paths_file()
        config['mock_parameter'] = "./some/../file"
        self.assertEqual(config['mock_parameter'], os.path.join(self.template_dir, "file"))

    """Does PathsFile automatically resolve lists of paths when queried? What if some list items are
    strings, others are mappings with a "path" key, and others are None/empty?"""

    def test_getitem_group(self):
        config = make_paths_file()
        config.groups = ('mock_parameter',)

        mapping_1 = {'path': "./some/../file", "other_key": "irrelevant"}
        mapping_2 = {'path': "../templates/another_file"}
        path_string = "../../START_HERE/reference_data/ram1.gff3"

        config['mock_parameter'] = [mapping_1, mapping_2, path_string, None, '']

        # Notice: only 'path' is modified in mappings and empty entries are returned unmodified
        expected = [
            dict(mapping_1, path=os.path.join(self.template_dir, "file")),
            dict(mapping_2, path=os.path.join(self.template_dir, "another_file")),
            os.path.normpath(os.path.join(self.template_dir, path_string)),
            None,
            ''
        ]

        self.assertListEqual(config['mock_parameter'], expected)

    """Does PathsFile check for required parameters?"""

    def test_validate_required_parameters(self):
        config = make_paths_file()
        for key in PathsFile.required:
            oldval = config[key]

            with self.assertRaisesRegex(AssertionError, r".*(parameters are required).*"):
                config[key] = ''
                config.validate_paths()

            with self.assertRaisesRegex(AssertionError, r".*(parameters are required).*"):
                config[key] = None
                config.validate_paths()

            with self.assertRaisesRegex(AssertionError, r".*(parameters are required).*"):
                del config.config[key]
                config.validate_paths()

            config[key] = oldval

    """Does PathsFile check that at least one GFF file has been provided?"""

    def test_validate_gff_files(self):
        config = make_paths_file()

        with self.assertRaisesRegex(AssertionError, r".*(At least one GFF).*"):
            config['gff_files'] = [{'path': "", 'alias': []}]
            config.validate_paths()

        with self.assertRaisesRegex(AssertionError, r".*(At least one GFF).*"):
            config['gff_files'] = [{'irrelevant': "value"}]
            config.validate_paths()

    """Does PathsFile check for missing files for single entry parameters?"""

    def test_validate_missing_file_single(self):
        config = make_paths_file()
        key = random.choice(PathsFile.single)
        config[key] = '/dev/null/file_dne'

        with self.assertRaisesRegex(AssertionError, f".*(file provided for {key}).*"):
            config.validate_paths()

    """Does PathsFile check for missing files under list-type parameters?"""

    def test_validate_missing_file_group(self):
        config = make_paths_file()
        bad_path = "/dev/null/file_dne"
        config.append_to('gff_files', {'path': bad_path, 'alias': []})

        # Config now has one good entry and one bad entry; make sure the bad entry was specifically reported
        with self.assertRaisesRegex(AssertionError, f".*(file provided under gff_files).*\n\t{bad_path}"):
            config.validate_paths()

    """Does PathsFile detect a backward compatibility issue?"""

    def test_backward_incompatibility(self):
        config = make_paths_file()
        del config.config['gff_files']

        with self.assertRaisesRegex(AssertionError, r".*(check the release notes).*"):
            config.check_backward_compatibility()

if __name__ == '__main__':
    unittest.main()
