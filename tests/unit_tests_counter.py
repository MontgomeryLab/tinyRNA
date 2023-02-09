import io
import os
import unittest

from unittest.mock import patch, mock_open
from collections import defaultdict

import unit_test_helpers as helpers
import tiny.rna.counter.counter as counter
from tiny.rna.util import from_here

resources = "./testdata/counter"


class CounterTests(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.gff_file = f"{resources}/identity_choice_test.gff3"
        self.short_gff_file = f"{resources}/single.gff3"
        self.short_gff = helpers.read(self.short_gff_file)

        self.sam_file = f"{resources}/identity_choice_test.sam"
        self.short_sam_file = f"{resources}/single.sam"
        self.short_sam = helpers.read(self.short_sam_file)

        self.strand = {'sense': tuple('+'), 'antisense': tuple('-'), 'both': ('+', '-')}
        self.csv = staticmethod(helpers.csv_factory)

        # Represents an unparsed Features Sheet row
        # Key is the user-facing column header
        self.csv_feat_row_dict = {
            'Key':       "Class",
            'Value':     "CSR",
            'Class':     "",
            'Filter_s':  "",
            'Filter_t':  "",
            'Hierarchy': "1",
            'Strand':    "antisense",
            "nt5end":    '"C,G,U"',  # Needs to be double-quoted due to commas
            'Length':    "all",
            'Overlap':   "Partial"
        }

        # Represents the parsed Features Sheet row above
        # Key is the internal short name
        _row = self.csv_feat_row_dict
        self.parsed_feat_rule = [{
            'Identity':  (_row['Key'], _row['Value']),
            'Class':     _row['Class'],
            'Filter_s':  _row['Filter_s'],
            'Filter_t':  _row['Filter_t'],
            'Hierarchy': int(_row['Hierarchy']),
            'Strand':    _row['Strand'],
            'nt5end':    _row["nt5end"].upper().translate({ord('U'): 'T'}),
            'Length':    _row['Length'],
            'Overlap':   _row['Overlap'].lower()
        }]

        # Represents an unparsed Samples Sheet row
        # Key is the user-facing column header
        self.csv_samp_row_dict = {
            'File':          "test_file.fastq",
            'Group':         "test_group",
            'Replicate':     "0",
            'Control':       "",
            'Normalization': ""
        }

        # This is the same Samples Sheet row above, but with internal names
        # It does NOT represent the parsed result of loading the Samples Sheet
        _row = self.csv_samp_row_dict
        self.parsed_samp_rule = {
            'File':          _row['File'],
            'Group':         _row['Group'],
            'Replicate':     _row['Replicate'],
            'Control':       _row['Control'],
            'Normalization': _row['Normalization']
        }

    # === HELPERS ===
    
    def get_loaded_samples_row(self, row, exp_file):
        return [{
            'Name': "_rep_".join(row[i] for i in ["Group", "Replicate"]),
            'File': exp_file,
            'Norm': row['Normalization']
        }]
        
    # === TESTS ===

    """Does load_samples correctly parse a single record samples.csv for command line invocation?"""

    def test_load_samples_single_cmd(self):
        mock_samp_sheet_path = '/dev/null'
        inp_file = "test.fastq"
        exp_file = from_here(mock_samp_sheet_path, "test_aligned_seqs.sam")

        row = dict(self.csv_samp_row_dict, **{'File': inp_file})
        csv = self.csv("samples.csv", [row])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            inputs_step = counter.load_samples(mock_samp_sheet_path, is_pipeline=False)

        expected_result = self.get_loaded_samples_row(row, exp_file)
        self.assertEqual(inputs_step, expected_result)

    """Does load_samples correctly parse a single record samples.csv for pipeline invocation?"""

    def test_load_samples_single_pipeline(self):
        mock_samp_sheet_path = '/dev/null'
        inp_file = "test.fastq"
        exp_file = "test_aligned_seqs.sam"

        row = dict(self.csv_samp_row_dict, **{'File': inp_file})
        csv = self.csv("samples.csv", [row])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            inputs_pipeline = counter.load_samples(mock_samp_sheet_path, is_pipeline=True)

        expected_result = self.get_loaded_samples_row(row, exp_file)
        self.assertEqual(inputs_pipeline, expected_result)

    """Does load_samples correctly handle duplicate samples? There should be no duplicates."""

    def test_load_samples_duplicate(self):
        row = self.csv_samp_row_dict.copy()
        csv = self.csv("samples.csv", [row, row])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            inputs = counter.load_samples(dummy_file, False)

        self.assertEqual(len(inputs), 1)

    """Does load_samples correctly handle SAM filenames?"""

    def test_load_samples_sam(self):
        sam_filename = "/fake/absolute/path/sample.sam"
        row = dict(self.csv_samp_row_dict, **{'File': sam_filename})
        csv = self.csv("samples.csv", [row])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            inputs = counter.load_samples(dummy_file, is_pipeline=False)

        expected_result = self.get_loaded_samples_row(row, sam_filename)
        self.assertEqual(inputs, expected_result)

    """Does load_samples throw ValueError if a non-absolute path to a SAM file is provided?"""

    def test_load_samples_nonabs_path(self):
        bad = "./dne.sam"
        row = dict(self.csv_samp_row_dict, **{'File': bad})
        csv = self.csv("samples.csv", [row])

        expected_error = "The following file must be expressed as an absolute path:\n" + bad

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            with self.assertRaisesRegex(ValueError, expected_error):
                dummy_file = '/dev/null'
                counter.load_samples(dummy_file, False)

    """Does load_samples throw ValueError if sample filename does not have a .fastq or .sam extension?"""

    def test_load_samples_bad_extension(self):
        bad = "./bad_extension.xyz"
        row = dict(self.csv_samp_row_dict, **{'File': bad})
        csv = self.csv("samples.csv", [row])

        expected_error = r"The filenames defined in your Samples Sheet must have a \.fastq\(\.gz\) or \.sam extension\.\n" \
                         r"The following filename contained neither\:\n" + bad

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            with self.assertRaisesRegex(ValueError, expected_error):
                dummy_file = '/dev/null'
                counter.load_samples(dummy_file, False)

    """Does load_config correctly parse a single-entry features.csv config file for command line invocation?"""

    def test_load_config_single_cmd(self):
        # Features CSV with a single rule/row
        row = self.csv_feat_row_dict.copy()
        csv = self.csv("features.csv", [row])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            ruleset = counter.load_config(dummy_file, is_pipeline=False)

        expected_ruleset = self.parsed_feat_rule
        self.assertEqual(ruleset, expected_ruleset)

    """Does load_config correctly parse a single-entry features.csv config file for pipeline invocation?"""

    def test_load_config_single_pipeline(self):
        # Features CSV with a single rule/row
        row = self.csv_feat_row_dict.copy()
        csv = self.csv("features.csv", [row])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            ruleset = counter.load_config(dummy_file, is_pipeline=True)

        expected_ruleset = self.parsed_feat_rule
        self.assertEqual(ruleset, expected_ruleset)

    """Does load_config correctly handle duplicate rules? Want: no duplicate rules and no duplicate Name Attributes."""

    def test_load_config_duplicate_rules(self):
        # Features CSV with two duplicate rules/rows
        row = self.csv_feat_row_dict.copy()
        csv = self.csv("features.csv", [row, row])  # Duplicate rows
        
        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_filename = '/dev/null'
            ruleset = counter.load_config(dummy_filename, False)

        expected_ruleset = self.parsed_feat_rule
        self.assertEqual(ruleset, expected_ruleset)

    """Does load_config convert uracil to thymine for proper matching with cDNA sequences?"""

    def test_load_config_rna_to_cDNA(self):
        row = self.csv_feat_row_dict.copy()
        row["nt5end"] = 'U'
        csv = self.csv("features.csv", [row])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            ruleset = counter.load_config(dummy_file, False)

        self.assertEqual(ruleset[0]['nt5end'], 'T')


if __name__ == '__main__':
    unittest.main()
