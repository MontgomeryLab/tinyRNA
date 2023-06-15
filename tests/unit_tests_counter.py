import os
import unittest

from unittest.mock import patch, mock_open

import unit_test_helpers as helpers
import tiny.rna.counter.counter as counter
from tiny.rna.util import ReadOnlyDict
from tiny.rna.configuration import ConfigBase

resources = "./testdata/counter"


class CounterTests(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.gff_file = f"{resources}/gff/identity_choice_test.gff3"
        self.short_gff_file = f"{resources}/gff/single.gff3"
        self.short_gff = helpers.read(self.short_gff_file)

        self.sam_file = f"{resources}/sam/identity_choice_test.sam"
        self.short_sam_file = f"{resources}/sam/single.sam"
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
            'Overlap':   "Partial",
            'Mismatch':  "",
            'Strand':    "antisense",
            "nt5end":    '"C,G,U"',  # Needs to be double-quoted due to commas
            'Length':    "all",
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
            'Overlap':   _row['Overlap'].lower(),
            'Mismatch':  _row['Mismatch'],
            'Strand':    _row['Strand'],
            'nt5end':    _row["nt5end"].upper().translate({ord('U'): 'T'}),
            'Length':    _row['Length'],
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

    def get_mock_samples_sheet_path(self):
        return os.path.split('/dev/null/samples.csv')

    def get_standalone_args(self):
        return ReadOnlyDict({'in_pipeline': False, 'autodoc_dir': '.'})

    def get_pipeline_args(self):
        return ReadOnlyDict({'in_pipeline': True})
    
    def get_loaded_samples_row(self, row, exp_file):
        return [{
            'Name': "_rep_".join(row[i] for i in ["Group", "Replicate"]),
            'File': exp_file,
            'Norm': row['Normalization']
        }]
        
    # === TESTS ===

    """Does load_samples correctly parse a single record samples.csv in standalone mode?"""

    def test_load_samples_single_cmd(self):
        inp_file = "test.sam"  # sam or bam
        row = dict(self.csv_samp_row_dict, File=inp_file)
        csv = self.csv("samples.csv", [row])

        mock_dir, mock_file = self.get_mock_samples_sheet_path()
        exp_file = ConfigBase.joinpath(mock_dir, inp_file)

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            with patch('tiny.rna.configuration.os.path.isfile', return_value=True):
                args = self.get_standalone_args()
                inputs_standalone = counter.load_samples(f"{mock_dir}/{mock_file}", args)

        expected_result = self.get_loaded_samples_row(row, exp_file)
        self.assertListEqual(inputs_standalone, expected_result)

    """Does load_samples correctly parse a single record samples.csv in pipeline mode?"""

    def test_load_samples_single_pipeline(self):
        inp_file = "test.fastq"
        row = dict(self.csv_samp_row_dict, File=inp_file)
        csv = self.csv("samples.csv", [row])

        mock_dir, mock_file = self.get_mock_samples_sheet_path()
        inp_root, _ = os.path.splitext(os.path.basename(inp_file))
        exp_file = f"{inp_root}_aligned_seqs.sam"

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            with patch('tiny.rna.configuration.os.path.isfile', return_value=True):
                args = self.get_pipeline_args()
                inputs_pipeline = counter.load_samples(f"{mock_dir}/{mock_file}", args)

        expected_result = self.get_loaded_samples_row(row, exp_file)
        self.assertEqual(inputs_pipeline, expected_result)

    """Does load_samples correctly handle duplicate samples? Duplicates are forbidden."""

    def test_load_samples_duplicate(self):
        # Same sample file but all other fields differ
        row1 = dict(self.csv_samp_row_dict, File="test.sam")
        row2 = dict(self.csv_samp_row_dict, File="test.sam", **{k: v + "_" for k, v in row1.items() if k != "File"})
        csv = self.csv("samples.csv", [row1, row2])

        expected_error = r"Alignment files cannot be listed more than once in .* \(row 2\)"

        with self.assertRaisesRegex(AssertionError, expected_error):
            with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
                with patch('tiny.rna.configuration.os.path.isfile', return_value=True):
                    dummy_file = '/dev/null'
                    args = self.get_standalone_args()
                    counter.load_samples(dummy_file, args)

    """Does load_samples throw ValueError if sample file does not have a .sam or .bam extension in standalone mode?"""

    def test_load_samples_bad_extension(self):
        bad = "./bad_extension.xyz"
        row = dict(self.csv_samp_row_dict, **{'File': bad})
        csv = self.csv("samples.csv", [row])

        expected_error = r"Files in .* must have a .sam or .bam extension \(row 1\)"

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            with patch('tiny.rna.configuration.os.path.isfile', return_value=True):
                with self.assertRaisesRegex(AssertionError, expected_error):
                    dummy_file = '/dev/null'
                    args = self.get_standalone_args()
                    counter.load_samples(dummy_file, args)

    """Does load_config write an autodoc copy of the Features Sheet to the Run Directory in standalone mode only?"""

    def test_load_config_autodoc(self):
        csv = self.csv("features.csv", [self.csv_feat_row_dict.copy()])
        dummy_file = '/dev/null'

        # As a pipeline step
        args = self.get_pipeline_args()
        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)), \
                patch('tiny.rna.configuration.shutil') as sh:
            ruleset1 = counter.load_config(dummy_file, args)
            sh.copyfile.assert_not_called()

        # As a standalone step
        args = self.get_standalone_args()
        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)), \
                patch('tiny.rna.configuration.shutil') as sh:
            ruleset2 = counter.load_config(dummy_file, args)
            sh.copyfile.assert_called_once()

        self.assertListEqual(ruleset1, ruleset2)

    """Does load_config correctly parse a single-entry features.csv config file in standalone mode?"""

    def test_load_config_single_cmd(self):
        # Features CSV with a single rule/row
        row = self.csv_feat_row_dict.copy()
        csv = self.csv("features.csv", [row])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            args = self.get_standalone_args()
            ruleset = counter.load_config(dummy_file, args)

        expected_ruleset = self.parsed_feat_rule
        self.assertEqual(ruleset, expected_ruleset)

    """Does load_config correctly parse a single-entry features.csv config file for pipeline mode?"""

    def test_load_config_single_pipeline(self):
        # Features CSV with a single rule/row
        row = self.csv_feat_row_dict.copy()
        csv = self.csv("features.csv", [row])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)) as conf:
            dummy_file = '/dev/null'
            args = self.get_pipeline_args()
            ruleset = counter.load_config(dummy_file, args)

        expected_ruleset = self.parsed_feat_rule
        self.assertEqual(ruleset, expected_ruleset)

    """Does load_config correctly handle duplicate rules? Want: no duplicate rules
    (including rules that differ only by their hierarchy value)."""

    def test_load_config_duplicate_rules(self):
        # Features CSV with two duplicate rules and one that differs only by hierarchy
        row = self.csv_feat_row_dict.copy()
        row_diff_hierarchy = dict(row, Hierarchy=int(row["Hierarchy"]) + 1)
        csv = self.csv("features.csv", [row, row, row_diff_hierarchy])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_filename = '/dev/null'
            args = self.get_standalone_args()
            ruleset = counter.load_config(dummy_filename, args)

        expected_ruleset = self.parsed_feat_rule
        self.assertEqual(ruleset, expected_ruleset)

    """Does load_config convert uracil to thymine for proper matching with cDNA sequences?"""

    def test_load_config_rna_to_cDNA(self):
        row = self.csv_feat_row_dict.copy()
        row["nt5end"] = 'U'
        csv = self.csv("features.csv", [row])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            args = self.get_standalone_args()
            ruleset = counter.load_config(dummy_file, args)

        self.assertEqual(ruleset[0]['nt5end'], 'T')

    """Does load_config write an autodoc copy of the Features Sheet to the Run Directory in standalone mode only?"""

    def test_load_config_autodoc(self):
        csv = self.csv("features.csv", [self.csv_feat_row_dict.copy()])
        dummy_file = '/dev/null'

        # As a pipeline step
        args = self.get_pipeline_args()
        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)), \
                patch('tiny.rna.configuration.shutil') as sh:
            ruleset1 = counter.load_config(dummy_file, args)
            sh.copyfile.assert_not_called()

        # As a standalone step
        args = self.get_standalone_args()
        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)), \
                patch('tiny.rna.configuration.shutil') as sh:
            ruleset2 = counter.load_config(dummy_file, args)
            sh.copyfile.assert_called_once()

        self.assertListEqual(ruleset1, ruleset2)

if __name__ == '__main__':
    unittest.main()
