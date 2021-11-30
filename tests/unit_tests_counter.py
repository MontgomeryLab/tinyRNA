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

        # ID, Key, Value, Hierarchy, Strand, nt5, Length, Match, Source
        self.csv_feat_row_dict = {'Key': "Class", 'Value': "CSR", 'Name': "Alias", 'Hierarchy': "1",
                                  'Strand': "antisense", 'nt5': '"C,G,U"', 'Length': "all", 'Match': "Partial",
                                  'Source': "./testdata/cel_ws279/c_elegans_WS279_chr1.gff3"}
                                   # nt5 needs to be double quoted since it contains commas

        # Identity, Hierarchy, Strand, nt5, Length, Strict
        row = self.csv_feat_row_dict
        self.feat_rule = [{
            'Identity': (row['Key'], row['Value']),
            'Hierarchy': int(row['Hierarchy']),
            'Strand': row['Strand'],
            'nt5end': row['nt5'].upper().translate({ord('U'): 'T'}).replace('"', ''),  # Undo csv comma quoting
            'Length': row['Length'],
            'Strict': row['Match'] == 'Full'
        }]

        self.csv_samp_row_dict = {'file': "test_file.fastq", 'group': "test_group", 'rep': "0"}

    # === HELPERS ===

    @staticmethod
    def csv(type, rows):
        header = "\uFEFF"
        if type == "features.csv":
            # header = "\uFEFFName Attribute,Attribute Key,Attribute Value,Hierarchy,Strand (sense/antisense/both),5' End Nucleotide,Length,Match,Feature Source"
            header = "\uFEFFSelect for...,with value...,Alias by...,Hierarchy,Strand (sense/antisense/both),5' End Nucleotide,Length,Match,Feature Source"
        elif type == "samples.csv":
            header = "\uFEFFInput FASTQ Files,Sample/Group Name,Replicate number,Control"

        return '\n'.join([header, *map(','.join, rows)])
    
    def feat_csv_test_row(self):
        return ','.join(self.csv_feat_row_dict.values())
        
    # === TESTS ===

    """Does load_samples correctly parse a single record samples.csv for command line invocation?"""

    def test_load_samples_single_cmd(self):
        dummy_file = '/dev/null'
        inp_file = "test.fastq"
        exp_file = from_here(dummy_file, "test_aligned_seqs.sam")

        row = {'File': inp_file, 'Group': "test_group", 'Rep': "0"}
        csv = self.csv("samples.csv", [row.values()])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            inputs_step = counter.load_samples(dummy_file, is_pipeline=False)

        expected_lib_name = f"{row['Group']}_rep_{row['Rep']}"
        expected_result = [{'File': exp_file, 'Name': expected_lib_name}]
        self.assertEqual(inputs_step, expected_result)

    """Does load_samples correctly parse a single record samples.csv for pipeline invocation?"""

    def test_load_samples_single_pipeline(self):
        dummy_file = '/dev/null'
        inp_file = "test.fastq"
        exp_file = "test_aligned_seqs.sam"

        row = {'File': inp_file, 'Group': "test_group", 'Rep': "0"}
        csv = self.csv("samples.csv", [row.values()])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            inputs_pipeline = counter.load_samples(dummy_file, is_pipeline=True)

        expected_lib_name = f"{row['Group']}_rep_{row['Rep']}"
        expected_result = [{'File': exp_file, 'Name': expected_lib_name}]
        self.assertEqual(inputs_pipeline, expected_result)

    """Does load_samples correctly handle duplicate samples? There should be no duplicates."""

    def test_load_samples_duplicate(self):
        row = {'File': "test.fastq", 'Group': "N/A", 'Rep': "N/A"}
        csv = self.csv("samples.csv", [row.values(), row.values()])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            inputs = counter.load_samples(dummy_file, False)

        self.assertEqual(len(inputs), 1)

    """Does load_samples correctly handle SAM filenames?"""

    def test_load_samples_sam(self):
        sam_filename = "/fake/absolute/path/sample.sam"
        row = {'File': sam_filename, 'Group': "test_group", 'Rep': "0"}
        csv = self.csv("samples.csv", [row.values()])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            inputs = counter.load_samples(dummy_file, False)

        expected_lib_name = f"{row['Group']}_rep_{row['Rep']}"
        expected_result = [{'File': sam_filename, 'Name': expected_lib_name}]

        self.assertEqual(inputs, expected_result)

    """Does load_samples throw ValueError if a non-absolute path to a SAM file is provided?"""

    def test_load_samples_nonabs_path(self):
        bad = "./dne.sam"
        row = [bad, "test_group", "0"]
        csv = self.csv("samples.csv", [row])

        expected_error = "The following file must be expressed as an absolute path:\n" + bad

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            with self.assertRaisesRegex(ValueError, expected_error):
                dummy_file = '/dev/null'
                counter.load_samples(dummy_file, False)

    """Does load_samples throw ValueError if sample filename does not have a .fastq or .sam extension?"""

    def test_load_samples_bad_extension(self):
        bad = "./bad_extension.xyz"
        row = [bad, "test_group", "0"]
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
        row = self.csv_feat_row_dict.values()
        csv = self.csv("features.csv", [row])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            ruleset, gff_files = counter.load_config(dummy_file, is_pipeline=False)

        r = self.csv_feat_row_dict
        expected_gff_file = from_here(dummy_file, r['Source'])
        expected_gff_ret = defaultdict(list, zip([expected_gff_file], [[r['Name']]]))
        expected_ruleset = self.feat_rule

        self.assertEqual(gff_files, expected_gff_ret)
        self.assertEqual(ruleset, expected_ruleset)

    """Does load_config correctly parse a single-entry features.csv config file for pipeline invocation?"""

    def test_load_config_single_pipeline(self):
        # Features CSV with a single rule/row
        row = self.csv_feat_row_dict.values()
        csv = self.csv("features.csv", [row])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            ruleset, gff_files = counter.load_config(dummy_file, is_pipeline=True)

        r = self.csv_feat_row_dict
        expected_gff_file = os.path.basename(r['Source'])
        expected_gff_ret = defaultdict(list, zip([expected_gff_file], [[r['Name']]]))
        expected_ruleset = self.feat_rule

        self.assertEqual(gff_files, expected_gff_ret)
        self.assertEqual(ruleset, expected_ruleset)

    """Does load_config correctly handle duplicate rules? Want: no duplicate rules and no duplicate Name Attributes."""

    def test_load_config_duplicate_rules(self):
        # Features CSV with two duplicate rules/rows
        row = self.csv_feat_row_dict.values()
        csv = self.csv("features.csv", [row, row])  # Duplicate rows
        
        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_filename = '/dev/null'
            ruleset, gff_files = counter.load_config(dummy_filename, False)

        r = self.csv_feat_row_dict
        expected_gff_file = from_here(dummy_filename, r['Source'])
        expected_gff_ret = defaultdict(list, zip([expected_gff_file], [[r['Name']]]))
        expected_ruleset = self.feat_rule

        self.assertEqual(gff_files, expected_gff_ret)
        self.assertEqual(ruleset, expected_ruleset)

    """Does load_config convert uracil to thymine for proper matching with cDNA sequences?"""

    def test_load_config_rna_to_cDNA(self):
        row = self.csv_feat_row_dict.copy()
        row['nt5'] = 'U'
        csv = self.csv("features.csv", [row.values()])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            ruleset, _ = counter.load_config(dummy_file, False)

        self.assertEqual(ruleset[0]['nt5end'], 'T')

    """Does load_config properly screen for "ID" Name Attributes?"""

    def test_load_config_id_name_attr(self):
        row = self.csv_feat_row_dict.copy()
        row['Name'] = 'ID'
        csv = self.csv("features.csv", [row.values()])

        with patch('tiny.rna.configuration.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            _, gff_files = counter.load_config(dummy_file, False)

        # Expect {file: [empty Name Attribute list]}
        from_dummy = from_here(dummy_file, row['Source'])
        expected = defaultdict(list, zip([from_dummy], [[]]))
        self.assertEqual(gff_files, expected)


if __name__ == '__main__':
    unittest.main()
