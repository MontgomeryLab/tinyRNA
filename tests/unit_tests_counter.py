import unittest
import HTSeq

from unittest.mock import patch, mock_open, call
from collections import defaultdict

import tests.unit_test_helpers as helpers
import aquatx.srna.counter.counter as counter
from aquatx.srna.counter.hts_parsing import Alignment, read_SAM
from aquatx.srna.util import from_here

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
        self.csv_feat_row_dict = {'Name': "Alias", 'Key': "Class", 'Value': "CSR", 'Hierarchy': "1",
                                  'Strand': "antisense", 'nt5': '"C,G,U"', 'Length': "all", 'Match': "Partial",
                                  'Source': "./testdata/cel_ws279/c_elegans.PRJNA13758.WS279.chr1.gff3"}
                                   # nt5 needs to be double quoted since it contains commas

        # Identity, Hierarchy, Strand, nt5, Length, Strict
        row = self.csv_feat_row_dict
        self.feat_rule = [{
            'Identity': (row['Key'], row['Value']),
            'Hierarchy': int(row['Hierarchy']),
            'Strand': self.strand[row['Strand']],
            'nt5': row['nt5'].upper().translate({ord('U'): 'T'}).replace('"', ''),  # Undo csv comma quoting
            'Length': row['Length'],
            'Strict': row['Match'] == 'Full'
        }]

        self.csv_samp_row_dict = {'file': "test_file.fastq", 'group': "test_group", 'rep': "0"}

    # === HELPERS ===

    @staticmethod
    def csv(type, rows):
        header = "\uFEFF"
        if type == "features.csv":
            header = "\uFEFFName Attribute,Attribute Key,Attribute Value,Hierarchy,Strand (sense/antisense/both),5' End Nucleotide,Length,Match,Feature Source"
        elif type == "samples.csv":
            header = "\uFEFFInput FASTQ Files,Sample/Group Name,Replicate number"

        return '\n'.join([header, *map(','.join, rows)])
    
    def feat_csv_test_row(self):
        return ','.join(self.csv_feat_row_dict.values())
        
    # === TESTS ===

    """Does load_samples correctly parse a single record samples.csv?"""

    def test_load_samples_single(self):
        dummy_file = '/dev/null'
        inp_file = "test.fastq"
        exp_file = from_here(dummy_file, "test_aligned_seqs.sam")

        row = {'File': inp_file, 'Group': "test_group", 'Rep': "0"}
        csv = self.csv("samples.csv", [row.values()])

        with patch('aquatx.srna.counter.counter.open', mock_open(read_data=csv)):
            inputs = counter.load_samples(dummy_file)

        expected_lib_name = f"{row['Group']}_rep_{row['Rep']}"
        expected_result = [{'File': exp_file, 'Name': expected_lib_name}]
        self.assertEqual(inputs, expected_result)

        """Does load_samples correctly handle duplicate samples? There should be no duplicates."""

    def test_load_samples_duplicate(self):
        row = {'File': "test.fastq", 'Group': "N/A", 'Rep': "N/A"}
        csv = self.csv("samples.csv", [row.values(), row.values()])

        with patch('aquatx.srna.counter.counter.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            inputs = counter.load_samples(dummy_file)

        self.assertEqual(len(inputs), 1)

    """Does load_samples correctly handle SAM filenames?"""

    def test_load_samples_sam(self):
        sam_filename = "/fake/absolute/path/sample.sam"
        row = {'File': sam_filename, 'Group': "test_group", 'Rep': "0"}
        csv = self.csv("samples.csv", [row.values()])

        with patch('aquatx.srna.counter.counter.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            inputs = counter.load_samples(dummy_file)

        expected_lib_name = f"{row['Group']}_rep_{row['Rep']}"
        expected_result = [{'File': sam_filename, 'Name': expected_lib_name}]

        self.assertEqual(inputs, expected_result)

    """Does load_samples throw ValueError if a non-absolute path to a SAM file is provided?"""

    def test_load_samples_nonabs_path(self):
        bad = "./dne.sam"
        row = [bad, "test_group", "0"]
        csv = self.csv("samples.csv", [row])

        expected_error = "The following file must be expressed as an absolute path:\n" + bad

        with patch('aquatx.srna.counter.counter.open', mock_open(read_data=csv)):
            with self.assertRaisesRegex(ValueError, expected_error):
                dummy_file = '/dev/null'
                counter.load_samples(dummy_file)

    """Does load_samples throw ValueError if sample filename does not have a .fastq or .sam extension?"""

    def test_load_samples_bad_extension(self):
        bad = "./bad_extension.xyz"
        row = [bad, "test_group", "0"]
        csv = self.csv("samples.csv", [row])

        expected_error = r"The filenames defined in your Samples Sheet must have a \.fastq\(\.gz\) or \.sam extension\.\n" \
                         r"The following filename contained neither\:\n" + bad

        with patch('aquatx.srna.counter.counter.open', mock_open(read_data=csv)):
            with self.assertRaisesRegex(ValueError, expected_error):
                dummy_file = '/dev/null'
                counter.load_samples(dummy_file)

    """Does load_config correctly parse a single-entry features.csv config file?"""

    def test_load_config_single(self):
        # Features CSV with a single rule/row
        row = self.csv_feat_row_dict.values()
        csv = self.csv("features.csv", [row])

        with patch('aquatx.srna.counter.counter.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            ruleset, gff_files = counter.load_config(dummy_file)

        r = self.csv_feat_row_dict
        expected_gff_file = from_here(dummy_file, r['Source'])
        expected_gff_ret = defaultdict(list, zip([expected_gff_file], [[r['Name']]]))
        expected_ruleset = self.feat_rule

        self.assertEqual(gff_files, expected_gff_ret)
        self.assertEqual(ruleset, expected_ruleset)

    """Does load_config correctly handle duplicate rules? Want: no duplicate rules and no duplicate Name Attributes."""

    def test_load_config_duplicate_rules(self):
        # Features CSV with two duplicate rules/rows
        row = self.csv_feat_row_dict.values()
        csv = self.csv("features.csv", [row, row])  # Duplicate rows
        
        with patch('aquatx.srna.counter.counter.open', mock_open(read_data=csv)):
            dummy_filename = '/dev/null'
            ruleset, gff_files = counter.load_config(dummy_filename)

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

        with patch('aquatx.srna.counter.counter.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            ruleset, _ = counter.load_config(dummy_file)

        self.assertEqual(ruleset[0]['nt5'], 'T')

    """Does load_config properly screen for "ID" Name Attributes?"""

    def test_load_config_id_name_attr(self):
        row = self.csv_feat_row_dict.copy()
        row['Name'] = 'ID'
        csv = self.csv("features.csv", [row.values()])

        with patch('aquatx.srna.counter.counter.open', mock_open(read_data=csv)):
            dummy_file = '/dev/null'
            _, gff_files = counter.load_config(dummy_file)

        # Expect {file: [empty Name Attribute list]}
        from_dummy = from_here(dummy_file, row['Source'])
        expected = defaultdict(list, zip([from_dummy], [[]]))
        self.assertEqual(gff_files, expected)

    """Do GenomicArraysOfSets slice to step intervals that overlap, even if by just one base?"""

    def test_HTSeq_iv_slice(self):
        gas = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        iva = HTSeq.GenomicInterval("I", 1, 10, "+")
        ivb = HTSeq.GenomicInterval("I", 5, 15, "+")
        ivc = HTSeq.GenomicInterval("I", 9, 20, "+")
        ivd = HTSeq.GenomicInterval("I", 2, 4, "+")
        gas[iva] += "TestA"
        gas[ivb] += "TestB"

        """
        iva:  1 |--TestA--| 10
        ivb:      5 |---TestB--| 15
        ivc:          9 |-----------| 20
        ivd:   2 |--| 4
                         ^ Single base overlap: iva ∩ ivc
        Expect:       9 |-|{B}-|-{}-| 20
                     [9, 10)   [15,20)
                         ^ {A ∩ B}
        """

        matches = list(gas[ivc].array[ivc.start:ivc.end].get_steps(values_only=True))
        matches_with_cooridnates = list(gas[ivc].steps())
        self.assertEqual(matches, [{"TestA", "TestB"}, {"TestB"}, set()])
        self.assertEqual(matches_with_cooridnates[0][0], HTSeq.GenomicInterval("I", 9,10,'+'))
        self.assertEqual(matches_with_cooridnates[1][0], HTSeq.GenomicInterval("I", 10,15,'+'))
        self.assertEqual(matches_with_cooridnates[2][0], HTSeq.GenomicInterval("I", 15,20,'+'))
        self.assertEqual(matches_with_cooridnates[2][0], HTSeq.GenomicInterval("I", 15,20,'+'))

    """What happens if we invert the previous test's strand? Do intervals behave differently?"""

    def test_HTSeq_antisense_iv_slice(self):
        gas = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        iva = HTSeq.GenomicInterval("I", 1, 10, '-')
        ivb = HTSeq.GenomicInterval("I", 5, 15, '-')
        ivc = HTSeq.GenomicInterval("I", 9, 20, '-')
        ivd = HTSeq.GenomicInterval("I", 2, 4, "-")
        gas[iva] += "TestA"
        gas[ivb] += "TestB"
        gas[ivd] += "TestD"

        """
        iva:  1 |--TestA--| 10
        ivb:      5 |---TestB--| 15
        ivc:          9 |-----------| 20
        ivd:   2 |--| 4
                         ^ Single base overlap: iva ∩ ivc
        Expect:       9 |-|{B}-|-{}-| 20
                     [9, 10)   [15,20)
                         ^ {A ∩ B}
        """

        matches = list(gas[ivc].array[ivc.start:ivc.end].get_steps(values_only=True))
        matches_with_cooridnates = list(gas[ivc].steps())
        self.assertEqual(matches, [{"TestA", "TestB"}, {"TestB"}, set()])
        self.assertEqual(matches_with_cooridnates[0][0], HTSeq.GenomicInterval("I", 9,10,'-'))
        self.assertEqual(matches_with_cooridnates[1][0], HTSeq.GenomicInterval("I", 10,15,'-'))
        self.assertEqual(matches_with_cooridnates[2][0], HTSeq.GenomicInterval("I", 15,20,'-'))
        self.assertEqual(matches_with_cooridnates[2][0], HTSeq.GenomicInterval("I", 15,20,'-'))

    """Does assign_features return features that overlap the query interval by a single base?"""

    def test_assign_features_single_base_overlap(self):
        features = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        iv_feat = HTSeq.GenomicInterval("I", 0, 2, "+")  # The "feature"
        iv_olap = HTSeq.GenomicInterval("I", 1, 2, "+")  # A single-base overlapping feature
        iv_none = HTSeq.GenomicInterval("I", 2, 3, "+")  # A non-overlapping interval

        features[iv_feat] += 'The "feature"'
        features[iv_none] += "Non-overlapping feature"

        olap_alignment = Alignment(iv_olap, "Single base overlap", b"A")

        with patch.object(counter, "selector") as selector:
            counter.features = features
            counter.assign_features(olap_alignment)

        expected_match_list = [(1, 2, {'The "feature"'})]
        selector.choose.assert_called_once_with(expected_match_list, olap_alignment)

    """Does assign_features correctly handle alignments with zero feature matches?"""

    def test_assign_features_no_match(self):
        features = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        iv_feat = HTSeq.GenomicInterval("I", 0, 2, "+")
        iv_none = HTSeq.GenomicInterval("I", 2, 3, "+")

        features[iv_feat] += "Should not match"
        none_alignment = Alignment(iv_none, "Non-overlap", b"A")

        with patch.object(counter, "selector") as selector:
            counter.features = features
            counter.assign_features(none_alignment)

        selector.choose.assert_not_called()

    """Does count_reads call the right functions when handling a single record library?"""

    @patch.object(counter, "LibraryStats")
    @patch.object(counter, "assign_features", return_value=({'mock_feat'}, 1))
    @patch.object(counter.queue, "Queue", autospec=True)
    @patch.object(counter.parser, "Alignment")
    def test_count_reads_generic(self, aln_obj, queue_obj, assign_features, stats_obj):
        library = {'File': self.short_sam_file, 'Name': 'short'}
        alignment = next(read_SAM(library['File']))     # Grab first alignment object from the SAM reader
        aln_obj.configure_mock(return_value=alignment)  # Patch parser.Alignment to always return same object
        bundle = stats_obj().count_bundle.return_value

        expected_calls_to_stats = [
            call(),  # Produced above when grabbing the return value of count_bundle
            call(library, None, False),  # Call to constructor
            call().count_bundle([alignment]),
            call().count_bundle_alignments(bundle, alignment, {'mock_feat'}),
            call().finalize_bundle(bundle)
        ]

        # CALL FUNCTION
        counter.count_reads(library, queue_obj)

        assign_features.assert_called_once()
        stats_obj.assert_has_calls(expected_calls_to_stats)
        queue_obj.put.assert_called_with(stats_obj())

    # Todo: write factory functions for in-memory GFF and SAM file testing rather than tons of resource files


if __name__ == '__main__':
    unittest.main()
