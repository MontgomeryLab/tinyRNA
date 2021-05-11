import multiprocessing as mp
import unittest
import HTSeq
import queue

from unittest.mock import patch, mock_open, call, MagicMock

import tests.unit_test_helpers as helpers
import aquatx.srna.counter as counter
from aquatx.srna.hts_parsing import Alignment, read_SAM

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

        # ID, Key, Value, Hierarchy, Strand, 5pnt, Length, Match, Source
        self.csv_feat_row_dict = {'id_attr': "Alias", 'at_key': "Class", 'at_val': "CSR", 'rank': "1",
                                  'strand': "antisense", 'nt5': '"C,G,U"', 'length': "all", 'match': "Partial",
                                  'gff': "./testdata/cel_ws279/c_elegans.PRJNA13758.WS279.chr1.gff3"}

        self.csv_samp_row_dict = {'file': "test_file.fastq", 'group': "test_group", 'rep': "0"}

    # === HELPERS ===

    @staticmethod
    def csv(type, rows):
        header = "\uFEFF"
        if type == "features.csv":
            header = "\uFEFFID Attribute,Attribute Key,Attribute Value,Hierarchy,Strand (sense/antisense/both),5' End Nucleotide,Length,Match,Feature Source"
        elif type == "samples.csv":
            header = "\uFEFFInput FastQ/A Files,Sample/Group Name,Replicate number"

        return '\n'.join([header, *map(','.join, rows)])
    
    def feat_csv_test_row(self):
        return ','.join(self.csv_feat_row_dict.values())
        
    # === TESTS ===

    """Does load_samples correctly parse a single record samples.csv?"""

    def test_load_samples_single(self):
        inp_file = "test.fastq"
        exp_file = "test_aligned_seqs.sam"

        row = {'File': inp_file, 'Group': "test_group", 'Rep': "0"}
        csv = self.csv("samples.csv", [row.values()])

        with patch('aquatx.srna.counter.open', mock_open(read_data=csv)):
            inputs = counter.load_samples('/dev/null')

        expected_lib_name = f"{row['Group']}_replicate_{row['Rep']}"
        expected_result = [{'File': exp_file, 'Name': expected_lib_name}]
        self.assertEqual(expected_result, inputs)

    """Does load_samples correctly handle duplicate samples? There should be no duplicates."""

    def test_load_samples_duplicate(self):
        row = {'File': "test.fastq", 'Group': "N/A", 'Rep': "N/A"}
        csv = self.csv("samples.csv", [row.values(), row.values()])

        with patch('aquatx.srna.counter.open', mock_open(read_data=csv)):
            inputs = counter.load_samples('/dev/null')

        self.assertEqual(1, len(inputs))

    """Does load_samples correctly handle SAM filenames?"""

    def test_load_samples_sam(self):
        sam_filename = "/fake/absolute/path/sample.sam"
        row = {'File': sam_filename, 'Group': "test_group", 'Rep': "0"}
        csv = self.csv("samples.csv", [row.values()])

        with patch('aquatx.srna.counter.open', mock_open(read_data=csv)):
            inputs = counter.load_samples('/dev/null')

        expected_lib_name = f"{row['Group']}_replicate_{row['Rep']}"
        expected_result = [{'File': sam_filename, 'Name': expected_lib_name}]

        self.assertEqual(expected_result, inputs)

    """Does load_samples throw ValueError if a non-absolute path to a SAM file is provided?"""

    def test_load_samples_nonabs_path(self):
        bad = "./dne.sam"
        row = [bad, "test_group", "0"]
        csv = self.csv("samples.csv", [row])

        expected_error = "The following file must be expressed as an absolute path:\n" + bad

        with patch('aquatx.srna.counter.open', mock_open(read_data=csv)):
            with self.assertRaisesRegex(ValueError, expected_error):
                counter.load_samples('/dev/null')

    """Does load_samples throw ValueError if sample filename does not have a .fastq or .sam extension?"""

    def test_load_samples_bad_extension(self):
        bad = "./bad_extension.xyz"
        row = [bad, "test_group", "0"]
        csv = self.csv("samples.csv", [row])

        expected_error = "The filenames defined in your samples CSV file must have a .fastq or .sam extension.\n" \
                         "The following filename contained neither:\n" + bad

        with patch('aquatx.srna.counter.open', mock_open(read_data=csv)):
            with self.assertRaisesRegex(ValueError, expected_error):
                counter.load_samples('/dev/null')

    """Does load_config correctly parse a single-entry features.csv config file?"""

    def test_load_config_single(self):
        # Features CSV with a single rule/row
        row = self.csv_feat_row_dict.values()
        csv = self.csv("features.csv", [row])

        with patch('aquatx.srna.counter.open', mock_open(read_data=csv)):
            ruleset, gff_files = counter.load_config('/dev/null')

        r = self.csv_feat_row_dict
        expected_gff_ret = {(r['gff'], r['id_attr'])}
        expected_ruleset = [{'Strand': self.strand[r['strand']], 'Hierarchy': int(r['rank']), 'nt5': 'C,G,T',
                             'Length': r['length'], 'Identity': (r['at_key'], r['at_val']), 'Strict': r['match'] == 'Full'}]

        self.assertEqual(expected_gff_ret, gff_files)
        self.assertEqual(expected_ruleset, ruleset)

    """Does load_config correctly handle duplicate rules? Only one rule should be returned (no duplicates)."""

    def test_load_config_duplicate_rules(self):
        # Features CSV with two duplicate rules/rows
        row = self.csv_feat_row_dict.values()
        csv = self.csv("features.csv", [row, row])  # Duplicate rows
        
        with patch('aquatx.srna.counter.open', mock_open(read_data=csv)):
            ruleset, gff_files = counter.load_config('/dev/null')

        r = self.csv_feat_row_dict
        expected_gff_ret = {(r['gff'], r['id_attr'])}
        expected_ruleset = [
            {'Strand': self.strand[r['strand']], 'Hierarchy': int(r['rank']), 'nt5': 'C,G,T',
             'Length': r['length'], 'Identity': (r['at_key'], r['at_val']), 'Strict': r['match'] == 'Full'}]

        self.assertEqual(expected_gff_ret, gff_files)
        self.assertEqual(expected_ruleset, ruleset)

    """Does load_config convert uracil to thymine for proper matching with cDNA sequences?"""

    def test_load_config_rna_to_cDNA(self):
        row = self.csv_feat_row_dict.copy()
        row['nt5'] = 'U'
        csv = self.csv("features.csv", [row.values()])

        with patch('aquatx.srna.counter.open', mock_open(read_data=csv)):
            ruleset, _ = counter.load_config('/dev/null')

        self.assertEqual(ruleset[0]['nt5'], 'T')

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

    """Does count_reads properly handle a single record library?"""

    @patch.object(counter, "LibraryStats", autospec=True)
    @patch.object(counter, "assign_features", return_value=({'mock_feat'}, 1))
    def test_count_reads_generic(self, assign_features, LibraryStats):
        library = {'File': self.short_sam_file, 'Name': 'short'}
        alignment = next(read_SAM(library['File']))
        ret_q = queue.Queue()

        expected_LibStat_count_bundle_call = call().count_bundle([alignment])
        bundle = LibraryStats.count_bundle.return_value = MagicMock()
        # LibraryStats.count_bundle = MagicMock()
        # LibraryStats.count_bundle_alignments = MagicMock()
        # LibraryStats.finalize_bundle = MagicMock

        expected_LibStats_calls = [
            call(library, None, False),
            expected_LibStat_count_bundle_call,
            call().count_bundle_alignments(bundle),
            call().finalize_bundle(bundle)
        ]

        counter.count_reads(library, ret_q)

        assign_features.assert_called_once()
        LibraryStats.assert_has_calls(expected_LibStats_calls)

    # Todo: write factory functions for in-memory GFF and SAM file testing rather than tons of resource files


if __name__ == '__main__':
    unittest.main()
