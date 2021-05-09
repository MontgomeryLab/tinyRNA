import unittest
from typing import List
from unittest.mock import patch, mock_open, call, Mock

import HTSeq

import aquatx.srna.counter as counter

class CounterTests(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        # Simply for convenience for loading files during setup
        def read(file, mode='r'):
            with open(file, mode) as f:
                return f.read()

        self.res = "./testdata/counter"

        self.gff_file = f"{self.res}/identity_choice_test.gff3"
        self.short_gff_file = f"{self.res}/single.gff3"
        self.short_gff = read(self.short_gff_file)

        self.sam_file = f"{self.res}/identity_choice_test.sam"
        self.short_sam_file = f"{self.res}/single.sam"
        self.short_sam = read(self.short_sam_file)

        self.strand = {'sense': tuple('+'), 'antisense': tuple('-'), 'both': ('+', '-')}

        # ID, Key, Value, Hierarchy, Strand, 5pnt, Length, Match, Source
        self.csv_feat_row_dict = {'id_attr': "Alias", 'at_key': "Class", 'at_val': "CSR", 'rank': "1",
                                  'strand': "antisense", 'nt5': '"C,G,U"', 'length': "all", 'match': "Partial",
                                  'gff': "./testdata/cel_ws279/c_elegans.PRJNA13758.WS279.chr1.gff3"}

        self.csv_samp_row_dict = {'file': self.short_sam_file, 'group': "test_group", 'rep': "0"}

    # === HELPERS ===

    def csv(self, type, rows):
        header = "\uFEFF"
        if type == "features.csv":
            header = "\uFEFFID Attribute,Attribute Key,Attribute Value,Hierarchy,Strand (sense/antisense/both),5' End Nucleotide,Length,Match,Feature Source"
        elif type == "samples.csv":
            header = "\uFEFFInput FastQ/A Files,Sample/Group Name,Replicate number"

        return '\n'.join([header, *map(','.join, rows)])
    
    def feat_csv_test_row(self):
        return ','.join(self.csv_feat_row_dict.values())
        
    # === TESTS ===

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

    """Does load_samples correctly parse a single record samples.csv?"""

    def test_load_samples_single(self):
        row = self.csv_samp_row_dict
        csv = self.csv("samples.csv", [row.values()])

        with patch('aquatx.srna.counter.open', mock_open(read_data=csv)):
            inputs = counter.load_samples('/dev/null')

        expected_lib_name = f"{row['group']}_replicate_{row['rep']}"
        expected_result = [(self.short_sam_file, expected_lib_name)]
        self.assertEqual(expected_result, inputs)


    """Does load_samples throw ValueError if a non-absolute path to a SAM file is provided?"""

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




    """Does assign_features """

    """DRAFT (however, this test works as expected.)"""

    # The components of each test:
    #  1. The GFF file to define a feature and its attributes at an interval
    #  2. The SAM file with a sequence alignment that overlaps a defined feature
    #  3. A selection rule (features.csv) which selects for attributes of the feature and/or read

    def test_counter_main(self):
        rules = [["Alias", "Class", "CSR", "1", "antisense", "all", "all", "Full", self.gff_file],
                 ["sequence_name", "Class", "piRNA", "2", "sense", "all", "all", "Partial", self.gff_file]]

        csv = self.csv("features.csv", rules)
        cmd = f"counter -i {self.sam_file} -c /dev/null -o test".split(" ")

        with patch("aquatx.srna.counter.open", mock_open(read_data=csv)):
            with patch("sys.argv", cmd):
                counter.main()

    # Todo: additional test for build_ref_tables: what happens when "Class" is and is not defined in features.csv?
    # Todo: write factory functions for in-memory GFF and SAM file testing rather than tons of resource files


if __name__ == '__main__':
    unittest.main()
