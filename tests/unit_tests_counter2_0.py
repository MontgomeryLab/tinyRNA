import unittest
from typing import List
from unittest.mock import patch, mock_open, call, Mock

import HTSeq

import aquatx.srna.counter2_0 as counter

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
        self.csv_row_dict = {'id_attr': "Alias", 'at_key': "Class", 'at_val': "CSR", 'rank': "1", 'strand': "antisense",
                             'nt5': '"C,G,U"', 'length': "all", 'match': "Partial",
                             'gff': "./testdata/cel_ws279/c_elegans.PRJNA13758.WS279.chr1.gff3"}


    # === HELPERS ===

    def get_gff_attr_string(self, file_content):
        return file_content.split('\t')[-1]

    def parse_gff_attr(self, gff_file_content):
        attr_str = self.get_gff_attr_string(gff_file_content)
        return counter.parse_GFF_attribute_string(attr_str)

    def features_csv(self, rules) -> str:
        header = "\uFEFFID Attribute,Attribute Key,Attribute Value,Hierarchy,Strand (sense/antisense/both),5' End Nucleotide,Length,Match,Feature Source"
        return '\n'.join([header, *map(','.join, rules)])
    
    def feat_csv_test_row(self):
        return ','.join(self.csv_row_dict.values())
        
    # === TESTS ===

    """Were only the correct attribute keys present in the parser result?"""

    def test_gff_attr_keys(self):
        attr = self.parse_gff_attr(self.short_gff)
        expected_keys = ["ID", "Name", "interpolated_map_position", "sequence_name",
                         "biotype", "so_term_name", "curie", "Alias", "Class"]

        self.assertEqual(list(attr.keys()), expected_keys)

    """Were keys and values of the correct type in the parser result?"""

    def test_gff_attr_kv_types(self):
        attr = self.parse_gff_attr(self.short_gff)

        # All attribute values should be tuples, all attribute keys should be strs
        self.assertTrue(all([type(val) == tuple for val in attr.values()]))
        self.assertTrue(all([type(key) == str for key in attr.keys()]))

    """Were list values, and non-list values, parsed as such?"""

    def test_gff_attr_list_vals(self):
        attr = self.parse_gff_attr(self.short_gff)

        # Comma-separated list values are supported for all attributes, not just the Class attribute
        self.assertTrue(all([len(val) == 1 for key, val in attr.items() if key != "Class"]))
        self.assertEqual(len(attr['Class']), 2)

    """Does load_config correctly parse a single-entry features.csv config file?"""

    def test_load_config_single(self):
        # Features CSV with a single rule/row
        row = self.csv_row_dict.values()
        csv = self.features_csv([row])

        with patch('aquatx.srna.counter2_0.open', mock_open(read_data=csv)):
            ruleset, gff_files = counter.load_config('/dev/null')

        r = self.csv_row_dict
        expected_gff_ret = {(r['gff'], r['id_attr'])}
        expected_ruleset = [{'Strand': self.strand[r['strand']], 'Hierarchy': int(r['rank']), 'nt5': 'C,G,T',
                             'Length': r['length'], 'Identity': (r['at_key'], r['at_val']), 'Strict': r['match'] == 'Full'}]

        self.assertEqual(expected_gff_ret, gff_files)
        self.assertEqual(expected_ruleset, ruleset)

    """Does load_config correctly handle duplicate rules? Only one rule should be returned (no duplicates)."""

    def test_load_config_duplicate_rules(self):
        # Features CSV with two duplicate rules/rows
        row = self.csv_row_dict.values()
        csv = self.features_csv([row, row])  # Duplicate rows
        
        with patch('aquatx.srna.counter2_0.open', mock_open(read_data=csv)):
            ruleset, gff_files = counter.load_config('/dev/null')

        r = self.csv_row_dict
        expected_gff_ret = {(r['gff'], r['id_attr'])}
        expected_ruleset = [
            {'Strand': self.strand[r['strand']], 'Hierarchy': int(r['rank']), 'nt5': 'C,G,T',
             'Length': r['length'], 'Identity': (r['at_key'], r['at_val']), 'Strict': r['match'] == 'Full'}]

        self.assertEqual(expected_gff_ret, gff_files)
        self.assertEqual(expected_ruleset, ruleset)

    """Does load_config convert uracil to thymine for proper matching with cDNA sequences?"""

    def test_load_config_rna_to_cDNA(self):
        row = self.csv_row_dict.copy()
        row['nt5'] = 'U'
        csv = self.features_csv([row.values()])

        with patch('aquatx.srna.counter2_0.open', mock_open(read_data=csv)):
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


    """Does build_reference_tables return the expected features, attributes, and alias for a single record GFF?"""

    def test_ref_tables_single_feature(self):
        feature_source = {(self.short_gff_file, "sequence_name")}
        iv = HTSeq.GenomicInterval("I", 3746, 3908, "-")
        selection_rules = [
            {'Identity': ("Class", "CSR"), 'Strand': "N/A", 'Hierarchy': "N/A", '5pnt': "N/A", 'Length': "N/A", 'Strict': "N/A"},
            {'Identity': ("biotype", "snoRNA"), 'Strand': "N/A", 'Hierarchy': "N/A", '5pnt': "N/A", 'Length': "N/A", 'Strict': "N/A"}
        ]

        feats, attrs, alias = counter.build_reference_tables(feature_source, selection_rules)
        steps = list(feats[iv].array[iv.start:iv.end].get_steps(values_only=True))

        self.assertEqual((type(feats), type(attrs), type(alias)), (HTSeq.GenomicArrayOfSets, dict, dict))
        self.assertEqual(steps, [{"Gene:WBGene00023193"}])
        self.assertEqual(attrs, {'Gene:WBGene00023193': [('Class', ("unknown", "additional_class")), ('biotype', ("snoRNA",))]})
        self.assertEqual(alias, {'Gene:WBGene00023193': ('Y74C9A.6',)})

    """Does build_reference_tables raise ValueError when an ID Attribute refers to a missing attribute?"""

    def test_ref_tables_missing_ID_attribute(self):
        bad = "bad_ID_attribute"
        feature_source = {(self.short_gff_file, bad)}
        selection_rules = []

        expected_err = f"Feature Gene:WBGene00023193 does not contain a '{bad}' attribute in {self.short_gff_file}"

        with self.assertRaisesRegex(ValueError, expected_err):
            counter.build_reference_tables(feature_source, selection_rules)

    """Does build_reference_tables raise ValueError when a selection rule refers to a missing attribute?"""

    def test_ref_tables_missing_identity(self):
        bad = "BAD_attribute_key"
        feature_source = {(self.short_gff_file, "ID")}
        selection_rules = [{'Identity': (bad, "BAD_attribute_value")}]

        expected_err = f"Feature Gene:WBGene00023193 does not contain a '{bad}' attribute in {self.short_gff_file}"

        with self.assertRaisesRegex(ValueError, expected_err):
            counter.build_reference_tables(feature_source, selection_rules)

    """Does build_reference_tables raise ValueError when it encounters a GFF entry without strand information?"""

    def test_ref_tables_unstranded(self):
        gff_file = f"{self.res}/unstranded.gff3"
        feature_source = {(gff_file, "ID")}
        selection_rules = []

        expected_err = f"Feature Gene:WBGene00023193 in {gff_file} has no strand information."

        with self.assertRaisesRegex(ValueError, expected_err):
            counter.build_reference_tables(feature_source, selection_rules)

    """Does build_reference_tables properly concatenate aliases if there is more than one alias for a feature?"""

    def test_ref_tables_alias_concat(self):
        feature_source = {(self.short_gff_file, "ID"), (self.short_gff_file, "sequence_name")}
        selection_rules = []

        expected_alias = {"Gene:WBGene00023193": ("Y74C9A.6",)}
        _, _, alias = counter.build_reference_tables(feature_source, selection_rules)

        self.assertDictEqual(expected_alias, alias)

    """Does build_reference_tables properly concatenate attributes if more than one GFF file defines a feature with different attributes?"""

    def test_ref_tables_attr_concat(self):
        feature_source = {(self.short_gff_file, "ID"), (f"{self.res}/single2.gff3", "ID")}
        selection_rules = [{'Identity': ("Name", "N/A"), 'Strand': "N/A", 'Hierarchy': "N/A", '5pnt': "N/A",
                            'Length': "N/A", 'Strict': "N/A"}]

        _, attrs, _ = counter.build_reference_tables(feature_source, selection_rules)
        expected_attrs = {"Gene:WBGene00023193": [('Name', ("WBGene00023193", "WBGene00023193b"))]}
        actual_attrs = attrs['Gene:WBGene00023193'][1][1]

        # The ordering of values for key 'Name' is non-deterministic since iteration of feature_source {set} is not
        #  This doesn't matter for our use, membership matching of ('Name', ("WBGene00023193")) and ('Name', ("WBGene00023193b"))
        #  will be performed against FeatureSelector's ruleset for Identities
        self.assertIn('WBGene00023193', actual_attrs)
        self.assertIn('WBGene00023193b', actual_attrs)
        self.assertEqual(2, len(actual_attrs))

    """Does assign_features """

    """DRAFT (however, this test works as expected.)"""

    # The components of each test:
    #  1. The GFF file to define a feature and its attributes at an interval
    #  2. The SAM file with a sequence alignment that overlaps a defined feature
    #  3. A selection rule (features.csv) which selects for attributes of the feature and/or read

    def test_counter_main(self):
        rules = [["Alias", "Class", "CSR", "1", "antisense", "all", "all", "Full", self.gff_file],
                 ["sequence_name", "Class", "piRNA", "2", "sense", "all", "all", "Partial", self.gff_file]]

        csv = self.features_csv(rules)
        cmd = f"counter -i {self.sam_file} -c /dev/null -o test".split(" ")

        with patch("aquatx.srna.counter2_0.open", mock_open(read_data=csv)):
            with patch("sys.argv", cmd):
                counter.main()

    # Todo: additional test for build_ref_tables: what happens when "Class" is and is not defined in features.csv?
    # Todo: write factory functions for in-memory GFF and SAM file testing rather than tons of resource files


if __name__ == '__main__':
    unittest.main()
