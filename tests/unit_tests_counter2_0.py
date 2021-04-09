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

        self.short_gff_file = "./testdata/counter/short.gff3"
        self.short_gff = read(self.short_gff_file)

        self.short_sam_file = "./testdata/counter/short.sam"
        self.short_sam = read(self.short_sam_file)

        self.strand = {'sense': tuple('+'), 'antisense': tuple('-'), 'both': ('+', '-')}

        # ID, Key, Value, Strand, Source, Hierarchy, 5pnt, Length
        self.csv_row_dict = {'id_attr': "Alias", "at_key": "Class", "at_val": "CSR", 'strand': "antisense",
                             'gff': "./testdata/cel_ws279/c_elegans.PRJNA13758.WS279.chr1.gff3", 'rank': "1",
                             'nt': '"C,G,U"', "length": "all"}



    # === HELPERS ===

    def get_gff_attr_string(self, file_content):
        return file_content.split('\t')[-1]

    def get_gff_attr_parsed(self, file_content):
        attr_str = self.get_gff_attr_string(file_content)
        return counter.parse_GFF_attribute_string(attr_str)

    def features_csv(self, rules) -> str:
        header = "\uFEFFID Attribute,Attribute Key,Attribute Value,Strand (sense/antisense/both),Feature Source,Hierarchy,5' End Nucleotide,Length"
        return '\n'.join([header, *map(','.join, rules)])
    
    def feat_csv_test_row(self):
        return ','.join(self.csv_row_dict.values())
        
    # === TESTS ===

    """Were only the correct attribute keys present in the parser result?"""

    def test_gff_attr_keys(self):
        attr = self.get_gff_attr_parsed(self.short_gff)
        expected_keys = ["ID", "Name", "interpolated_map_position", "sequence_name",
                         "biotype", "so_term_name", "curie", "Alias", "Class"]

        self.assertEqual(list(attr.keys()), expected_keys)

    """Were keys and values of the correct type in the parser result?"""

    def test_gff_attr_kv_types(self):
        attr = self.get_gff_attr_parsed(self.short_gff)

        # All attribute values should be tuples, all attribute keys should be strs
        self.assertTrue(all([type(val) == tuple for val in attr.values()]))
        self.assertTrue(all([type(key) == str for key in attr.keys()]))

    """Were list values, and non-list values, parsed as such?"""

    def test_gff_attr_list_vals(self):
        attr = self.get_gff_attr_parsed(self.short_gff)

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
        expected_return = {(r['gff'], r['id_attr'])}
        expected_ruleset = [{'Strand': self.strand[r['strand']], 'Hierarchy': int(r['rank']), '5pnt': r['nt'].strip('"'),
                             'Length': r['length'], 'Identity': (r['at_key'], r['at_val'])}]

        self.assertEqual(expected_return, gff_files)
        self.assertEqual(expected_ruleset, ruleset)

    """Does load_config correctly handle duplicate rules? Only one rule should be returned (no duplicates)."""

    def test_load_config_duplicate_rules(self):
        # Features CSV with two duplicate rules/rows
        row = self.csv_row_dict.values()
        csv = self.features_csv([row, row])
        
        with patch('aquatx.srna.counter2_0.open', mock_open(read_data=csv)):
            ruleset, gff_files = counter.load_config('/dev/null')

        r = self.csv_row_dict
        expected_return = {(r['gff'], r['id_attr'])}
        expected_ruleset = [
            {'Strand': self.strand[r['strand']], 'Hierarchy': int(r['rank']), '5pnt': r['nt'].strip('"'),
             'Length': r['length'], 'Identity': (r['at_key'], r['at_val'])}]

        self.assertEqual(expected_return, gff_files)
        self.assertEqual(expected_ruleset, ruleset)

    """DRAFT"""
    def test_heavy(self):
        sam = "./testdata/counter/identity_choice_test.sam"
        gff = "./testdata/counter/identity_choice_test.gff3"
        rules = [["Alias", "Class", "CSR", "antisense", gff, "1", "all", "all"],
                 ["ID", "Class", "piRNA", "sense", gff, "2", "all", "all"]]

        csv = self.features_csv(rules)
        cmd = f"counter -i {sam} -c /dev/null -o test".split(" ")

        with patch("aquatx.srna.counter2_0.open", mock_open(read_data=csv)):
            with patch("sys.argv", cmd):
                counter.main()




    def test_ref_tables_(self):
        pass

    def test_BAM_reader(self):
        samfile = "./testdata/counter/short.sam"
        read = HTSeq.BAM_Reader(samfile)

        for rec in read:
            print(rec)


if __name__ == '__main__':
    unittest.main()
