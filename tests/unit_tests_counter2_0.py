import unittest
from typing import List
from unittest.mock import patch, mock_open, call, Mock

import HTSeq

import aquatx.srna.counter2_0 as counter

class MyTestCase(unittest.TestCase):
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


    # === HELPERS ===

    def get_gff_attr_string(self, file_content):
        return file_content.split('\t')[-1]

    def get_gff_attr_parsed(self, file_content):
        attr_str = self.get_gff_attr_string(file_content)
        return counter.parse_GFF_attribute_string(attr_str)

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

    """Were list values, and non-list values, properly parsed?"""

    def test_gff_attr_list_vals(self):
        attr = self.get_gff_attr_parsed(self.short_gff)

        # Comma-separated list values are supported for all attributes, not just the Class attribute
        self.assertTrue(all([len(val) == 1 for key, val in attr.items() if key != "Class"]))
        self.assertEqual(len(attr['Class']), 2)

    def features_csv(self, rules: List[list]) -> str:
        header = "\uFEFFID Attribute,Attribute Key,Attribute Value,Strand (sense/antisense/both),Feature Source,Hierarchy,5' End Nucleotide,Length"
        return '\n'.join([header, *map(','.join, rules)])

    """Does load_config correctly parse a single-entry features.csv config file?"""

    @patch('aquatx.srna.counter2_0.FeatureSelector')
    def test_load_config_single(self, feat_select):
        # Features CSV with a single rule/row
        id_attr, at_key, at_val = "Alias", "Class", "CSR"
        strand, rank, nt, length = "antisense", "1", '"C,G,U"', "all"
        gff = './testdata/cel_ws279/c_elegans.PRJNA13758.WS279.chr1.gff3'
        csv = self.features_csv([[id_attr, at_key, at_val, strand, gff, rank, nt, length]])

        with patch('aquatx.srna.counter2_0.open', mock_open(read_data=csv)):
            gff_files = counter.load_config('/dev/null')

        expected_return = {(gff, id_attr)}
        expected_ruleset = [{'Strand': self.strand[strand], 'Hierarchy': int(rank), '5pnt': nt.strip('"'), 'Length': length, 'Identity': (at_key, at_val)}]

        self.assertEqual(expected_return, gff_files)
        self.assertEqual(call(expected_ruleset), feat_select.call_args)



    def test_BAM_reader(self):
        samfile = "./testdata/counter/short.sam"
        read = HTSeq.BAM_Reader(samfile)

        for rec in read:
            print(rec)


if __name__ == '__main__':
    unittest.main()
