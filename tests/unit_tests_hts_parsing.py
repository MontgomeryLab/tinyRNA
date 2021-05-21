import unittest

from aquatx.srna.hts_parsing import *
from tests.unit_tests_counter import resources
import tests.unit_test_helpers as helpers


class MyTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.gff_file = f"{resources}/identity_choice_test.gff3"
        self.short_gff_file = f"{resources}/single.gff3"
        self.short_gff = helpers.read(self.short_gff_file)

        self.sam_file = f"{resources}/identity_choice_test.sam"
        self.short_sam_file = f"{resources}/single.sam"
        self.short_sam = helpers.read(self.short_sam_file)

    # === HELPERS ===

    def get_gff_attr_string(self, file_content):
        return file_content.split('\t')[-1]

    def parse_gff_attr(self, gff_file_content):
        attr_str = self.get_gff_attr_string(gff_file_content)
        return parse_GFF_attribute_string(attr_str)

    # === TESTS ===

    """Did SAM_reader correctly skip header values and parse all pertinent info from a single record SAM file?"""

    def test_sam_reader(self):
        sam_record = next(read_SAM(self.short_sam_file))

        self.assertEqual(sam_record.iv, HTSeq.GenomicInterval("I", 15064569, 15064590, '-'))
        self.assertEqual(sam_record.read.name, "read_id")
        self.assertEqual(sam_record.read.seq, b"CAAGACAGAGCTTCACCGTTC")
        self.assertEqual(sam_record.read.len, 21)

    """Does the alignment object construct and retain expected attributes and structure?"""

    def test_alignment_obj(self):
        iv = HTSeq.GenomicInterval("I", 15064569, 15064590, "-")
        seq = b"CAAGACAGAGCTTCACCGTTC"
        name = "test_aln"

        aln = Alignment(iv, name, seq)

        # The following object structure is expected by HTSeq, StatsCollector, and FeatureSelector
        self.assertEqual(aln.iv, iv)
        self.assertEqual(aln.iv.strand, "-")
        self.assertEqual(aln.read.seq, seq)
        self.assertEqual(aln.read.name, name)
        self.assertEqual(aln.read.len, len(seq))

    """Does our custom SAM parser produce the same pertinent info as HTSeq's BAM_reader?
    
    A note on SAM files: reads are always stored 5' to 3', so antisense reads are actually
    recorded in reverse complement. HTSeq automatically performs this conversion, but we
    are only really concerned about a sequence's 5' end NT, so our Alignment class performs
    this conversion more surgically for only the 5' end NT at construction time.
    """

    # Todo: move all referenced .sam files to a more appropriate (and available) testdata folder
    def test_sam_parser_comparison(self):
        file = "./run_directory/Lib304_test_aligned_seqs.sam"
        ours = read_SAM(file)
        theirs = HTSeq.BAM_Reader(file)

        for our, their in zip(ours, theirs):
            self.assertEqual(our.iv, their.iv)
            self.assertEqual(our.iv.strand, their.iv.strand)
            self.assertEqual(our.read.nt5, chr(their.read.seq[0]))  # See note above
            self.assertEqual(our.read.name, their.read.name)
            self.assertEqual(len(our.read), len(their.read))

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

    """Does build_reference_tables return the expected features, attributes, and alias for a single record GFF?"""

    def test_ref_tables_single_feature(self):
        feature_source = {(self.short_gff_file, "sequence_name")}
        iv = HTSeq.GenomicInterval("I", 3746, 3908, "-")
        selection_rules = [
            {'Identity': ("Class", "CSR"), 'Strand': "N/A", 'Hierarchy': "N/A", '5pnt': "N/A", 'Length': "N/A",
             'Strict': "N/A"},
            {'Identity': ("biotype", "snoRNA"), 'Strand': "N/A", 'Hierarchy': "N/A", '5pnt': "N/A", 'Length': "N/A",
             'Strict': "N/A"}
        ]

        feats, attrs, alias = build_reference_tables(feature_source, selection_rules)
        steps = list(feats[iv].array[iv.start:iv.end].get_steps(values_only=True))

        self.assertEqual((type(feats), type(attrs), type(alias)), (HTSeq.GenomicArrayOfSets, dict, dict))
        self.assertEqual(steps, [{"Gene:WBGene00023193"}])
        self.assertEqual(attrs, {
            'Gene:WBGene00023193': [('Class', ("unknown", "additional_class")), ('biotype', ("snoRNA",))]})
        self.assertEqual(alias, {'Gene:WBGene00023193': ('Y74C9A.6',)})

    """Does build_reference_tables raise ValueError when an ID Attribute refers to a missing attribute?"""

    def test_ref_tables_missing_ID_attribute(self):
        bad = "bad_ID_attribute"
        feature_source = {(self.short_gff_file, bad)}
        selection_rules = []

        expected_err = f"Feature Gene:WBGene00023193 does not contain a '{bad}' attribute in {self.short_gff_file}"

        with self.assertRaisesRegex(ValueError, expected_err):
            build_reference_tables(feature_source, selection_rules)

    """Does build_reference_tables raise ValueError when a selection rule refers to a missing attribute?"""

    def test_ref_tables_missing_identity(self):
        bad = "BAD_attribute_key"
        feature_source = {(self.short_gff_file, "ID")}
        selection_rules = [{'Identity': (bad, "BAD_attribute_value")}]

        expected_err = f"Feature Gene:WBGene00023193 does not contain a '{bad}' attribute in {self.short_gff_file}"

        with self.assertRaisesRegex(ValueError, expected_err):
            build_reference_tables(feature_source, selection_rules)

    """Does build_reference_tables raise ValueError when it encounters a GFF entry without strand information?"""

    def test_ref_tables_unstranded(self):
        gff_file = f"{resources}/unstranded.gff3"
        feature_source = {(gff_file, "ID")}
        selection_rules = []

        expected_err = f"Feature Gene:WBGene00023193 in {gff_file} has no strand information."

        with self.assertRaisesRegex(ValueError, expected_err):
            build_reference_tables(feature_source, selection_rules)

    """Does build_reference_tables properly concatenate aliases if there is more than one alias for a feature?"""

    def test_ref_tables_alias_concat(self):
        feature_source = {(self.short_gff_file, "ID"), (self.short_gff_file, "sequence_name")}
        selection_rules = []

        expected_alias = {"Gene:WBGene00023193": ("Y74C9A.6",)}
        _, _, alias = build_reference_tables(feature_source, selection_rules)

        self.assertDictEqual(expected_alias, alias)

    """Does build_reference_tables properly concatenate attributes if more than one GFF file defines a feature with different attributes?"""

    def test_ref_tables_attr_concat(self):
        feature_source = {(self.short_gff_file, "ID"), (f"{resources}/single2.gff3", "ID")}
        selection_rules = [{'Identity': ("Name", "N/A"), 'Strand': "N/A", 'Hierarchy': "N/A", '5pnt': "N/A",
                            'Length': "N/A", 'Strict': "N/A"}]

        _, attrs, _ = build_reference_tables(feature_source, selection_rules)
        expected_attrs = {"Gene:WBGene00023193": [('Name', ("WBGene00023193", "WBGene00023193b"))]}
        actual_attrs = attrs['Gene:WBGene00023193'][1][1]

        # The ordering of values for key 'Name' is non-deterministic since iteration of feature_source {set} is not
        #  This doesn't matter for our use, membership matching of ('Name', ("WBGene00023193")) and ('Name', ("WBGene00023193b"))
        #  will be performed against FeatureSelector's ruleset for Identities
        self.assertIn('WBGene00023193', actual_attrs)
        self.assertIn('WBGene00023193b', actual_attrs)
        self.assertEqual(2, len(actual_attrs))

    # Todo: additional test for build_ref_tables: what happens when "Class" is and is not defined in features.csv?

if __name__ == '__main__':
    unittest.main()
