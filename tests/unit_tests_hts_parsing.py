import unittest
from unittest.mock import patch

from rna.counter.features import FeatureSelector
from tiny.rna.counter.hts_parsing import *
# from tests.unit_tests_counter import resources
import unit_test_helpers as helpers
resources = "./testdata/counter"


class MockFeatureSelector:
    def __init__(self, rules_table):
        self.rules_table = FeatureSelector.build_selectors(rules_table)
        self.inv_ident = FeatureSelector.build_inverted_identities(rules_table)


class MyTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.gff_file = f"{resources}/identity_choice_test.gff3"
        self.short_gff_file = f"{resources}/single.gff3"
        self.short_gff = helpers.read(self.short_gff_file)

        self.sam_file = f"{resources}/identity_choice_test.sam"
        self.short_sam_file = f"{resources}/single.sam"
        self.short_sam = helpers.read(self.short_sam_file)

        self.rules_template = [{'Identity': ("Name", "N/A"), 'Strand': "+", 'Hierarchy': "0", 'nt5end': "N",
                                'Length': "0", 'Strict': True}]

    # === HELPERS ===

    def get_gff_attr_string(self, gff_line):
        return gff_line.split('\t')[-1]

    def parse_gff_attr(self, gff_file_content):
        attr_str = self.get_gff_attr_string(gff_file_content)
        return parse_GFF_attribute_string(attr_str)

    def make_single_sam(self, name="read_id", flag="16", chrom="I", pos="15064570", seq="CAAGACAGAGCTTCACCGTTC"):
        length = str(len(seq))
        header = '\t'.join(["@SQ", "SN:%s", "LN:%s"]) % (chrom, length)
        record = '\t'.join([
            name, flag, chrom, pos, "255", length + "M", "*", "0", "0", seq,
            "IIIIIIIIIIIIIIIIIIIII", "XA:i:0",	"MD:Z:" + length, "NM:i:0", "XM:i:2"])

        return header + '\n' + record

    # === TESTS ===

    def test_cystep(self):
        import CyStep._stepvector as StepVector
        setattr(HTSeq._HTSeq, "StepVector", StepVector)

        # gas = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        gas = HTSeq.GenomicArray(chroms="auto", stranded=True, typecode='O')
        iva = HTSeq.GenomicInterval("I", 1, 10, '-')
        ivb = HTSeq.GenomicInterval("I", 5, 15, '-')
        ivc = HTSeq.GenomicInterval("I", 9, 20, '-')
        ivd = HTSeq.GenomicInterval("I", 2, 4, "-")
        gas[iva] += {"TestA"}
        gas[ivb] += {"TestB"}
        gas[ivd] += {"TestD"}

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

        # list(gas[ivc].array.get_steps())
        # Mine:  [(0, 2, {'TestA'}), (2, 4, {'TestA', 'TestD'}), (4, 5, {'TestA'}), (5, 10, {'TestB', 'TestA'}), (10, 9223372036854775807, {'TestB'})]
        # Their: [(0, 2, {'TestA'}), (2, 4, {'TestA', 'TestD'}), (4, 5, {'TestA'}), (5, 10, {'TestA', 'TestB'}), (10, 9223372036854775807, {'TestB'})]

        # list(gas[ivc].array[ivc.start:ivc.end].get_steps())
        # Mine:  [(9, 10, {'TestA', 'TestB'}), (10, 15, {'TestB'}), (15, 20, set())]
        # Their: [(9, 10, {'TestB', 'TestA'}), (10, 15, {'TestB'}), (15, 20, set())]

        matches = list(gas[ivc].array[ivc.start:ivc.end].get_steps(values_only=True))
        matches_with_cooridnates = list(gas[ivc].steps())
        self.assertEqual(matches, [{"TestA", "TestB"}, {"TestB"}, set()])
        self.assertEqual(matches_with_cooridnates[0][0], HTSeq.GenomicInterval("I", 9, 10, '-'))
        self.assertEqual(matches_with_cooridnates[1][0], HTSeq.GenomicInterval("I", 10, 15, '-'))
        self.assertEqual(matches_with_cooridnates[2][0], HTSeq.GenomicInterval("I", 15, 20, '-'))
        self.assertEqual(matches_with_cooridnates[2][0], HTSeq.GenomicInterval("I", 15, 20, '-'))


    """Did SAM_reader correctly skip header values and parse all pertinent info from a single record SAM file?"""

    def test_sam_reader(self):
        sam_bundle = next(read_SAM(self.short_sam_file))
        sam_record = sam_bundle[0]

        self.assertEqual(sam_record['chrom'], "I")
        self.assertEqual(sam_record['start'], 15064569)
        self.assertEqual(sam_record['end'], 15064590)
        self.assertEqual(sam_record['strand'], '-')
        self.assertEqual(sam_record['name'], "read_id")
        self.assertEqual(sam_record['seq'], b"CAAGACAGAGCTTCACCGTTC")
        self.assertEqual(sam_record['len'], 21)
        self.assertEqual(sam_record['nt5'], 'G')

    """Does our custom SAM parser produce the same pertinent info as HTSeq's BAM_reader?
    
    A note on SAM files: reads are always stored 5' to 3', so antisense reads are actually
    recorded in reverse complement. HTSeq automatically performs this conversion, but we
    are only really concerned about a sequence's 5' end NT, so our alignment dicts performs
    this conversion more surgically for only the 5' end NT at construction time.
    """

    def test_sam_parser_comparison(self):
        file = f"{resources}/Lib304_test.sam"
        ours = read_SAM(file)
        theirs = HTSeq.bundle_multiple_alignments(HTSeq.BAM_Reader(file))

        for our_bundle, their_bundle in zip(ours, theirs):
            self.assertEqual(len(our_bundle), len(their_bundle))
            for our, their in zip(our_bundle, their_bundle):
                self.assertEqual(our['chrom'], their.iv.chrom)
                self.assertEqual(our['start'], their.iv.start)
                self.assertEqual(our['end'], their.iv.end)
                self.assertEqual(our['strand'], their.iv.strand)
                self.assertEqual(our['nt5'], chr(their.read.seq[0]))  # See note above
                self.assertEqual(our['name'], their.read.name)
                self.assertEqual(our['seq'], their.read.seq)

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

    """Does ReferenceTables.get() return the expected features, identities, aliases, and classes for a single record GFF?"""

    def test_ref_tables_single_feature(self):
        feature_source = {self.short_gff_file: ["sequence_name"]}
        iv = HTSeq.GenomicInterval("I", 3746, 3908, "-")
        mock_selector = MockFeatureSelector([
            {'Identity': ("Class", "CSR"), 'Strand': "+", 'Hierarchy': 1, 'nt5end': "N/A", 'Length': "20",
             'Strict': True},
            {'Identity': ("biotype", "snoRNA"), 'Strand': "-", 'Hierarchy': 2, 'nt5end': "N/A", 'Length': "30",
             'Strict': False}
        ])
        kwargs = {'all_features': True}

        feats, alias, classes = ReferenceTables(feature_source, mock_selector, **kwargs).get()
        steps = list(feats[iv].array[iv.start:iv.end].get_steps(values_only=True))

        self.assertEqual((type(feats), type(alias), type(classes)), (HTSeq.GenomicArrayOfSets, dict, dict))
        self.assertEqual(steps, [{("Gene:WBGene00023193", '-', ((1, 2, False),))}])
        self.assertEqual(alias, {'Gene:WBGene00023193': ('Y74C9A.6',)})
        self.assertEqual(classes, {'Gene:WBGene00023193': ('additional_class', 'unknown')})

    """Does ReferenceTables.get() raise ValueError when a Name Attribute refers to a missing attribute?"""

    def test_ref_tables_missing_name_attribute(self):
        bad = "bad_name_attribute"
        feature_source = {self.short_gff_file: [bad]}
        selection_rules = MockFeatureSelector([])
        kwargs = {'all_features': True}

        expected_err = f"Feature Gene:WBGene00023193 does not contain a '{bad}' attribute." + '\n'
        expected_err += f"Error occurred on line 1 of {self.short_gff_file}"

        with self.assertRaisesRegex(ValueError, expected_err):
            ReferenceTables(feature_source, selection_rules, **kwargs).get()

    """Does ReferenceTables.get() raise ValueError when a selection rule refers to a missing attribute?"""

    def test_ref_tables_missing_identity(self):
        bad = "BAD_attribute_key"
        feature_source = {self.short_gff_file: ["ID"]}
        selection_rules = [{'Identity': (bad, "BAD_attribute_value")}]
        kwargs = {'all_features': True}

        expected_err = f"Feature Gene:WBGene00023193 does not contain a '{bad}' attribute." + '\n'
        expected_err += f"Error occurred on line 1 of {self.short_gff_file}"

        with self.assertRaisesRegex(ValueError, expected_err):
            ReferenceTables(feature_source, selection_rules, **kwargs).get()

    """Does ReferenceTables.get() properly concatenate aliases if there is more than one alias for a feature?"""
    """Does ReferenceTables.get() properly concatenate aliases when Name Attribute refers to a list-type alias?"""
    # 2 for 1!

    def test_ref_tables_alias_concat(self):
        feature_source = {self.short_gff_file: ["ID", "Class"]}
        kwargs = {'all_features': True}

        # Notice: screening for "ID" name attribute happens earlier in counter.load_config()
        expected_alias = {"Gene:WBGene00023193": ("Gene:WBGene00023193", "additional_class", "unknown")}
        _, alias, _ = ReferenceTables(feature_source, MockFeatureSelector([]), **kwargs).get()

        self.assertDictEqual(alias, expected_alias)

    """Does ReferenceTables.get() properly concatenate identities if more than one GFF file defines a feature with different identities?"""
    # todo
    def test_ref_tables_idents_concat(self):
        feature_source = {self.short_gff_file: ["ID"], f"{resources}/single2.gff3": ["ID"]}
        selection_rules = self.rules_template
        kwargs = {'all_features': True}

        _, idents, _, _, _ = ReferenceTables(feature_source, selection_rules, **kwargs).get()

        # Notice: the 'Class' attribute is included first by default for Feature Class column of Feature Counts outfile
        expected_idents = [('Class', ('unknown', 'additional_class')), ('Name', ("WBGene00023193", "WBGene00023193b"))]
        actual_idents = idents['Gene:WBGene00023193']

        for act_attr, exp_attr in zip(actual_idents, expected_idents):
            self.assertEqual(act_attr[0], exp_attr[0])

            # Notice: ordering of attribute values is non-deterministic. This is due to an intermediary set
            # that is used when concatenating attributes to prevent duplicate values.
            act_sorted = sorted(act_attr[1])
            exp_sorted = sorted(exp_attr[1])
            self.assertEqual(act_sorted, exp_sorted)

    """Does ReferenceTables.get() properly handle aliases for discontinuous features?"""

    def test_ref_tables_discontinuous_aliases(self):
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        mock_selector = MockFeatureSelector(self.rules_template)
        kwargs = {'all_features': True}

        _, alias, _ = ReferenceTables(feature_source, mock_selector, **kwargs).get()

        # Ancestor depth of 1, distinct aliases
        self.assertEqual(alias['Parent2'], ('Child2Name', 'Parent2Name'))
        # Ancestor depth >1, shared aliases
        self.assertEqual(alias['GrandParent'], ('SharedName',))
        # Siblings, distinct aliases
        self.assertEqual(alias['Sibling'], ('Sibling1', 'Sibling2', 'Sibling3'))

    """Does ReferenceTables.get() properly handle intervals for discontinous features?"""

    def test_ref_tables_discontinuous_intervals(self):
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        selection_rules = self.rules_template
        kwargs = {'all_features': True}

        _, _, _, intervals, _ = ReferenceTables(feature_source, selection_rules, **kwargs).get()

        grandparent_iv = HTSeq.GenomicInterval('I', 0, 10, '-')
        parent_w_p_iv = HTSeq.GenomicInterval('I', 9, 20, '-')
        child_w_gp_iv = HTSeq.GenomicInterval('I', 29, 40, '-')
        parent_2 = HTSeq.GenomicInterval('I', 19, 30, '-')
        child_2 = HTSeq.GenomicInterval('I', 39, 50, '-')
        sib_1 = HTSeq.GenomicInterval('I', 99, 110, '-')
        sib_2 = HTSeq.GenomicInterval('I', 110, 120, '-')
        sib_3 = HTSeq.GenomicInterval('I', 139, 150, '-')

        # Ancestor depth of 1
        self.assertEqual(intervals['GrandParent'], [grandparent_iv, parent_w_p_iv, child_w_gp_iv])
        # Ancestor depth >1
        self.assertEqual(intervals['Parent2'], [parent_2, child_2])
        # Siblings
        self.assertEqual(intervals['Sibling'], [sib_1, sib_2, sib_3])

    """Does ReferenceTables.get() properly handle attributes for discontinuous features?"""
    # todo
    def test_ref_tables_discontinuous_attributes(self):
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        selection_rules = self.rules_template
        kwargs = {'all_features': True}

        _, attrs, _, _, _ = ReferenceTables(feature_source, selection_rules, **kwargs).get()

        for multival in ['Sibling', 'Parent2']:
            for i, attr in enumerate(attrs[multival]):
                attrs[multival][i] = (attr[0], tuple(sorted(attr[1])))

        self.assertEqual(attrs['GrandParent'], [('Class', ('NA',)), ('Name', ('SharedName',))])
        self.assertEqual(attrs['Parent2'], [('Class', ('NA',)), ('Name', ('Child2Name', 'Parent2Name'))])
        self.assertEqual(attrs['Sibling'], [('Class', ('Class1', 'Class2', 'Class3')), ('Name', ('Sibling1', 'Sibling2', 'Sibling3'))])

    """Does ReferenceTables.get() properly build a GenomicArrayOfSets for discontinuous features?"""

    def test_ref_tables_discontinuous_features(self):
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        selection_rules = self.rules_template
        kwargs = {'all_features': True}

        feats, _, _, _, _ = ReferenceTables(feature_source, selection_rules, **kwargs).get()

        expected = [(0, 19, {'GrandParent'}),
                    (19, 20, {'GrandParent', 'Parent2'}),
                    (20, 29, {'Parent2'}),
                    (29, 30, {'GrandParent', 'Parent2'}),
                    (30, 39, {'GrandParent'}),
                    (39, 40, {'GrandParent', 'Parent2'}),
                    (40, 50, {'Parent2'}),
                    (50, 99, set()),
                    (99, 120, {'Sibling'}),
                    (120, 139, set()),
                    (139, 150, {'Sibling'}),
                    (150, sys.maxsize, set())]

        for act, exp in zip(feats.chrom_vectors["I"]["-"].array.get_steps(), expected):
            self.assertEqual(act, exp)

    """Does ReferenceTables.get() properly handle source filters for discontinuous features?"""
    # todo
    def test_ref_tables_source_filter(self):

        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        selection_rules = self.rules_template
        kwargs = {'all_features': True}

        rt = ReferenceTables(feature_source, selection_rules, source_filter=["Source2Name"], **kwargs)
        feats, idents, alias, intervals, classes = rt.get()

        exp_alias = {'Child2': ('Child2Name',)}
        exp_ident = {'Child2': ('Child2Name', 'Name')}
        exp_feats = [(0, 39, set()), (39, 50, {'Child2'}), (50, sys.maxsize, set())]
        exp_intervals = {'Child2': [HTSeq.GenomicInterval('I', 39, 50, '-')]}
        exp_classes = {}
        exp_filtered = {"GrandParent", "ParentWithGrandparent", "Parent2", "Child1", "Sibling"}
        exp_parents = {'ParentWithGrandparent': 'GrandParent', 'Child1': 'ParentWithGrandparent', 'Child2': 'Parent2'}

        self.assertEqual(alias, exp_alias)
        self.assertEqual(idents, exp_ident)
        self.assertEqual(intervals, exp_intervals)
        self.assertEqual(rt.parents, exp_parents)
        self.assertEqual(rt.filtered, exp_filtered)
        self.assertEqual(rt.classes, exp_classes)
        self.assertEqual(list(feats.chrom_vectors['I']['-'].array.get_steps()), exp_feats)
        self.clear_filters()

    """Does ReferenceTables.get() properly handle type filters for discontinuous features?"""

    def test_ref_tables_type_filter(self):

        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        selection_rules = self.rules_template

        rt = ReferenceTables(feature_source, selection_rules, type_filter=["CDS"])
        feats, ident, alias, intervals, classes = rt.get()

        exp_alias = {'Child1': ['SharedName']}
        exp_attrs = {'Child1': [('Class', ('NA',)), ('Name', ('SharedName',))]}
        exp_feats = [(0, 29, set()), (29, 40, {'Child1'}), (40, sys.maxsize, set())]
        exp_intervals = {'Child1': [HTSeq.GenomicInterval('I', 29, 40, '-')]}
        exp_filtered = {"GrandParent", "ParentWithGrandparent", "Parent2", "Child2", "Sibling"}
        exp_parents = {'ParentWithGrandparent': 'GrandParent', 'Child1': 'ParentWithGrandparent', 'Child2': 'Parent2'}

        self.assertEqual(alias, exp_alias)
        self.assertEqual(ident, exp_attrs)
        self.assertEqual(intervals, exp_intervals)
        self.assertEqual(rt.parents, exp_parents)
        self.assertEqual(rt.filtered, exp_filtered)
        self.assertEqual(list(feats.chrom_vectors['I']['-'].array.get_steps()), exp_feats)
        self.clear_filters()

    """Does ReferenceTables.get() properly handle both source and type filters for discontinuous features?"""

    def test_ref_tables_both_filter(self):

        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        selection_rules = MockFeatureSelector(self.rules_template)

        rt = ReferenceTables(feature_source, selection_rules, source_filter=["SourceName"], type_filter=["gene"])
        feats, alias, classes = rt.get()

        self.assertEqual(rt.filtered, {'Child1', 'Child2'})
        self.assertEqual(rt.parents, {'ParentWithGrandparent': 'GrandParent', 'Child1': 'ParentWithGrandparent', 'Child2': 'Parent2'})
        # self.assertEqual(list(attrs.keys()), ['GrandParent', 'Parent2', 'Sibling'])
        # self.assertEqual(list(intervals.keys()), ['GrandParent', 'Parent2', 'Sibling'])
        self.assertEqual(list(alias.keys()), ['GrandParent', 'Parent2', 'Sibling'])
        self.assertEqual(len(list(feats.chrom_vectors['I']['-'].array.get_steps())), 8)
        self.clear_filters()

    def clear_filters(self):
        """Since the filters in ReferenceTables are class attributes, they must be cleared.
        Otherwise they will interfere with subsequent tests."""

        ReferenceTables.source_filter = []
        ReferenceTables.type_filter = []

if __name__ == '__main__':
    unittest.main()
