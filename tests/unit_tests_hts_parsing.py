import collections
import unittest
from copy import deepcopy
from unittest.mock import patch, mock_open, call

from rna.counter.features import FeatureSelector
from tiny.rna.counter.hts_parsing import *
# from tests.unit_tests_counter import resources
import unit_test_helpers as helpers
resources = "./testdata/counter"


class MockFeatureSelector:
    def __init__(self, rules_table):
        self.rules_table = FeatureSelector.build_selectors(rules_table)
        self.inv_ident = FeatureSelector.build_inverted_identities(self.rules_table)


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

    def get_gff_attr_string(self, gff_line):
        return gff_line.split('\t')[-1]

    def parse_gff_attr(self, gff_file_content):
        attr_str = self.get_gff_attr_string(gff_file_content)
        return parse_GFF_attribute_string(attr_str)

    def selector_with_rules(self, updates_list):
        """Returns a MockFeatureSelector with the specified updates to the default rule template"""

        rules = [deepcopy(helpers.rules_template[0]) for _ in range(len(updates_list))]
        for changes, template in zip(updates_list, rules):
            template.update(changes)
        return MockFeatureSelector(rules)

    def exhaust_iterator(self, it):
        collections.deque(it, maxlen=0)

    # === TESTS ===

    """Did SAM_reader correctly skip header values and parse all pertinent info from a single record SAM file?"""

    def test_sam_reader(self):
        sam_bundle = next(SAM_reader().bundle_multi_alignments(self.short_sam_file))
        sam_record = sam_bundle[0]

        self.assertEqual(sam_record['chrom'], "I")
        self.assertEqual(sam_record['start'], 15064569)
        self.assertEqual(sam_record['end'], 15064590)
        self.assertEqual(sam_record['strand'], '-')
        self.assertEqual(sam_record['name'], "0_count=5")
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
        ours = SAM_reader().bundle_multi_alignments(file)
        theirs = HTSeq.bundle_multiple_alignments(HTSeq.BAM_Reader(file))

        for our_bundle, their_bundle in zip(ours, theirs):
            self.assertEqual(len(our_bundle), len(their_bundle))
            for our, their in zip(our_bundle, their_bundle):
                self.assertEqual(our['chrom'], their.iv.chrom)
                self.assertEqual(our['start'], their.iv.start)
                self.assertEqual(our['end'], their.iv.end)
                self.assertEqual(our['name'], their.read.name)
                self.assertEqual(our['nt5'], chr(their.read.seq[0]))  # See note above
                self.assertEqual(our['strand'], their.iv.strand)
                if our['strand'] == '-':                              # See note above
                    self.assertEqual(our['seq'][::-1].translate(helpers.complement), their.read.seq)
                else:
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

    """Does ReferenceTables.get() return the expected features, aliases, and classes for a single record GFF?"""

    def test_ref_tables_single_feature(self):
        feature_source = {self.short_gff_file: ["sequence_name"]}
        feature_selector = self.selector_with_rules([
            {'Identity': ("Class", "CSR"), 'Strand': "+", 'Hierarchy': 1, 'nt5end': "N/A", 'Length': "20",
             'Strict': True},
            {'Identity': ("biotype", "snoRNA"), 'Strand': "-", 'Hierarchy': 2, 'nt5end': "N/A", 'Length': "30",
             'Strict': False}
        ])
        iv = HTSeq.GenomicInterval("I", 3746, 3908, "-")
        kwargs = {'all_features': True}

        feats, alias, classes = ReferenceTables(feature_source, feature_selector, **kwargs).get()
        steps = list(feats[iv].array[iv.start:iv.end].get_steps(values_only=True))

        self.assertEqual((type(feats), type(alias), type(classes)), (HTSeq.GenomicArrayOfSets, dict, dict))
        self.assertEqual(steps, [{("Gene:WBGene00023193", 3746, 3909, '-', ((1, 2, False),))}])
        self.assertEqual(alias, {'Gene:WBGene00023193': ('Y74C9A.6',)})
        self.assertEqual(classes, {'Gene:WBGene00023193': ('additional_class', 'unknown')})

    """Repeating the previous test with all_features=False should produce the same result for this test."""

    def test_ref_tables_single_feature_all_features_false(self):
        kwargs = {'all_features': False}
        feature_source = {self.short_gff_file: ["sequence_name"]}
        feature_selector = self.selector_with_rules([
            {'Identity': ("Class", "CSR"), 'Strand': "+", 'Hierarchy': 1, 'nt5end': "N/A", 'Length': "20",
             'Strict': True},
            {'Identity': ("biotype", "snoRNA"), 'Strand': "-", 'Hierarchy': 2, 'nt5end': "N/A", 'Length': "30",
             'Strict': False}
        ])
        iv = HTSeq.GenomicInterval("I", 3746, 3908, "-")

        feats, alias, classes = ReferenceTables(feature_source, feature_selector, **kwargs).get()
        steps = list(feats[iv].array[iv.start:iv.end].get_steps(values_only=True))

        self.assertEqual((type(feats), type(alias), type(classes)), (HTSeq.GenomicArrayOfSets, dict, dict))
        self.assertEqual(steps, [{("Gene:WBGene00023193", 3746, 3909, '-', ((1, 2, False),))}])
        self.assertEqual(alias, {'Gene:WBGene00023193': ('Y74C9A.6',)})
        self.assertEqual(classes, {'Gene:WBGene00023193': ('additional_class', 'unknown')})

    """Does ReferenceTables.get() raise ValueError when a Name Attribute refers to a missing attribute?"""

    def test_ref_tables_missing_name_attribute(self):
        bad = "bad_name_attribute"
        feature_source = {self.short_gff_file: [bad]}
        feature_selector = MockFeatureSelector([])
        kwargs = {'all_features': True}

        expected_err = f"Feature Gene:WBGene00023193 does not contain a '{bad}' attribute." + '\n'
        expected_err += f"Error occurred on line 1 of {self.short_gff_file}"

        with self.assertRaisesRegex(ValueError, expected_err):
            ReferenceTables(feature_source, feature_selector, **kwargs).get()

    """Repeating previous test with all_features=False as this yields different results"""

    def test_ref_tables_missing_name_attribute_all_features_false(self):
        kwargs = {'all_features': False}
        bad = "bad_name_attribute"
        feature_source = {self.short_gff_file: [bad]}
        feature_selector = MockFeatureSelector([])

        expected_err = "No features or classes were retained while parsing your GFF file.\n" \
                       "This may be due to a lack of features matching 'Select for...with value...'"

        # Since all_features is False and there are no identity matches, the main loop in
        # ReferenceTables.get() skips the steps for recording the feature's alias.
        # Instead, a different exception is raised due to reference tables being empty
        with self.assertRaisesRegex(ValueError, expected_err):
            ReferenceTables(feature_source, feature_selector, **kwargs).get()

    """Does ReferenceTables.get() raise ValueError when a feature lacks an ID attribute?"""

    def test_ref_tables_missing_id_attribute(self):
        feature_source = {self.short_gff_file: ["ID"]}
        feature_selector = self.selector_with_rules(helpers.rules_template)
        kwargs = {'all_features': True}

        gff_row_without_id = helpers.read(self.short_gff_file).replace('ID=Gene:WBGene00023193;', '')
        mock_reader = mock_open(read_data=gff_row_without_id)

        expected_err = f"Feature WBGene00023193 does not contain a 'ID' attribute.\n"
        expected_err += f"Error occurred on line 1 of {self.short_gff_file}"

        with patch('tiny.rna.counter.hts_parsing.HTSeq.utils.open', new=mock_reader):
            with self.assertRaisesRegex(ValueError, expected_err):
                _ = ReferenceTables(feature_source, feature_selector, **kwargs).get()

    """Does ReferenceTables.get() properly concatenate aliases if there is more than one alias for a feature?"""
    """Does ReferenceTables.get() properly concatenate aliases when Name Attribute refers to a list-type alias?"""
    # 2 for 1!

    def test_ref_tables_alias_multisource_concat(self):
        feature_source = {self.short_gff_file: ["ID", "Class"]}
        kwargs = {'all_features': True}

        # Notice: screening for "ID" name attribute happens earlier in counter.load_config()
        expected_alias = {"Gene:WBGene00023193": ("Gene:WBGene00023193", "additional_class", "unknown")}
        _, alias, _ = ReferenceTables(feature_source, MockFeatureSelector([]), **kwargs).get()

        self.assertDictEqual(alias, expected_alias)

    """Repeating previous test with all_features=False as this yields different results"""

    def test_ref_tables_alias_multisource_concat_all_features_false(self):
        feature_source = {self.short_gff_file: ["ID", "Class"]}
        kwargs = {'all_features': False}

        expected_err = "No features or classes were retained while parsing your GFF file.\n" \
                       "This may be due to a lack of features matching 'Select for...with value...'"

        with self.assertRaisesRegex(ValueError, expected_err):
            # No aliases saved due to all_features=False and the lack of identity matches
            _, alias, _ = ReferenceTables(feature_source, MockFeatureSelector([]), **kwargs).get()

    """Does ReferenceTables.get() properly concatenate identity match tuples when multiple GFF files define
    matches for a feature?"""

    def test_ref_tables_identity_matches_multisource_concat(self):
        feature_source = {
            self.short_gff_file: ["ID"],
            f"{resources}/single2.gff3": ["ID"]
        }

        kwargs = {'all_features': True}
        feature_selector = self.selector_with_rules([
            {'Identity': ('Name', 'WBGene00023193b'), 'Hierarchy': 1},
            {'Identity': ('Name', 'WBGene00023193'), 'Hierarchy': 2},
            {'Identity': ('biotype', 'snoRNA'), 'Hierarchy': 3}
        ])

        expected_matches = [
            set(),
            {('Gene:WBGene00023193', 3746, 3909, '-', ((0, 1, True), (1, 2, True), (2, 3, True)))},
            set()
        ]

        feats, _, _ = ReferenceTables(feature_source, feature_selector, **kwargs).get()

        actual_idents = list(feats.chrom_vectors['I']['.'].array.get_steps(values_only=True))
        for act_attr, exp_attr in zip(actual_idents, expected_matches):
            self.assertEqual(act_attr, exp_attr)

    """Does ReferenceTables.get() properly handle aliases for discontinuous features?"""

    def test_ref_tables_discontinuous_aliases(self):
        kwargs = {'all_features': True}
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        mock_selector = self.selector_with_rules(helpers.rules_template)

        _, alias, _ = ReferenceTables(feature_source, mock_selector, **kwargs).get()

        # Ancestor depth of 1, distinct aliases
        self.assertEqual(alias['Parent2'], ('Child2Name', 'Parent2Name'))
        # Ancestor depth >1, shared aliases
        self.assertEqual(alias['GrandParent'], ('SharedName',))
        # Siblings, distinct aliases
        self.assertEqual(alias['Sibling'], ('Sibling1', 'Sibling2', 'Sibling3'))

    """If all_features=False and there are no identity matches, are discontinuous features correctly omitted?"""

    def test_ref_tables_discontinuous_no_match_all_features_false(self):
        kwargs = {'all_features': False}
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        mock_selector = self.selector_with_rules(helpers.rules_template)

        expected_err = "No features or classes were retained while parsing your GFF file.\n" \
                       "This may be due to a lack of features matching 'Select for...with value...'"

        with self.assertRaisesRegex(ValueError, expected_err):
            ReferenceTables(feature_source, mock_selector, **kwargs).get()

    """Does ReferenceTables.get() properly handle intervals for discontinous features?"""

    def test_ref_tables_discontinuous_intervals(self):
        kwargs = {'all_features': True}
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        feature_selector = self.selector_with_rules(helpers.rules_template)

        grandparent_iv = HTSeq.GenomicInterval('I', 0, 10, '-')
        parent_w_p_iv = HTSeq.GenomicInterval('I', 9, 20, '-')
        child_w_gp_iv = HTSeq.GenomicInterval('I', 29, 40, '-')
        parent_2 = HTSeq.GenomicInterval('I', 19, 30, '-')
        child_2 = HTSeq.GenomicInterval('I', 39, 50, '-')
        sib_1 = HTSeq.GenomicInterval('I', 99, 110, '-')
        sib_2 = HTSeq.GenomicInterval('I', 110, 120, '-')
        sib_3 = HTSeq.GenomicInterval('I', 139, 150, '-')

        RT_instance = ReferenceTables(feature_source, feature_selector, **kwargs)
        _ = RT_instance.get()

        # Ancestor depth of 1
        self.assertEqual(RT_instance.intervals['GrandParent'], [grandparent_iv, parent_w_p_iv, child_w_gp_iv])
        # Ancestor depth >1
        self.assertEqual(RT_instance.intervals['Parent2'], [parent_2, child_2])
        # Siblings
        self.assertEqual(RT_instance.intervals['Sibling'], [sib_1, sib_2, sib_3])

    """Does ReferenceTables.get() properly merge identity matches of discontinuous features with the root feature?"""

    def test_ref_tables_discontinuous_identity_matches(self):

        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        feature_selector = self.selector_with_rules([
            {'Identity': ('Class', 'NA'), 'Hierarchy': 2},                  # Rule 1
            {'Identity': ('Name', 'Sibling3'), 'Hierarchy': 3},             # Rule 2
            {'Identity': ('UniqueAttr', 'FirstSibling'), 'Hierarchy': 0}    # Rule 3
        ])
        rt_kwargs = {'all_features': True}
        
        rule1_match = (0, 2, True)
        rule2_match = (1, 3, True)
        rule3_match = (2, 0, True)

        expected = [{('GrandParent', 0, 20, '-', (rule1_match,))},
                    {('GrandParent', 0, 20, '-', (rule1_match,)), ('Parent2', 19, 30, '-', (rule1_match,))},
                    {('Parent2', 19, 30, '-', (rule1_match,))},
                    {('Parent2', 19, 30, '-', (rule1_match,)), ('GrandParent', 29, 40, '-', (rule1_match,))},
                    {('GrandParent', 29, 40, '-', (rule1_match,))},
                    {('GrandParent', 29, 40, '-', (rule1_match,)), ('Parent2', 39, 50, '-', (rule1_match,))},
                    {('Parent2', 39, 50, '-', (rule1_match,))},
                    set(),
                    {('Sibling', 99, 110, '-', (rule3_match, rule2_match))},  # Note: sorted by rank, not rule index
                    {('Sibling', 110, 120, '-', (rule3_match, rule2_match))},
                    set(),
                    {('Sibling', 139, 150, '-', (rule3_match, rule2_match))},
                    set()]

        feats, _, _ = ReferenceTables(feature_source, feature_selector, **rt_kwargs).get()

        for act, exp in zip(feats.chrom_vectors["I"]["."].array.get_steps(values_only=True), expected):
            # For this test we are only interested in the match tuples for each feature
            self.assertEqual(act, exp)

    """Does ReferenceTables.get() properly build a GenomicArrayOfSets for discontinuous features?"""

    def test_ref_tables_discontinuous_features(self):

        kwargs = {'all_features': True}
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        feature_selector = self.selector_with_rules(helpers.rules_template)

        expected = [{('GrandParent', 0, 20, '-', ())},
                    {('Parent2', 19, 30, '-', ()), ('GrandParent', 0, 20, '-', ())},
                    {('Parent2', 19, 30, '-', ())},
                    {('Parent2', 19, 30, '-', ()), ('GrandParent', 29, 40, '-', ())},
                    {('GrandParent', 29, 40, '-', ())},
                    {('Parent2', 39, 50, '-', ()), ('GrandParent', 29, 40, '-', ())},
                    {('Parent2', 39, 50, '-', ())},
                    set(),
                    {('Sibling', 99, 110, '-', ())},
                    {('Sibling', 110, 120, '-', ())},
                    set(),
                    {('Sibling', 139, 150, '-', ())},
                    set()]

        feats, _, _, = ReferenceTables(feature_source, feature_selector, **kwargs).get()

        for act, exp in zip(feats.chrom_vectors["I"]["."].array.get_steps(values_only=True), expected):
            self.assertEqual(act, exp)

    """Does ReferenceTables.get() properly handle source filters for discontinuous features?"""

    def test_ref_tables_source_filter(self):

        kwargs = {'source_filter': ["Source2Name"], 'all_features': False}
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        selection_rule = dict(helpers.rules_template[0], Identity=('Parent', 'Parent2'))
        feature_selector = self.selector_with_rules([selection_rule])

        exp_alias = {'Child2': ('Child2Name',)}
        exp_feats = [set(), {('Child2', 39, 50, '-', ((0, 0, True),))}, set()]
        exp_intervals = {'Child2': [HTSeq.GenomicInterval('I', 39, 50, '-')]}
        exp_classes = {'Child2': ('NA',)}
        exp_filtered = {"GrandParent", "ParentWithGrandparent", "Parent2", "Child1", "Sibling"}
        exp_parents = {'ParentWithGrandparent': 'GrandParent', 'Child1': 'ParentWithGrandparent', 'Child2': 'Parent2'}

        rt = ReferenceTables(feature_source, feature_selector, **kwargs)
        feats, alias, classes = rt.get()

        self.assertEqual(alias, exp_alias)
        self.assertEqual(rt.intervals, exp_intervals)
        self.assertEqual(rt.parents, exp_parents)
        self.assertEqual(rt.filtered, exp_filtered)
        self.assertEqual(classes, exp_classes)
        self.assertEqual(list(feats.chrom_vectors['I']['.'].array.get_steps(values_only=True)), exp_feats)
        self.clear_filters()

    """Does ReferenceTables.get() properly handle type filters for discontinuous features?"""

    def test_ref_tables_type_filter(self):

        kwargs = {'type_filter': ["CDS"], 'all_features': False}
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        selection_rule = dict(helpers.rules_template[0], Identity=('Parent', 'ParentWithGrandparent'))
        feature_selector = self.selector_with_rules([selection_rule])

        exp_alias = {'Child1': ('SharedName',)}
        exp_feats = [set(), {('Child1', 29, 40, '-', ((0, 0, True),))}, set()]
        exp_intervals = {'Child1': [HTSeq.GenomicInterval('I', 29, 40, '-')]}
        exp_filtered = {"GrandParent", "ParentWithGrandparent", "Parent2", "Child2", "Sibling"}
        exp_parents = {'ParentWithGrandparent': 'GrandParent', 'Child1': 'ParentWithGrandparent', 'Child2': 'Parent2'}

        rt = ReferenceTables(feature_source, feature_selector, **kwargs)
        feats, alias, classes = rt.get()

        self.assertEqual(alias, exp_alias)
        self.assertEqual(rt.intervals, exp_intervals)
        self.assertEqual(rt.parents, exp_parents)
        self.assertEqual(rt.filtered, exp_filtered)
        self.assertEqual(list(feats.chrom_vectors['I']['.'].array.get_steps(values_only=True)), exp_feats)
        self.clear_filters()

    """Does ReferenceTables.get() properly handle both source and type filters for discontinuous features?"""

    def test_ref_tables_both_filter(self):

        kwargs = {'source_filter': ["SourceName"], 'type_filter': ["gene"], 'all_features': True}
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        feature_selector = self.selector_with_rules(helpers.rules_template)

        rt = ReferenceTables(feature_source, feature_selector, **kwargs)
        feats, alias, classes = rt.get()

        self.assertEqual(rt.filtered, {'Child1', 'Child2'})
        self.assertEqual(rt.parents, {'ParentWithGrandparent': 'GrandParent', 'Child1': 'ParentWithGrandparent', 'Child2': 'Parent2'})
        self.assertEqual(list(classes.keys()), ['GrandParent', 'Parent2', 'Sibling'])
        self.assertEqual(list(alias.keys()), ['GrandParent', 'Parent2', 'Sibling'])
        self.assertEqual(len(list(feats.chrom_vectors['I']['.'].array.get_steps(values_only=True))), 9)
        self.clear_filters()

    def clear_filters(self):
        """Since the filters in ReferenceTables are class attributes, they must be cleared.
        Otherwise they will interfere with subsequent tests."""

        ReferenceTables.source_filter = []
        ReferenceTables.type_filter = []

    """Does SAM_reader._get_decollapsed_filename() create an appropriate filename?"""

    def test_SAM_reader_get_decollapsed_filename(self):
        reader = SAM_reader()
        reader.file = "~/path/to/input/sam_file.sam"

        sam_out = reader._get_decollapsed_filename()

        self.assertEqual(sam_out, "sam_file_decollapsed.sam")

    """Does SAM_reader._read_thru_header() correctly identify header lines and write them to the decollapsed file?"""

    def test_SAM_reader_read_thru_header(self):
        reader = SAM_reader(decollapse=True)
        reader._decollapsed_filename = "mock_outfile_name.sam"

        with open(self.short_sam_file, 'rb') as sam_in:
            with patch('builtins.open', mock_open()) as sam_out:
                line = reader._read_thru_header(sam_in)

        expected_writelines = [
            call('mock_outfile_name.sam', 'w'),
            call().__enter__(),
            call().writelines(["@SQ	SN:I	LN:21\n"]),
            call().__exit__(None, None, None)
        ]

        sam_out.assert_has_calls(expected_writelines)
        self.assertTrue(len(reader._headers) == 1)

    """Does SAM_reader._write_decollapsed_sam() write the correct number of duplicates to the decollapsed file?"""

    def test_SAM_reader_write_decollapsed_sam(self):
        reader = SAM_reader(decollapse=True)
        reader._decollapsed_reads = [(b"0_count=5", b"mock line from SAM file")]
        reader._decollapsed_filename = "mock_outfile_name.sam"

        expected_writelines = [
            call('mock_outfile_name.sam', 'ab'),
            call().__enter__(),
            call().writelines([b"mock line from SAM file"] * 5),
            call().__exit__(None, None, None)
        ]

        with patch('builtins.open', mock_open()) as outfile:
            reader._write_decollapsed_sam()

        outfile.assert_has_calls(expected_writelines)
        self.assertTrue(len(reader._decollapsed_reads) == 0)

    """Does SAM_reader._parse_alignments() save lines and write them to the decollapsed file when appropriate?"""

    def test_SAM_reader_parse_alignments_decollapse(self):
        with patch.object(SAM_reader, "_write_decollapsed_sam") as write_fn, \
                patch('tiny.rna.counter.hts_parsing.open', new_callable=mock_open) as mopen:

            reader = SAM_reader(decollapse=True)
            reader._decollapsed_reads = [0] * 99999     # At 100,001, buffer will be written
            reader.file = self.short_sam_file           # File with single alignment

            with open(self.short_sam_file, 'rb') as sam_in:
                self.exhaust_iterator(reader._parse_alignments(sam_in))
                write_fn.assert_not_called()

                # Rewind and add one more alignment to push it over threshold
                sam_in.seek(0)
                self.exhaust_iterator(reader._parse_alignments(sam_in))
                write_fn.assert_called_once()




if __name__ == '__main__':
    unittest.main()
