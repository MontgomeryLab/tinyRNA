import collections
import contextlib
import unittest
import HTSeq
import io

from copy import deepcopy
from random import randint
from unittest.mock import patch, mock_open, call

from tiny.rna.counter.features import FeatureSelector
from tiny.rna.counter.matching import IntervalPartialMatch
from tiny.rna.counter.statistics import LibraryStats
from tiny.rna.counter.hts_parsing import *

import unit_test_helpers as helpers

resources = "./testdata/counter"


class MyTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.gff_file = f"{resources}/identity_choice_test.gff3"
        self.short_gff_file = f"{resources}/single.gff3"
        self.short_gff = helpers.read(self.short_gff_file)

        self.sam_file = f"{resources}/identity_choice_test.sam"
        self.short_sam_file = f"{resources}/single.sam"
        self.short_sam = helpers.read(self.short_sam_file)

        self.maxDiff = None

    # === HELPERS ===

    def get_gff_attr_string(self, gff_line):
        return gff_line.split('\t')[-1]

    def parse_gff_attr(self, gff_file_content):
        attr_str = self.get_gff_attr_string(gff_file_content)
        return parse_GFF_attribute_string(attr_str)

    def selector_with_template(self, updates_list):
        """Returns a MockFeatureSelector with the specified updates to the default rule template"""

        rules = [deepcopy(helpers.rules_template[0]) for _ in range(len(updates_list))]
        for changes, template in zip(updates_list, rules):
            template.update(changes)
        return FeatureSelector(rules)

    def exhaust_iterator(self, it):
        collections.deque(it, maxlen=0)

    # === TESTS ===

    """Did SAM_reader correctly skip header values and parse all pertinent info from a single record SAM file?"""

    def test_sam_reader(self):
        sam_bundle, read_count = next(SAM_reader().bundle_multi_alignments(self.short_sam_file))
        sam_record = sam_bundle[0]

        self.assertEqual(sam_record['Chrom'], "I")
        self.assertEqual(sam_record['Start'], 15064569)
        self.assertEqual(sam_record['End'], 15064590)
        self.assertEqual(sam_record['Strand'], False)
        self.assertEqual(sam_record['Name'], b"0_count=5")
        self.assertEqual(sam_record['Seq'], b"CAAGACAGAGCTTCACCGTTC")
        self.assertEqual(sam_record['Length'], 21)
        self.assertEqual(sam_record['nt5end'], 'G')

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

        for (our_bundle, _), their_bundle in zip(ours, theirs):
            self.assertEqual(len(our_bundle), len(their_bundle))
            for our, their in zip(our_bundle, their_bundle):
                self.assertEqual(our['Chrom'], their.iv.chrom)
                self.assertEqual(our['Start'], their.iv.start)
                self.assertEqual(our['End'], their.iv.end)
                self.assertEqual(our['Name'].decode(), their.read.name)
                self.assertEqual(our['nt5end'], chr(their.read.seq[0]))  # See note above
                self.assertEqual(our['Strand'], helpers.strand_to_bool(their.iv.strand))
                if our['Strand'] is False:                               # See note above
                    self.assertEqual(our['Seq'][::-1].translate(helpers.complement), their.read.seq)
                else:
                    self.assertEqual(our['Seq'], their.read.seq)

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

    """Does ReferenceFeatures.get() return the expected features, aliases, and classes for a single record GFF?"""

    def test_ref_tables_single_feature(self):
        feature_source = {self.short_gff_file: ["sequence_name"]}
        feature_selector = self.selector_with_template([
            # Fails to match due to Identity selector
            {'Identity': ("Class", "CSR"), 'Strand': "sense", 'Hierarchy': 1, 'Class': 'none', 'nt5end': "all",
             'Overlap': 'nested', 'Length': "20"},
            # Match
            {'Identity': ("biotype", "snoRNA"), 'Strand': "antisense", 'Hierarchy': 2, 'Class': 'tag', 'nt5end': "all",
             'Overlap': 'partial', 'Length': "30"}
        ])
        iv = HTSeq.GenomicInterval("I", 3746, 3909, "-")
        kwargs = {'all_features': True}

        feats, alias, tags = ReferenceFeatures(feature_source, **kwargs).get(feature_selector)
        steps = list(feats[iv].array[iv.start:iv.end].get_steps(values_only=True))

        self.assertEqual((type(feats), type(alias), type(tags)), (HTSeq.GenomicArrayOfSets, dict, defaultdict))
        self.assertEqual(steps, [{(("Gene:WBGene00023193", 'tag'), False, ((1, 2, IntervalPartialMatch(iv)),))}])
        self.assertEqual(alias, {"Gene:WBGene00023193": ('Y74C9A.6',)})
        self.assertDictEqual(tags, {"Gene:WBGene00023193": {('Gene:WBGene00023193', 'tag')}})

    """Repeating the previous test with all_features=False should produce the same result for this test."""

    def test_ref_tables_single_feature_all_features_false(self):
        kwargs = {'all_features': False}
        feature_source = {self.short_gff_file: ["sequence_name"]}
        feature_selector = self.selector_with_template([
            # Fails to match due to Identity selector
            {'Identity': ("Class", "CSR"), 'Strand': "sense", 'Hierarchy': 1, 'Class': 'none', 'nt5end': "all",
             'Overlap': 'nested', 'Length': "20"},
            # Match
            {'Identity': ("biotype", "snoRNA"), 'Strand': "antisense", 'Hierarchy': 2, 'Class': 'tag', 'nt5end': "all",
             'Overlap': 'partial', 'Length': "30"}
        ])
        iv = HTSeq.GenomicInterval("I", 3746, 3909, "-")
        kwargs = {'all_features': True}

        feats, alias, tags = ReferenceFeatures(feature_source, **kwargs).get(feature_selector)
        steps = list(feats[iv].array[iv.start:iv.end].get_steps(values_only=True))

        self.assertEqual((type(feats), type(alias), type(tags)), (HTSeq.GenomicArrayOfSets, dict, defaultdict))
        self.assertEqual(steps, [{(("Gene:WBGene00023193", 'tag'), False, ((1, 2, IntervalPartialMatch(iv)),))}])
        self.assertDictEqual(alias, {"Gene:WBGene00023193": ('Y74C9A.6',)})
        self.assertDictEqual(tags, {"Gene:WBGene00023193": {('Gene:WBGene00023193', 'tag')}})

    """Repeating previous test with all_features=False as this yields different results"""

    def test_ref_tables_missing_name_attribute_all_features_false(self):
        kwargs = {'all_features': False}
        bad = "bad_name_attribute"
        feature_source = {self.short_gff_file: [bad]}
        feature_selector = FeatureSelector([])

        expected_err = "No features were retained while parsing your GFF file.\n" \
                       "This may be due to a lack of features matching 'Select for...with value...'"

        # Since all_features is False and there are no identity matches, the main loop in
        # ReferenceFeatures.get() skips the steps for recording the feature's alias.
        # Instead, a different exception is raised due to reference tables being empty
        with self.assertRaisesRegex(ValueError, expected_err):
            ReferenceFeatures(feature_source, **kwargs).get(feature_selector)

    """Does ReferenceFeatures.get() raise ValueError when a feature lacks an ID attribute?"""

    def test_ref_tables_missing_id_attribute(self):
        feature_source = {self.short_gff_file: ["ID"]}
        feature_selector = self.selector_with_template(helpers.rules_template)
        kwargs = {'all_features': True}

        gff_row_without_id = helpers.read(self.short_gff_file).replace('ID=Gene:WBGene00023193;', '')
        mock_reader = mock_open(read_data=gff_row_without_id)

        expected_err = f"Feature WBGene00023193 does not have an ID attribute.\n"
        expected_err += f"Error occurred on line 1 of {self.short_gff_file}"

        with patch('tiny.rna.counter.hts_parsing.HTSeq.utils.open', new=mock_reader):
            with self.assertRaisesRegex(ValueError, expected_err):
                _ = ReferenceFeatures(feature_source, **kwargs).get(feature_selector)

    """Does ReferenceFeatures.get() properly concatenate aliases if there is more than one alias for a feature?"""
    """Does ReferenceFeatures.get() properly concatenate aliases when Name Attribute refers to a list-type alias?"""
    # 2 for 1!

    def test_ref_tables_alias_multisource_concat(self):
        feature_source = {self.short_gff_file: ["ID", "Class"]}
        kwargs = {'all_features': True}

        # Notice: screening for "ID" name attribute happens earlier in counter.load_config()
        expected_alias = {"Gene:WBGene00023193": ("additional_class", "Gene:WBGene00023193", "unknown")}
        _, alias, _ = ReferenceFeatures(feature_source, **kwargs).get(FeatureSelector([]))

        self.assertDictEqual(alias, expected_alias)

    """Repeating previous test with all_features=False as this yields different results"""

    def test_ref_tables_alias_multisource_concat_all_features_false(self):
        feature_source = {self.short_gff_file: ["ID", "Class"]}
        kwargs = {'all_features': False}

        expected_err = "No features were retained while parsing your GFF file.\n" \
                       "This may be due to a lack of features matching 'Select for...with value...'"

        with self.assertRaisesRegex(ValueError, expected_err):
            # No aliases saved due to all_features=False and the lack of identity matches
            _, alias, _ = ReferenceFeatures(feature_source, **kwargs).get(FeatureSelector([]))

    """Does ReferenceFeatures.get() properly concatenate identity match tuples when multiple GFF files define
    matches for a feature?"""

    def test_ref_tables_identity_matches_multisource_concat(self):
        feature_source = {
            self.short_gff_file: ["ID"],
            f"{resources}/single2.gff3": ["ID"]
        }

        kwargs = {'all_features': True}
        feature_selector = self.selector_with_template([
            {'Identity': ('Name', 'WBGene00023193b'), 'Hierarchy': 1},
            {'Identity': ('Name', 'WBGene00023193'), 'Hierarchy': 2},
            {'Identity': ('biotype', 'snoRNA'), 'Hierarchy': 3}
        ])

        ivm = IntervalPartialMatch(HTSeq.GenomicInterval('n/a', 3746, 3909))

        expected_matches = [
            set(),
            {(('Gene:WBGene00023193', ''), False, ((0, 1, ivm), (1, 2, ivm), (2, 3, ivm)))},
            set()
        ]

        feats, _, _ = ReferenceFeatures(feature_source, **kwargs).get(feature_selector)

        actual_matches = list(feats.chrom_vectors['I']['.'].array.get_steps(values_only=True))
        self.assertListEqual(actual_matches, expected_matches)

    """Does ReferenceFeatures.get() properly handle aliases for discontinuous features?"""

    def test_ref_tables_discontinuous_aliases(self):
        kwargs = {'all_features': True}
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        mock_selector = self.selector_with_template(helpers.rules_template)

        _, alias, _ = ReferenceFeatures(feature_source, **kwargs).get(mock_selector)

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
        mock_selector = self.selector_with_template([{'Identity': ('No', 'Match')}])

        expected_err = "No features were retained while parsing your GFF file.\n" \
                       "This may be due to a lack of features matching 'Select for...with value...'"

        with self.assertRaisesRegex(ValueError, expected_err):
            ReferenceFeatures(feature_source, **kwargs).get(mock_selector)

    """Does ReferenceFeatures.get() properly build a GenomicArrayOfSets for discontinuous features
        where no features match but all_features is True? If the feature didn't match we only want to
        display it in the Feature ID column of counts table with 0 counts. We DON'T want it to pop
        up as a candidate in Stage 2 selection."""

    def test_ref_tables_discontinuous_features(self):

        kwargs = {'all_features': True}
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        feature_selector = self.selector_with_template([{'Identity': ('No', 'Match')}])

        # Features that fail to match on identity are not added to the StepVector,
        # EVEN if all_features = True.
        expected = [set()]

        feats, _, _ = ReferenceFeatures(feature_source, **kwargs).get(feature_selector)
        actual = list(feats.chrom_vectors["I"]["."].array.get_steps(values_only=True))
        self.assertListEqual(actual, expected)

    """Does ReferenceFeatures.get() properly handle intervals for discontinous features?"""

    def test_ref_tables_discontinuous_intervals(self):
        kwargs = {'all_features': True}
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        feature_selector = self.selector_with_template(helpers.rules_template)

        grandparent_iv = HTSeq.GenomicInterval('I', 0, 10, '-')
        parent_w_p_iv = HTSeq.GenomicInterval('I', 9, 20, '-')
        child_w_gp_iv = HTSeq.GenomicInterval('I', 29, 40, '-')
        parent_2 = HTSeq.GenomicInterval('I', 19, 30, '-')
        child_2 = HTSeq.GenomicInterval('I', 39, 50, '-')
        sib_1 = HTSeq.GenomicInterval('I', 99, 110, '-')
        sib_2 = HTSeq.GenomicInterval('I', 110, 120, '-')
        sib_3 = HTSeq.GenomicInterval('I', 139, 150, '-')

        RT_instance = ReferenceFeatures(feature_source, **kwargs)
        _ = RT_instance.get(feature_selector)

        # Ancestor depth of 1
        self.assertEqual(RT_instance.intervals['GrandParent'], [grandparent_iv, parent_w_p_iv, child_w_gp_iv])
        # Ancestor depth >1
        self.assertEqual(RT_instance.intervals['Parent2'], [parent_2, child_2])
        # Siblings
        self.assertEqual(RT_instance.intervals['Sibling'], [sib_1, sib_2, sib_3])

    """Does ReferenceFeatures.get() properly merge identity matches of discontinuous features with the root feature?
    Identity match tuples now also contain the corresponding rule's IntervalSelector, so extra bookkeeping must be
    performed for intervals in this test."""

    def test_ref_tables_discontinuous_identity_matches(self):
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        feature_selector = self.selector_with_template([
            {'Identity': ('Class', 'NA'), 'Hierarchy': 2},                  # Rule 1
            {'Identity': ('Name', 'Sibling3'), 'Hierarchy': 3},             # Rule 2
            {'Identity': ('UniqueAttr', 'FirstSibling'), 'Hierarchy': 0}    # Rule 3
        ])
        rt_kwargs = {'all_features': True}

        gp_ivs =  [HTSeq.GenomicInterval('I', 0, 20, '-'), HTSeq.GenomicInterval('I', 29, 40, '-')]  # [0,10)+[9,20)->[0,20)
        p2_ivs =  [HTSeq.GenomicInterval('I', 19, 30, '-'), HTSeq.GenomicInterval('I', 39, 50, '-')]
        sib_ivs = [HTSeq.GenomicInterval('I', 99, 110, '-'), HTSeq.GenomicInterval('I', 110, 120, '-'),
                   HTSeq.GenomicInterval('I', 139, 150, '-')]

        rule1_gp =  {f"{iv.start}:{iv.end}": ((0, 2, IntervalPartialMatch(iv)),) for iv in gp_ivs}
        rule1_p2 =  {f"{iv.start}:{iv.end}": ((0, 2, IntervalPartialMatch(iv)),) for iv in p2_ivs}
        rule2_sib = {f"{iv.start}:{iv.end}": (1, 3, IntervalPartialMatch(iv))    for iv in sib_ivs}
        rule3_sib = {f"{iv.start}:{iv.end}": (2, 0, IntervalPartialMatch(iv))    for iv in sib_ivs}

        # For tables that store features in tagged form
        GrandParent, Parent2, Sibling = ('GrandParent',''), ('Parent2',''), ('Sibling','')

        expected = [{(GrandParent, False, rule1_gp['0:20'])},
                    {(GrandParent, False, rule1_gp['0:20']),  (Parent2,     False, rule1_p2["19:30"])},
                    {(Parent2,     False, rule1_p2["19:30"])},
                    {(Parent2,     False, rule1_p2["19:30"]), (GrandParent, False, rule1_gp['29:40'])},
                    {(GrandParent, False, rule1_gp['29:40'])},
                    {(GrandParent, False, rule1_gp['29:40']), (Parent2,     False, rule1_p2['39:50'])},
                    {(Parent2,     False, rule1_p2['39:50'])},
                    set(),
                    {(Sibling,     False, (rule3_sib['99:110'],  rule2_sib['99:110']))},  # Note: sorted by rank, not rule index
                    {(Sibling,     False, (rule3_sib['110:120'], rule2_sib['110:120']))},
                    set(),
                    {(Sibling,     False, (rule3_sib['139:150'], rule2_sib['139:150']))},
                    set()]

        feats, _, _ = ReferenceFeatures(feature_source, **rt_kwargs).get(feature_selector)
        actual_steps = list(feats.chrom_vectors["I"]["."].array.get_steps(values_only=True))
        self.assertListEqual(actual_steps, expected)

    """Does ReferenceFeatures.get() properly handle source filters for discontinuous features?"""

    def test_ref_tables_source_filter(self):

        kwargs = {'all_features': False}
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        feature_selector = self.selector_with_template([{'Filter_s': "Source2Name"}])

        child2_iv =     HTSeq.GenomicInterval('I', 39, 50, '-')
        exp_alias =     {'Child2': ('Child2Name',)}
        exp_feats =     [set(), {(('Child2', ''), False, ((0, 0, IntervalPartialMatch(child2_iv)),))}, set()]
        exp_intervals = {'Child2': [child2_iv]}
        exp_classes =   {'Child2': ('NA',)}
        exp_filtered =  {"GrandParent", "ParentWithGrandparent", "Parent2", "Child1", "Sibling"}
        exp_parents =   {'ParentWithGrandparent': 'GrandParent', 'Child1': 'ParentWithGrandparent', 'Child2': 'Parent2'}

        rt = ReferenceFeatures(feature_source, **kwargs)
        feats, alias, _ = rt.get(feature_selector)

        self.assertEqual(alias, exp_alias)
        self.assertEqual(rt.parents, exp_parents)
        self.assertEqual(rt.filtered, exp_filtered)
        self.assertEqual(list(feats.chrom_vectors['I']['.'].array.get_steps(values_only=True)), exp_feats)
        self.assertDictEqual(rt.intervals, exp_intervals)

    """Does ReferenceFeatures.get() properly handle type filters for discontinuous features?"""

    def test_ref_tables_type_filter(self):

        kwargs = {'all_features': False}
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        feature_selector = self.selector_with_template([{'Filter_t': "CDS"}])

        child1_iv =     HTSeq.GenomicInterval('I', 29, 40, '-')
        exp_alias =     {'Child1': ('SharedName',)}
        exp_feats =     [set(), {(('Child1', ''), False, ((0, 0, IntervalPartialMatch(child1_iv)),))}, set()]
        exp_intervals = {'Child1': [child1_iv]}
        exp_classes =   {'Child1': ('NA',)}
        exp_filtered =  {"GrandParent", "ParentWithGrandparent", "Parent2", "Child2", "Sibling"}
        exp_parents =   {'ParentWithGrandparent': 'GrandParent', 'Child1': 'ParentWithGrandparent', 'Child2': 'Parent2'}

        rt = ReferenceFeatures(feature_source, **kwargs)
        feats, alias, _ = rt.get(feature_selector)

        self.assertEqual(alias, exp_alias)
        self.assertEqual(rt.intervals, exp_intervals)
        self.assertEqual(rt.parents, exp_parents)
        self.assertEqual(rt.filtered, exp_filtered)
        self.assertEqual(list(feats.chrom_vectors['I']['.'].array.get_steps(values_only=True)), exp_feats)

    """Does ReferenceFeatures.get() properly handle both source and type filters for discontinuous features?"""

    def test_ref_tables_both_filter(self):

        kwargs = {'all_features': False}
        feature_source = {f"{resources}/discontinuous.gff3": ["Name"]}
        feature_selector = self.selector_with_template([{'Filter_s': "SourceName", 'Filter_t': "gene"}])

        rt = ReferenceFeatures(feature_source, **kwargs)
        feats, alias, _ = rt.get(feature_selector)

        self.assertEqual(rt.filtered, {'Child1', 'Child2'})
        self.assertEqual(rt.parents, {'ParentWithGrandparent': 'GrandParent', 'Child1': 'ParentWithGrandparent', 'Child2': 'Parent2'})
        self.assertEqual(list(alias.keys()), ['GrandParent', 'Parent2', 'Sibling'])
        self.assertEqual(len(list(feats.chrom_vectors['I']['.'].array.get_steps(values_only=True))), 9)

    """Does ReferenceFeatures.get() maintain correct records for a single feature matching tagged rules?"""

    def test_ref_tables_tagged_match_single(self):
        kwargs = {'all_features': False}
        feat_id = "Gene:WBGene00023193"
        feature_source = {f"{resources}/single.gff3": ["sequence_name"]}
        feature_selector = self.selector_with_template([
            {'Identity': ("ID", feat_id), 'Class': "tagged_match", 'Hierarchy': 1},
            {'Identity': ("ID", feat_id), 'Class': "",             'Hierarchy': 2}
        ])

        expected_tags = {feat_id: {(feat_id, "tagged_match"), (feat_id, '')}}
        expected_aliases = {feat_id: ('Y74C9A.6',)}
        iv = IntervalPartialMatch(HTSeq.GenomicInterval('n/a', 3746, 3909))
        expected_feats = [
            set(), {
                ((feat_id, 'tagged_match'), False, ((0, 1, iv),)),
                ((feat_id, ''),             False, ((1, 2, iv),))
            },
            set()
        ]

        feats, aliases, tags = ReferenceFeatures(feature_source, **kwargs).get(feature_selector)

        actual_feats = list(feats.chrom_vectors['I']['.'].array.get_steps(values_only=True))
        self.assertListEqual(actual_feats, expected_feats)
        self.assertDictEqual(aliases, expected_aliases)
        self.assertDictEqual(tags, expected_tags)

    """Does ReferenceFeatures.get() correctly merge records for discontinuous features matching multiple tagged rules?"""

    def test_ref_tables_tagged_match_merging(self):
        feature_source = {f"{resources}/discontinuous.gff3": ['Name']}

        # All rules match the same root feature
        feature_selector = self.selector_with_template([
            {'Identity': ("UniqueAttr", "child2"), 'Class': "shared", 'Hierarchy': 1},
            {'Identity': ("ID", "Parent2"),        'Class': "shared", 'Hierarchy': 2},
            {'Identity': ("ID", "Child2"),         'Class': "",       'Hierarchy': 3}
        ])

        expected_tags = {'Parent2': {('Parent2', 'shared'), ('Parent2', '')}}
        expected_aliases = {'Parent2': ('Child2Name', 'Parent2Name')}

        Parent2_iv = IntervalPartialMatch(HTSeq.GenomicInterval('n/a', 19, 30))
        Child2_iv = IntervalPartialMatch(HTSeq.GenomicInterval('n/a', 39, 50))
        expected_feats = [
            set(), {
                (('Parent2', 'shared'), False, ((0, 1, Parent2_iv), (1, 2, Parent2_iv))),
                (('Parent2', ''),       False, ((2, 3, Parent2_iv),)),
            },
            set(), {
                (('Parent2', 'shared'), False, ((0, 1, Child2_iv), (1, 2, Child2_iv))),
                (('Parent2', ''),       False, ((2, 3, Child2_iv),))
            },
            set()
        ]

        feats, aliases, tags = ReferenceFeatures(feature_source).get(feature_selector)

        stepvec = list(feats.chrom_vectors['I']['.'].array.get_steps(values_only=True))
        self.assertListEqual(stepvec, expected_feats)
        self.assertDictEqual(aliases, expected_aliases)
        self.assertDictEqual(tags, expected_tags)

    """Does SAM_reader._get_decollapsed_filename() create an appropriate filename?"""

    def test_SAM_reader_get_decollapsed_filename(self):
        reader = SAM_reader()
        reader.file = "~/path/to/input/sam_file.sam"

        sam_out = reader._get_decollapsed_filename()

        self.assertEqual(sam_out, "sam_file_decollapsed.sam")

    """Does SAM_reader._read_to_first_aln() correctly identify header lines and write them to the decollapsed file?"""

    def test_SAM_reader_read_thru_header(self):
        reader = SAM_reader(decollapse=True)
        reader._decollapsed_filename = "mock_outfile_name.sam"

        with open(self.short_sam_file, 'rb') as sam_in:
            with patch('builtins.open', mock_open()) as sam_out:
                line = reader._read_to_first_aln(sam_in)

        expected_writelines = [
            call('mock_outfile_name.sam', 'w'),
            call().__enter__(),
            call().writelines(["@SQ	SN:I	LN:21\n"]),
            call().__exit__(None, None, None)
        ]

        sam_out.assert_has_calls(expected_writelines)
        self.assertTrue(len(reader._header_lines) == 1)

    """Does SAM_reader._write_decollapsed_sam() write the correct number of duplicates to the decollapsed file?"""

    def test_SAM_reader_write_decollapsed_sam(self):
        reader = SAM_reader(decollapse=True)
        reader.collapser_type = "tiny-collapse"
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

    """Does SAM_reader report a single read count for non-collapsed SAM records?"""

    def test_SAM_reader_single_readcount_non_collapsed_SAM(self):
        # Read non-collapsed.sam but duplicate its single record twice
        with open(f"{resources}/non-collapsed.sam", 'rb') as f:
            sam_lines = f.readlines()
            sam_lines.extend([sam_lines[1]] * 2)
            mock_file = mock_open(read_data=b''.join(sam_lines))

        with patch('tiny.rna.counter.hts_parsing.open', new=mock_file):
            reader = SAM_reader()
            bundle, read_count = next(reader.bundle_multi_alignments('mock_file'))

        self.assertEqual(bundle[0]['Name'], b'NON_COLLAPSED_QNAME')
        self.assertEqual(len(bundle), 3)
        self.assertEqual(read_count, 1)

    """Are decollapsed outputs skipped when non-collapsed SAM files are supplied?"""

    def test_SAM_reader_no_decollapse_non_collapsed_SAM_files(self):
        stdout_capture = io.StringIO()
        with patch.object(SAM_reader, "_write_decollapsed_sam") as write_sam, \
                patch.object(SAM_reader, "_write_header_for_decollapsed_sam") as write_header:

            with contextlib.redirect_stderr(stdout_capture):
                reader = SAM_reader(decollapse=True)
                records = reader.bundle_multi_alignments(f"{resources}/non-collapsed.sam")
                self.exhaust_iterator(records)

        write_sam.assert_not_called()
        write_header.assert_not_called()
        self.assertEqual(reader.collapser_type, None)
        self.assertEqual(stdout_capture.getvalue(),
                         "Alignments do not appear to be derived from a supported collapser input. "
                         "Decollapsed SAM files will therefore not be produced.\n")

    """Does CaseInsensitiveAttrs correctly store, check membership, and retrieve?"""

    def test_CaseInsensitiveAttrs(self):
        def rand_case(string):
            return ''.join([l.upper() if randint(0,1) else l for l in string.lower()])

        cia = CaseInsensitiveAttrs()
        cia["AtTrKeY"] = ("AtTrVaL1", "AtTrVaL2")

        self.assertIn("attrkey", cia)
        self.assertIn(rand_case("attrkey"), cia)

        self.assertEqual(cia['attrkey'], ("AtTrVaL1", "AtTrVaL2"))
        self.assertEqual(cia[rand_case('attrkey')], ("AtTrVaL1", "AtTrVaL2"))

        self.assertEqual(cia.get('attrkey'), ("AtTrVaL1", "AtTrVaL2"))
        self.assertEqual(cia.get(rand_case('attrkey')), ("AtTrVaL1", "AtTrVaL2"))

        self.assertEqual(cia.get(rand_case('attrkey'), None), ("AtTrVaL1", "AtTrVaL2"))
        self.assertEqual(cia.get('badkey', "altval"), "altval")

        self.assertTrue(cia.contains_ident(('attrkey', 'attrval1')))
        self.assertTrue(cia.contains_ident((rand_case('attrkey'), rand_case('attrval2'))))

    """Does CaseInsensitiveAttrs correctly support iteration of original-case keys and values?"""

    def test_CaseInsensitiveAttrs_iteration(self):
        cia = CaseInsensitiveAttrs()

        cia["AtTrKeY"] = ("AtTrVaL1", "AtTrVaL2")
        cia["AtTrKeY2"] = ("AtTrVaL3",)
        cia["AtTrKeY3"] = ("AtTrVaL4", "AtTrVaL5", "AtTrVaL6")

        self.assertListEqual(list(cia.keys()), [
            "AtTrKeY",
            "AtTrKeY2",
            "AtTrKeY3"
        ])
        self.assertListEqual(list(cia.values()), [
            ("AtTrVaL1", "AtTrVaL2"),
            ("AtTrVaL3",),
            ("AtTrVaL4", "AtTrVaL5", "AtTrVaL6")
        ])
        self.assertListEqual(list(cia.items()), [
            ("AtTrKeY", ("AtTrVaL1", "AtTrVaL2")),
            ("AtTrKeY2", ("AtTrVaL3",)),
            ("AtTrKeY3", ("AtTrVaL4", "AtTrVaL5", "AtTrVaL6"))
        ])

    """Does CaseInsensitiveAttrs throw KeyError when there's a case-independent key mismatch on lookup?"""

    def test_CaseInsensitiveAttrs_KeyError(self):
        cia = CaseInsensitiveAttrs()
        cia['attrkey'] = ('N/A',)

        with self.assertRaises(KeyError):
            _ = cia['badkey']

        self.assertEqual(cia.get('badkey'), None)

    """Does CaseInsensitiveAttrs support the setdefault() method?"""

    def test_CaseInsensitiveAttrs_setdefault(self):
        cia = CaseInsensitiveAttrs()
        cia['AtTrKeY'] = ("AtTrVaL1", "AtTrVaL2")

        exists = cia.setdefault('attrkey', ('attrval',))
        dne = cia.setdefault('OtHeRkEy', ('AtTrVaL2',))

        self.assertEqual(exists, ("AtTrVaL1", "AtTrVaL2"))
        self.assertEqual(dne, ('AtTrVaL2',))
        self.assertEqual(cia['otherkey'], ('AtTrVaL2',))
        self.assertEqual(cia['attrkey'], ("AtTrVaL1", "AtTrVaL2"))

    """Does CaseInsensitiveAttrs.contains_ident() properly handle wildcard queries?"""

    def test_CaseInsensitiveAttrs_contains_ident_wildcard(self):
        cia = CaseInsensitiveAttrs()

        cia["AtTrKeY"] = ("AtTrVaL1", "AtTrVaL2")
        cia["AtTrKeY2"] = ("AtTrVaL3",)
        cia["AtTrKeY3"] = ("AtTrVaL4", "AtTrVaL5", "AtTrVaL6")

        self.assertTrue(cia.contains_ident(("attrkey3", "attrval5")))
        self.assertTrue(cia.contains_ident(("attrkey2", Wildcard())))
        self.assertTrue(cia.contains_ident((Wildcard(), "attrval6")))
        self.assertTrue(cia.contains_ident((Wildcard(), Wildcard())))

        self.assertFalse(cia.contains_ident(("attrkey4", "attrval7")))
        self.assertFalse(cia.contains_ident(("attrkey4", Wildcard())))
        self.assertFalse(cia.contains_ident((Wildcard(), "attrval7")))


@unittest.skip("Long-running test, execute manually")
class GenomeParsingTests(unittest.TestCase):
    """Runs full-scale, unmodified GFF3/GTF genomes for select species through the ReferenceFeatures class"""

    @classmethod
    def setUpClass(self):
        import requests
        release = "54"
        baseurl = "http://ftp.ensemblgenomes.org/pub/release-"
        self.data_dir = "./testdata/local_only/gff/"
        if not os.path.exists(self.data_dir): os.makedirs(self.data_dir)

        self.urls = {
            'Arabidopsis':
                baseurl + '%s/plants/{}/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.%s.{}.gz' % (release, release),
            'C. elegans':
                baseurl + '%s/metazoa/{}/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.%s.{}.gz' % (release, release)
        }

        self.genomes = {
            species: os.path.join(self.data_dir, os.path.basename(url))
            for species, url in self.urls.items()
        }

        # Download genome files if they aren't already in the data dir
        for ftype in ('gff3', 'gtf'):
            for species, file_template in self.genomes.items():
                file = file_template.format(ftype)
                if os.path.isfile(file): continue

                print(f"Downloading {ftype} genome for {species} from Ensembl...")
                url = self.urls[species].format(ftype, ftype)
                with requests.get(url, stream=True) as r:
                    r.raise_for_status()
                    with open(file, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=16384):
                            f.write(chunk)

    """Does ReferenceFeatures.get() process all genomes without throwing any errors?"""
    def test_gff_megazord(self):
        print("Running GFF Megazord test. This will take a long time...")

        # Single rule with all wildcard selectors, but only Identity is actually relevant within ReferenceFeatures
        rules = [{'Identity': ('', ''), 'Class': '', 'Filter_s': '', 'Filter_t': '', 'Hierarchy': 0,
                  'Overlap': 'partial', 'Strand': '', 'nt5end': '', 'Length': ''}]
        files = {gff.format(ftype): [] for gff in self.genomes.values() for ftype in ('gff3', 'gtf')}

        fs = FeatureSelector(rules)
        rt = ReferenceFeatures(files)

        # The test is passed if this command
        # completes without throwing errors.
        rt.get(fs)

if __name__ == '__main__':
    unittest.main()
