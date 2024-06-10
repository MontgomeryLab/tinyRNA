import collections
import contextlib
import unittest
import HTSeq
import io

from copy import deepcopy
from random import randint
from unittest.mock import patch, mock_open, call

from tiny.rna.counter.features import FeatureSelector
from tiny.rna.counter.matching import *
from tiny.rna.counter.hts_parsing import *
from tiny.rna.counter.parsing.alignments import _validate_alignment

import unit_test_helpers as helpers

resources = "./testdata/counter"
wc = Wildcard()  # used in many tests as a 'n/a' value

# To run all test suites
if __name__ == '__main__':
    unittest.main()


class AlignmentReaderTests(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.sam_file = f"{resources}/sam/identity_choice_test.sam"
        self.short_sam_file = f"{resources}/sam/single.sam"
        self.empty_sam_file = f"{resources}/sam/empty.sam"
        self.short_sam = helpers.read(self.short_sam_file)
        self.short_bam_file = f"{resources}/bam/single.bam"

    @staticmethod
    def exhaust_iterator(it):
        collections.deque(it, maxlen=0)

    # === TESTS ===

    """Did AlignmentReader correctly skip header values and parse all pertinent info from a single record SAM file?"""

    def test_AlignmentReader_single_sam(self):
        sam_bundle, read_count = next(AlignmentReader().bundle_multi_alignments(self.short_sam_file))
        sam_record = sam_bundle[0]

        self.assertEqual(sam_record['Chrom'], "I")
        self.assertEqual(sam_record['Start'], 15064569)
        self.assertEqual(sam_record['End'], 15064590)
        self.assertEqual(sam_record['Strand'], False)
        self.assertEqual(sam_record['Name'], "0_count=5")
        self.assertEqual(sam_record['Seq'], "CAAGACAGAGCTTCACCGTTC")
        self.assertEqual(sam_record['Length'], 21)
        self.assertEqual(sam_record['nt5end'], 'G')

    """Did AlignmentReader correctly skip header values and parse all pertinent info from a single record BAM file?"""

    def test_AlignmentReader_single_bam(self):
        bam_bundle, read_count = next(AlignmentReader().bundle_multi_alignments(self.short_bam_file))
        bam_record = bam_bundle[0]

        self.assertEqual(bam_record['Chrom'], "I")
        self.assertEqual(bam_record['Start'], 15064569)
        self.assertEqual(bam_record['End'], 15064590)
        self.assertEqual(bam_record['Strand'], False)
        self.assertEqual(bam_record['Name'], "0_count=5")
        self.assertEqual(bam_record['Seq'], "CAAGACAGAGCTTCACCGTTC")
        self.assertEqual(bam_record['Length'], 21)
        self.assertEqual(bam_record['nt5end'], 'G')

    """Does our AlignmentReader produce the same pertinent info from a SAM file as HTSeq's BAM_reader?

    A note on SAM files: reads are always stored 5' to 3', so antisense reads are actually
    recorded in reverse complement. HTSeq automatically performs this conversion, but we
    are only really concerned about a sequence's 5' end NT, so our alignment dicts performs
    this conversion more surgically for only the 5' end NT at construction time.
    """

    def test_sam_parser_comparison(self):
        file = f"{resources}/sam/Lib304_test.sam"
        ours = AlignmentReader().bundle_multi_alignments(file)
        theirs = HTSeq.bundle_multiple_alignments(HTSeq.BAM_Reader(file))

        for (our_bundle, _), their_bundle in zip(ours, theirs):
            self.assertEqual(len(our_bundle), len(their_bundle))
            for our, their in zip(our_bundle, their_bundle):
                self.assertEqual(our['Chrom'], their.iv.chrom)
                self.assertEqual(our['Start'], their.iv.start)
                self.assertEqual(our['End'], their.iv.end)
                self.assertEqual(our['Name'], their.read.name)
                self.assertEqual(our['nt5end'], chr(their.read.seq[0]))  # See note above
                self.assertEqual(our['Strand'], helpers.strand_to_bool(their.iv.strand))
                if our['Strand'] is False:  # See note above
                    self.assertEqual(our['Seq'][::-1].translate(helpers.complement), their.read.seq.decode())
                else:
                    self.assertEqual(our['Seq'], their.read.seq.decode())

    """Does our AlignmentReader produce the same pertinent info from a BAM file as HTSeq's BAM_reader?"""

    def test_bam_parser_comparison(self):
        file = f"{resources}/bam/Lib304_test.bam"
        ours = AlignmentReader().bundle_multi_alignments(file)
        theirs = HTSeq.bundle_multiple_alignments(HTSeq.BAM_Reader(file))

        for (our_bundle, _), their_bundle in zip(ours, theirs):
            self.assertEqual(len(our_bundle), len(their_bundle))
            for our, their in zip(our_bundle, their_bundle):
                self.assertEqual(our['Chrom'], their.iv.chrom)
                self.assertEqual(our['Start'], their.iv.start)
                self.assertEqual(our['End'], their.iv.end)
                self.assertEqual(our['Name'], their.read.name)
                self.assertEqual(our['nt5end'], chr(their.read.seq[0]))  # See note above
                self.assertEqual(our['Strand'], helpers.strand_to_bool(their.iv.strand))
                if our['Strand'] is False:  # See note above
                    self.assertEqual(our['Seq'][::-1].translate(helpers.complement), their.read.seq.decode())
                else:
                    self.assertEqual(our['Seq'], their.read.seq.decode())

    """Does AlignmentReader._get_decollapsed_filename() create an appropriate filename?"""

    def test_AlignmentReader_get_decollapsed_filename(self):
        reader = AlignmentReader()
        reader.file = "~/path/to/input/sam_file.sam"

        sam_out = reader._get_decollapsed_filename()

        self.assertEqual(sam_out, "sam_file_decollapsed.sam")

    """Does AlignmentReader._new_bundle report the correct read count for different collapser types?"""

    def test_AlignmentReader_new_bundle(self):
        qnames = ["0_count=3", "seq0_x5", "non-collapsed"]
        counts = [3, 5, 1]

        reader = AlignmentReader()
        reader._header_dict = {'HD': {'SO': 'queryname'}}

        for qname, expected in zip(qnames, counts):
            reader._determine_collapser_type(qname)
            _, read_count = reader._new_bundle({'Name': qname})
            self.assertEqual(read_count, expected)

    """Does AlignmentReader._gather_metadata() correctly report SAM files that lack alignments?"""

    def test_AlignmentReader_empty_sam(self):
        reader = AlignmentReader()
        reader._assign_library(self.empty_sam_file)
        sam_in = pysam.AlignmentFile(self.empty_sam_file)

        message = rf"Alignment file is empty \({os.path.basename(self.empty_sam_file)}\)\."
        with self.assertRaisesRegex(ValueError, message):
            reader._gather_metadata(sam_in)

    """Does AlignmentReader._gather_metadata() correctly identify metadata and write the decollapsed file header?"""

    def test_AlignmentReader_gather_metadata(self):
        reader = AlignmentReader(decollapse=True)
        reader._decollapsed_filename = "mock_outfile_name.sam"
        reader._assign_library(self.short_sam_file)
        sam_in = pysam.AlignmentFile(self.short_sam_file)

        with patch('builtins.open', mock_open()) as sam_out:
            reader._gather_metadata(sam_in)

        expected_writelines = [
            call('mock_outfile_name.sam', 'w'),
            call().__enter__(),
            call().write("@HD\tSO:unsorted\n@SQ\tSN:I\tLN:21\n@PG\tID:bowtie\n"),
            call().__exit__(None, None, None)
        ]

        expected_header_dict = {
            'HD': {'SO': 'unsorted'},
            'SQ': [{'LN': 21, 'SN': 'I'}],
            'PG': [{'ID': 'bowtie'}]
        }

        sam_out.assert_has_calls(expected_writelines)
        self.assertEqual(reader.collapser_type, 'tiny-collapse')
        self.assertDictEqual(reader._header_dict, expected_header_dict)
        self.assertEqual(reader.references, ('I',))
        self.assertIn("NM", reader.expected_tags)

    """Does AlignmentReader._write_decollapsed_sam() write the correct number of duplicates to the decollapsed file?"""

    def test_AlignmentReader_write_decollapsed_sam(self):
        header = pysam.AlignmentHeader()
        alignment_in = pysam.AlignedSegment(header)
        alignment_in.query_name = "0_count=5"
        alignment_out = pysam.AlignedSegment(header)
        alignment_out.query_name = "0_count"

        reader = AlignmentReader(decollapse=True)
        reader.collapser_type = "tiny-collapse"
        reader._collapser_token = "="
        reader._decollapsed_reads = [alignment_in]
        reader._decollapsed_filename = "mock_outfile_name.sam"

        expected_writelines = [
            call('mock_outfile_name.sam', 'a'),
            call().__enter__(),
            call().writelines([alignment_out.to_string() + '\n'] * 5),
            call().__exit__(None, None, None)
        ]

        with patch('builtins.open', mock_open()) as outfile:
            reader._write_decollapsed_sam()

        outfile.assert_has_calls(expected_writelines)
        self.assertTrue(len(reader._decollapsed_reads) == 0)

    """Does AlignmentReader._parse_alignments() save lines and write them to the decollapsed file when appropriate?"""

    def test_AlignmentReader_parse_alignments_decollapse(self):
        with patch.object(AlignmentReader, "_write_decollapsed_sam") as write_fn:
            # Set up AlignmentReader class
            reader = AlignmentReader(decollapse=True)
            reader.collapser_type = "tiny-collapse"
            reader._decollapsed_reads = buffer = [0] * 99999  # At 100,001, expect buffer to be written
            reader.file = self.short_sam_file                 # File with single alignment

            # Set up AlignmentIter class
            sam_in = pysam.AlignmentFile(reader.file)
            callback = reader._write_decollapsed_sam
            first_aln_offset = sam_in.tell()
            expected_tags = ("NM",)
            aln_iter = AlignmentIter(sam_in, expected_tags, callback, buffer)

            # Add 100,000th alignment to the buffer
            self.exhaust_iterator(aln_iter)
            write_fn.assert_not_called()

            # Rewind and add one more alignment to push it over threshold
            sam_in.seek(first_aln_offset)
            self.exhaust_iterator(aln_iter)
            write_fn.assert_called_once()

    """Are decollapsed outputs skipped when non-collapsed SAM files are supplied?"""

    def test_AlignmentReader_no_decollapse_non_collapsed_SAM_files(self):
        stdout_capture = io.StringIO()
        with patch.object(AlignmentReader, "_write_decollapsed_sam") as write_sam, \
                patch.object(AlignmentReader, "_write_header_for_decollapsed_sam") as write_header:
            with contextlib.redirect_stderr(stdout_capture):
                reader = AlignmentReader(decollapse=True)
                records = reader.bundle_multi_alignments(f"{resources}/sam/non-collapsed.sam")
                self.exhaust_iterator(records)

        write_sam.assert_not_called()
        write_header.assert_not_called()
        self.assertEqual(reader.collapser_type, None)
        self.assertEqual(stdout_capture.getvalue(),
                         "Alignments do not appear to be derived from a supported collapser input. "
                         "Decollapsed SAM files will therefore not be produced.\n")

    """Is incompatible alignment file ordering correctly identified from @HD header values?"""

    def test_AlignmentReader_incompatible_HD_header(self):
        # Valid SO values: ["unknown", "unsorted", "queryname", "coordinate"]
        # Valid GO values: ["none", "query", "reference"]
        
        strictly_compatible = [
            {'HD': {'SO': "queryname"}},
            {'HD': {'GO': "query"}},
        ]

        # Should not throw error
        for header in strictly_compatible:
            reader = AlignmentReader()
            reader._header_dict = header
            reader._assign_library("mock_infile.sam")
            reader._check_for_incompatible_order()

        strictly_incompatible = [
            ({'HD': {'SO': "coordinate"}}, "by coordinate"),
            ({'HD': {'GO': "reference"}},  "by reference"),
            ({'HD': {}},                   "sorting/grouping couldn't be determined"),
            ({},                           "sorting/grouping couldn't be determined"),
        ]

        # Should throw error
        for header, message in strictly_incompatible:
            reader = AlignmentReader()
            reader._header_dict = header
            reader._assign_library("mock_infile.sam")
            with self.assertRaisesRegex(ValueError, message):
                reader._check_for_incompatible_order()

    """Is incompatible alignment file ordering correctly identified from @PG header values?"""

    def test_AlignmentReader_incompatible_PG_header(self):
        # SO = ["unknown", "unsorted", "queryname", "coordinate"]
        # GO = ["none", "query", "reference"]

        self.assertEqual(AlignmentReader.compatible_unordered, ("bowtie", "bowtie2", "star"))
        SO_unordered = [{'SO': "unknown"}, {'SO': "unsorted"}]
        GO_unordered = [{'GO': "none"}]
        compatible = [
            {'PG': [{'ID': tool}], 'HD': un_so, 'GO': un_go}
            for tool in AlignmentReader.compatible_unordered
            for un_so in SO_unordered
            for un_go in GO_unordered
        ]

        # Only the last reported PG ID matters
        multiple_tools = [{'ID': "INCOMPATIBLE"}, {'ID': "bowtie"}]
        compatible += [{'PG': multiple_tools, 'HD': {'SO': "unsorted"}}]

        # Should not throw error
        for header in compatible:
            reader = AlignmentReader()
            reader._header_dict = header
            reader._assign_library("mock_infile.sam")
            reader._check_for_incompatible_order()

        incompatible = [
            {'PG': [{'ID': "INCOMPATIBLE"}], 'HD': un_so, 'GO': un_go}
            for un_so in SO_unordered
            for un_go in GO_unordered
        ]

        # Should throw error
        expected_error = "adjacency couldn't be determined"
        for header in incompatible:
            reader = AlignmentReader()
            reader._header_dict = header
            reader._assign_library("mock_infile.sam")
            with self.assertRaisesRegex(ValueError, expected_error):
                reader._check_for_incompatible_order()

    """Does AlignmentIter reject alignments that contain unsupported CIGAR operators?"""

    def test_AlignmentIter_incompatible_cigar_ops(self):
        good = ["1M", "1D", "1I", "1=", "1X"]
        bad = ["1N", "1S", "1H", "1P"]

        header = pysam.AlignmentHeader()
        test_aln = pysam.AlignedSegment(header)

        for cigar in good:
            test_aln.cigarstring = cigar
            _validate_alignment(test_aln)

        for cigar in bad:
            test_aln.cigarstring = cigar
            with self.assertRaisesRegex(ValueError, "not supported at this time"):
                _validate_alignment(test_aln)


class ReferenceFeaturesTests(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.gff_file = f"{resources}/gff/identity_choice_test.gff3"
        self.short_gff_file = f"{resources}/gff/single.gff3"
        self.short_gff = helpers.read(self.short_gff_file)

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

        self.assertEqual((type(feats), type(alias), type(tags)), (HTSeq.GenomicArray, dict, defaultdict))
        self.assertEqual(steps, [{(("Gene:WBGene00023193", 'tag'), False, ((1, 2, IntervalPartialMatch(iv), wc),))}])
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

        self.assertEqual((type(feats), type(alias), type(tags)), (HTSeq.GenomicArray, dict, defaultdict))
        self.assertEqual(steps, [{(("Gene:WBGene00023193", 'tag'), False, ((1, 2, IntervalPartialMatch(iv), wc),))}])
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
            f"{resources}/gff/single2.gff3": ["ID"]
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
            {(('Gene:WBGene00023193', ''), False, ((0, 1, ivm, wc), (1, 2, ivm, wc), (2, 3, ivm, wc)))},
            set()
        ]

        feats, _, _ = ReferenceFeatures(feature_source, **kwargs).get(feature_selector)

        actual_matches = list(feats.chrom_vectors['I']['.'].array.get_steps(values_only=True))
        self.assertListEqual(actual_matches, expected_matches)

    """Does ReferenceFeatures.get() properly handle aliases for discontinuous features?"""

    def test_ref_tables_discontinuous_aliases(self):
        kwargs = {'all_features': True}
        feature_source = {f"{resources}/gff/discontinuous.gff3": ["Name"]}
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
        feature_source = {f"{resources}/gff/discontinuous.gff3": ["Name"]}
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
        feature_source = {f"{resources}/gff/discontinuous.gff3": ["Name"]}
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
        feature_source = {f"{resources}/gff/discontinuous.gff3": ["Name"]}
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
        feature_source = {f"{resources}/gff/discontinuous.gff3": ["Name"]}
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

        rule1_gp =  {f"{iv.start}:{iv.end}": ((0, 2, IntervalPartialMatch(iv), wc),) for iv in gp_ivs}
        rule1_p2 =  {f"{iv.start}:{iv.end}": ((0, 2, IntervalPartialMatch(iv), wc),) for iv in p2_ivs}
        rule2_sib = {f"{iv.start}:{iv.end}": (1, 3, IntervalPartialMatch(iv), wc)    for iv in sib_ivs}
        rule3_sib = {f"{iv.start}:{iv.end}": (2, 0, IntervalPartialMatch(iv), wc)    for iv in sib_ivs}

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
        feature_source = {f"{resources}/gff/discontinuous.gff3": ["Name"]}
        feature_selector = self.selector_with_template([{'Filter_s': "Source2Name"}])

        child2_iv =     HTSeq.GenomicInterval('I', 39, 50, '-')
        exp_alias =     {'Child2': ('Child2Name',)}
        exp_feats =     [set(), {(('Child2', ''), False, ((0, 0, IntervalPartialMatch(child2_iv), wc),))}, set()]
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
        feature_source = {f"{resources}/gff/discontinuous.gff3": ["Name"]}
        feature_selector = self.selector_with_template([{'Filter_t': "CDS"}])

        child1_iv =     HTSeq.GenomicInterval('I', 29, 40, '-')
        exp_alias =     {'Child1': ('SharedName',)}
        exp_feats =     [set(), {(('Child1', ''), False, ((0, 0, IntervalPartialMatch(child1_iv), wc),))}, set()]
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
        feature_source = {f"{resources}/gff/discontinuous.gff3": ["Name"]}
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
        feature_source = {f"{resources}/gff/single.gff3": ["sequence_name"]}
        feature_selector = self.selector_with_template([
            {'Identity': ("ID", feat_id), 'Class': "tagged_match", 'Hierarchy': 1},
            {'Identity': ("ID", feat_id), 'Class': "",             'Hierarchy': 2}
        ])

        expected_tags = {feat_id: {(feat_id, "tagged_match"), (feat_id, '')}}
        expected_aliases = {feat_id: ('Y74C9A.6',)}
        iv = IntervalPartialMatch(HTSeq.GenomicInterval('n/a', 3746, 3909))
        expected_feats = [
            set(), {
                ((feat_id, 'tagged_match'), False, ((0, 1, iv, wc),)),
                ((feat_id, ''),             False, ((1, 2, iv, wc),))
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
        feature_source = {f"{resources}/gff/discontinuous.gff3": ['Name']}

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
                (('Parent2', 'shared'), False, ((0, 1, Parent2_iv, wc), (1, 2, Parent2_iv, wc))),
                (('Parent2', ''),       False, ((2, 3, Parent2_iv, wc),)),
            },
            set(), {
                (('Parent2', 'shared'), False, ((0, 1, Child2_iv, wc), (1, 2, Child2_iv, wc))),
                (('Parent2', ''),       False, ((2, 3, Child2_iv, wc),))
            },
            set()
        ]

        feats, aliases, tags = ReferenceFeatures(feature_source).get(feature_selector)

        stepvec = list(feats.chrom_vectors['I']['.'].array.get_steps(values_only=True))
        self.assertListEqual(stepvec, expected_feats)
        self.assertDictEqual(aliases, expected_aliases)
        self.assertDictEqual(tags, expected_tags)


class ReferenceSequencesTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        cls.maxDiff = None

    # === HELPERS ===

    @staticmethod
    def ReferenceSeqs_add(seq_id, seq_len, matches):
        rs = ReferenceSeqs({seq_id, seq_len})
        rs.selector = FeatureSelector([])
        rs.add_reference_seq(seq_id, seq_len, matches)
        return rs

    @staticmethod
    def get_steps(rs, chrom):
        return {(step[0], step[1]): step[2]
                  for step in rs.feats[chrom]['.'].array.get_steps()}

    # === TESTS ===

    """Does add_reference_seq() produce the expected GenomicArray for an untagged rule matching a single feature?"""

    def test_add_reference_seq_single(self):
        seq_id = "seq"
        seq_len = 10
        matches = {'': [(0, 0, "partial", wc)]}

        rs = self.ReferenceSeqs_add(seq_id, seq_len, matches)
        actual = self.get_steps(rs, seq_id)

        # The record_tuples should have intervals on both strands
        # and the overlap selector should have the same interval.
        # For these selectors, same interval on both strands.
        iv = HTSeq.GenomicInterval(seq_id, 0, seq_len)
        match_tuple = ((0, 0, IntervalPartialMatch(iv), wc),)

        expected = {
            (0, 10): {
                ((seq_id, ''), True,  match_tuple), # sense
                ((seq_id, ''), False, match_tuple)  # antisense
            },
            (10, sys.maxsize): set()
        }

        self.assertDictEqual(actual, expected)

    """Does add_reference_seq() produce the expected GenomicArray when rules share a classifier?"""

    def test_add_reference_seq_shared_classifier(self):
        seq_id = "seq"
        seq_len = 10
        matches = {'shared': [(0, 0, "partial", wc), (1, 1, "nested", wc)]}

        rs = self.ReferenceSeqs_add(seq_id, seq_len, matches)
        actual = self.get_steps(rs, seq_id)

        # The record_tuples will have intervals on both strands
        # and the overlap selector should have the same interval.
        # For these selectors, same interval on both strands.
        iv = HTSeq.GenomicInterval(seq_id, 0, seq_len)
        match_tuples = (0, 0, IntervalPartialMatch(iv), wc), (1, 1, IntervalNestedMatch(iv), wc)

        expected = {
            (0, 10): {
                ((seq_id, "shared"), True,  match_tuples),  # sense
                ((seq_id, "shared"), False, match_tuples)   # antisense
            },
            (10, sys.maxsize): set()
        }

        self.assertDictEqual(actual, expected)

    """Does add_reference_seq() produce the expected GenomicArray when shift parameters
    from different rules produce the same interval?"""

    def test_add_reference_seq_shared_iv(self):
        seq_id = "seq"
        seq_len = 10
        matches = {'exact': [(0, 0, "exact, 2, -2", wc)], 'nested': [(1, 1, "nested, 2, -2", wc)]}

        rs = self.ReferenceSeqs_add(seq_id, seq_len, matches)
        actual = self.get_steps(rs, seq_id)

        # The record_tuples will have intervals on both strands
        # and the overlap selector should have the same interval.
        # For these selectors, same interval on both strands.
        iv = HTSeq.GenomicInterval(seq_id, 2, seq_len - 2)
        match_exact =  ((0, 0, IntervalExactMatch(iv), wc),)
        match_nested = ((1, 1, IntervalNestedMatch(iv), wc),)

        expected = {
            (0, 2): set(),
            (2, 8): {
                ((seq_id, "exact"),  True,  match_exact),   # sense
                ((seq_id, "exact"),  False, match_exact),   # antisense
                ((seq_id, "nested"), True,  match_nested),  # sense
                ((seq_id, "nested"), False, match_nested)   # antisense
            },
            (8, sys.maxsize): set()
        }

        self.assertDictEqual(actual, expected)
    
    """Does add_reference_seq() produce the expected GenomicArray when there are
    multiple tagged rules with different shift intervals, including asymmetric shifts 
    and an overlap selector that produces different intervals on each strand?"""

    def test_add_reference_seq_complex(self):
        seq_id = "seq"
        seq_len = 10
        matches = {
            "class1": [(0, 0, "nested, 1, -1", wc), (0, 0, "exact, 5, 0", wc)],
            "class2": [(0, 0, "5' anchored, 5, 0", wc)]
        }

        rs = self.ReferenceSeqs_add(seq_id, seq_len, matches)
        actual = self.get_steps(rs, seq_id)

        # Since nested shift is symmetric, iv is same on both strands
        iv_n = HTSeq.GenomicInterval(seq_id, 1, seq_len-1)
        match_nested = ((0, 0, IntervalNestedMatch(iv_n), wc),)

        # Since exact and anchored shift is asymmetric and by the same
        # amount, iv differs per strand but is shared by both selectors.
        # If they both had the same classifier then these match tuples
        # would share the same record tuple on both strands
        iv_e5_s = HTSeq.GenomicInterval(seq_id, 5, seq_len, '+')
        match_exact_sense = ((0, 0, IntervalExactMatch(iv_e5_s), wc),)
        match_5anch_sense = ((0, 0, Interval5pMatch(iv_e5_s), wc),)

        iv_e5_a = HTSeq.GenomicInterval(seq_id, 0, seq_len-5, '-')
        match_exact_antis = ((0, 0, IntervalExactMatch(iv_e5_a), wc),)
        match_5anch_antis = ((0, 0, Interval5pMatch(iv_e5_a), wc),)

        expected = {
            (0, 1): {
                (('seq', 'class2'), False, match_5anch_antis),
                (('seq', 'class1'), False, match_exact_antis),
            },
            (1, 5): {
                (('seq', 'class1'), True,  match_nested),
                (('seq', 'class2'), False, match_5anch_antis),
                (('seq', 'class1'), False, match_exact_antis),
                (('seq', 'class1'), False, match_nested)
            },
            (5, 9): {
                (('seq', 'class1'), True,  match_exact_sense),
                (('seq', 'class1'), True,  match_nested),
                (('seq', 'class2'), True,  match_5anch_sense),
                (('seq', 'class1'), False, match_nested)
            },
            (9, 10): {
                (('seq', 'class1'), True,  match_exact_sense),
                (('seq', 'class2'), True,  match_5anch_sense)
            },
            (10, sys.maxsize): set()
        }
        
        self.assertDictEqual(actual, expected)


class CaseInsensitiveAttrsTests(unittest.TestCase):

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
