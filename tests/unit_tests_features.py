import unittest
from copy import deepcopy

from unittest.mock import patch, call
from tiny.rna.counter.features import *
from unit_test_helpers import read, make_parsed_sam_record, rules_template, strand_to_bool

resources = "./testdata/counter"

class FeaturesTests(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.gff_file = f"{resources}/identity_choice_test.gff3"
        self.short_gff_file = f"{resources}/single.gff3"
        self.short_gff = read(self.short_gff_file)

        self.sam_file = f"{resources}/identity_choice_test.sam"
        self.short_sam_file = f"{resources}/single.sam"
        self.short_sam = read(self.short_sam_file)

    """Helper functions"""

    def make_feature_for_interval_test(self, iv_rule, feat_id, chrom, strand: str, start, stop):
        feat_iv = HTSeq.GenomicInterval(chrom, start, stop, strand)
        rules = [dict(deepcopy(rules_template[0]), Overlap=iv_rule, Identity=('*', '*'), nt5end='*', Length='*')]
        fs = FeatureSelector(rules, LibraryStats(normalize_by_hits=True))

        # Feature with specified coordinates, matching rule 0 with hierarchy 0 and the appropriate selector for iv_rule
        selectors = fs.build_interval_selectors(feat_iv, [(0, 0, iv_rule)])
        match_tuple = (selectors[feat_iv][0],)
        feat = {(feat_id, strand_to_bool(strand), match_tuple)}

        return feat, fs

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
        self.assertEqual(matches_with_cooridnates[0][0], HTSeq.GenomicInterval("I", 9, 10, '+'))
        self.assertEqual(matches_with_cooridnates[1][0], HTSeq.GenomicInterval("I", 10, 15, '+'))
        self.assertEqual(matches_with_cooridnates[2][0], HTSeq.GenomicInterval("I", 15, 20, '+'))
        self.assertEqual(matches_with_cooridnates[2][0], HTSeq.GenomicInterval("I", 15, 20, '+'))

    """What happens if we invert the previous test's strand? Do intervals behave differently?"""

    def test_HTSeq_antisense_iv_slice(self):
        gas = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        iva = HTSeq.GenomicInterval("I", 1, 10, '-')
        ivb = HTSeq.GenomicInterval("I", 5, 15, '-')
        ivc = HTSeq.GenomicInterval("I", 9, 20, '-')
        ivd = HTSeq.GenomicInterval("I", 2, 4, "-")
        gas[iva] += "TestA"
        gas[ivb] += "TestB"
        gas[ivd] += "TestD"

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
        self.assertEqual(matches_with_cooridnates[0][0], HTSeq.GenomicInterval("I", 9, 10, '-'))
        self.assertEqual(matches_with_cooridnates[1][0], HTSeq.GenomicInterval("I", 10, 15, '-'))
        self.assertEqual(matches_with_cooridnates[2][0], HTSeq.GenomicInterval("I", 15, 20, '-'))
        self.assertEqual(matches_with_cooridnates[2][0], HTSeq.GenomicInterval("I", 15, 20, '-'))

    """Does assign_features correctly handle alignments with zero feature matches?"""

    def test_assign_features_no_match(self):
        htsgas = HTSeq.GenomicArrayOfSets("auto", stranded=False)
        Features.chrom_vectors = htsgas.chrom_vectors
        chrom, strand = "I", "+"

        # Add test feature and interval to the GenomicArray
        iv_none = HTSeq.GenomicInterval(chrom, 0, 2, strand)
        htsgas[iv_none] += "Should not match"

        # Create mock SAM alignment with non-overlapping interval
        sam_aln = make_parsed_sam_record(**{'Start': 2, 'Chrom': chrom, 'Strand': strand_to_bool(strand)})

        """
        iv_none: 0 |--| 2
        sam_aln:    2 |-- ... --|
        """

        with patch("tiny.rna.counter.features.FeatureCounter") as mock:
            instance = mock.return_value
            FeatureCounter.assign_features(instance, sam_aln)

        instance.choose.assert_not_called()
        instance.stats.chrom_misses.assert_not_called()

    """Does assign_features return features that overlap the query interval by a single base?"""

    def test_assign_features_single_base_overlap(self):
        htsgas = HTSeq.GenomicArrayOfSets("auto", stranded=False)
        Features.chrom_vectors = htsgas.chrom_vectors
        chrom, strand = "I", "+"
        
        # Add test features and intervals to the Genomic Array
        iv_olap = HTSeq.GenomicInterval(chrom, 1, 2, strand)
        iv_none = HTSeq.GenomicInterval(chrom, 2, 3, strand)
        htsgas[iv_olap] += "Overlaps alignment by one base"
        htsgas[iv_none] += "Non-overlapping feature"

        # Create mock SAM alignment which overlaps iv_olap by one base
        sam_aln = make_parsed_sam_record(**{'Chrom': chrom, 'Strand': strand_to_bool(strand), 'Start': 0, 'Seq': 'AT'})

        """
        iv_none:    2 |-| 3
        iv_olap:  1 |-| 2
        sam_aln: 0 |--| 2
        """

        with patch("tiny.rna.counter.features.FeatureCounter") as mock:
            instance = mock.return_value
            FeatureCounter.assign_features(instance, sam_aln)

        expected_match_list = {'Overlaps alignment by one base'}
        instance.selector.choose.assert_called_once_with(expected_match_list, sam_aln)
        instance.stats.chrom_misses.assert_not_called()

    """Does count_reads call the right functions when handling a single record library?"""

    @patch("tiny.rna.counter.features.FeatureCounter")
    def test_count_reads_generic(self, mock):
        instance = mock.return_value
        library = {'File': self.short_sam_file, 'Name': 'short'}
        bundle = instance.stats.count_bundle.return_value

        mock_bundle = ["mock_alignment"]
        mock_read_count = 1

        instance.sam_reader.bundle_multi_alignments.return_value = iter([(mock_bundle, mock_read_count)])
        instance.assign_features.return_value = ({"mock_feat"}, 1)

        expected_calls_to_stats = [
            call.assign_library(library),
            call.count_bundle(mock_bundle, mock_read_count),
            call.count_bundle_assignments(bundle, mock_bundle[0], {'mock_feat'}, len({'mock_feat'})),
            call.finalize_bundle(bundle)
        ]

        # CALL FUNCTION
        FeatureCounter.count_reads(instance, library)

        instance.stats.assign_library.assert_called_once()
        instance.assign_features.assert_called_once_with("mock_alignment")
        instance.stats.assert_has_calls(expected_calls_to_stats)

    """Does FeatureSelector.choose() correctly select features defining `full` interval matching rules?"""

    def test_feature_selector_full_interval(self):
        iv_rule = 'full'
        chrom, strand, start, stop = "n/a", "+", 5, 10
        feat, fs = self.make_feature_for_interval_test(iv_rule, "Full Overlap", chrom, strand, start, stop)

        aln_base = {'Seq': 'ATGC', 'Chrom': chrom, 'Strand': strand_to_bool(strand)}
        aln_spill_lo = make_parsed_sam_record(**dict(aln_base, Start=start - 1, Name="spill"))
        aln_spill_hi = make_parsed_sam_record(**dict(aln_base, Start=start + 2, Name="spill"))
        aln_contained = make_parsed_sam_record(**dict(aln_base, Seq="N", Start=7, Name="contained"))
        aln_contained_lo = make_parsed_sam_record(**dict(aln_base, Start=start, Name="contained"))
        aln_contained_hi = make_parsed_sam_record(**dict(aln_base, Start=start + 1, Name="contained"))

        """
        aln:                  |ATGC|
        feat:               5 |-----| 10
        aln_spill_lo:      4 |ATGC| 8
        aln_spill_hi:         7 |ATGC| 11
        aln_contained:        7 |N| 8       # Fully contained
        aln_contained_lo:   5 |ATGC| 9      # Shared start position
        aln_contained_hi:    6 |ATGC| 10    # Shared end position
        """

        self.assertEqual(fs.choose(feat, aln_spill_lo), {})
        self.assertEqual(fs.choose(feat, aln_spill_hi), {})
        self.assertEqual(fs.choose(feat, aln_contained), {"Full Overlap": {0}})
        self.assertEqual(fs.choose(feat, aln_contained_lo), {"Full Overlap": {0}})
        self.assertEqual(fs.choose(feat, aln_contained_hi), {"Full Overlap": {0}})

    """Does FeatureSelector.choose() correctly select features with `partial` interval matching rules?"""

    def test_feature_selector_partial_interval(self):
        iv_rule = "partial"
        chrom, strand, start, stop = "n/a", "+", 5, 10
        feat, fs = self.make_feature_for_interval_test(iv_rule, "Partial Overlap", chrom, strand, start, stop)

        aln_base = {'Seq': 'ATGC', 'Chrom': chrom, 'Strand': strand_to_bool(strand)}
        aln_spill_lo = make_parsed_sam_record(**dict(aln_base, Start=start - 1, Name="spill"))
        aln_spill_hi = make_parsed_sam_record(**dict(aln_base, Start=start + 2, Name="spill"))
        aln_contained_lo = make_parsed_sam_record(**dict(aln_base, Start=start, Name="contained"))
        aln_contained_hi = make_parsed_sam_record(**dict(aln_base, Start=start + 1, Name="contained"))

        """
        aln:                  |ATGC|
        feat:               5 |-----| 10
        aln_spill_lo:      4 |ATGC| 8
        aln_spill_hi:         7 |ATGC| 11
        aln_contained_lo:   5 |ATGC| 9      Shared start position
        aln_contained_hi:    6 |ATGC| 10    Shared end position
        """

        self.assertEqual(fs.choose(feat, aln_spill_lo), {"Partial Overlap": {0}})
        self.assertEqual(fs.choose(feat, aln_spill_hi), {"Partial Overlap": {0}})
        self.assertEqual(fs.choose(feat, aln_contained_lo), {"Partial Overlap": {0}})
        self.assertEqual(fs.choose(feat, aln_contained_hi), {"Partial Overlap": {0}})

    """Does FeatureSelector.choose() correctly select features with `exact` interval matching rules?"""

    def test_feature_selector_exact_interval(self):
        iv_rule = "exact"
        chrom, strand, start, stop = "n/a", "+", 5, 9
        feat, fs = self.make_feature_for_interval_test(iv_rule, "Exact Overlap", chrom, strand, start, stop)

        aln_match = {'Seq': 'ATGC', 'Chrom': chrom, 'Strand': strand_to_bool(strand)}
        aln_short = {'Seq': 'NNN', 'Chrom': chrom, 'Strand': strand_to_bool(strand)}
        aln_exact = make_parsed_sam_record(**dict(aln_match, Start=start, Name="exact"))
        aln_short_lo = make_parsed_sam_record(**dict(aln_short, Start=start, Name="match lo"))
        aln_short_hi = make_parsed_sam_record(**dict(aln_short, Start=start + 1, Name="match hi"))
        aln_spill_lo = make_parsed_sam_record(**dict(aln_match, Start=start - 1, Name="spill lo"))
        aln_spill_hi = make_parsed_sam_record(**dict(aln_match, Start=start + 1, Name="spill hi"))

        """
        aln_match:            |ATGC|
        aln_short:            |NNN|
        feat:               5 |----| 9
        aln_exact:          5 |ATGC| 9
        aln_short_lo:       5 |NNN| 8
        aln_short_hi:        6 |NNN| 9
        aln_spill_lo:      4 |ATGC| 8
        aln_spill_hi:        6 |ATGC| 10
        """

        self.assertEqual(fs.choose(feat, aln_exact), {"Exact Overlap": {0}})
        self.assertEqual(fs.choose(feat, aln_short_lo), {})
        self.assertEqual(fs.choose(feat, aln_short_hi), {})
        self.assertEqual(fs.choose(feat, aln_spill_lo), {})
        self.assertEqual(fs.choose(feat, aln_spill_hi), {})

    """Does FeatureSelector.choose() correctly select features with `5' anchored` interval matching rules?"""

    def test_feature_selector_5p_interval(self):
        iv_rule = "5' anchored"
        chrom, start, end = "n/a", 5, 10

        """
                No match    | 6 -------->| 10     aln_none
                   Match  5 |------------|--> 11  aln_long
                   Match  5 |----------->| 10     aln_exact
                   Match  5 |--------> 9 |        aln_short
        (+) 5' -------------|==feat_A===>|-----------> 3'
                  start = 5                end = 10
        (-) 3' <------------|<===feat_B==|------------ 5'
                            | 6 <--------| 10  Match
                          5 |<-----------| 10  Match
                       4 <--|------------| 10  Match
                          5 |<-------- 9 |     No match
        """

        # Test feat_A on (+) strand
        feat_A, fs = self.make_feature_for_interval_test(iv_rule, "5' Anchored Overlap (+)", chrom, '+', start, end)
        aln_base = {'Start': start, 'Chrom': chrom, 'Strand': True}
        aln = {
            'aln_none': make_parsed_sam_record(**dict(aln_base, Start=start + 1, Seq="ATGC")),
            'aln_long': make_parsed_sam_record(**dict(aln_base, Seq="ATGCNN")),
            'aln_exact': make_parsed_sam_record(**dict(aln_base, Seq="ATGCN")),
            'aln_short': make_parsed_sam_record(**dict(aln_base, Seq="ATGC")),
        }

        self.assertEqual(fs.choose(feat_A, aln['aln_none']), {})
        self.assertEqual(fs.choose(feat_A, aln['aln_long']), {"5' Anchored Overlap (+)": {0}})
        self.assertEqual(fs.choose(feat_A, aln['aln_exact']), {"5' Anchored Overlap (+)": {0}})
        self.assertEqual(fs.choose(feat_A, aln['aln_short']), {"5' Anchored Overlap (+)": {0}})

        # Test feat_B on (-) strand
        feat_B, fs = self.make_feature_for_interval_test(iv_rule, "5' Anchored Overlap (-)", chrom, '-', start, end)
        aln['aln_short'].update({'Start': 6, 'End': 10, 'Strand': False})
        aln['aln_exact'].update({'Start': 5, 'End': 10, 'Strand': False})
        aln['aln_long'].update({'Start': 4, 'End': 10, 'Strand': False})
        aln['aln_none'].update({'Start': 5, 'End': 9, 'Strand': False})

        self.assertEqual(fs.choose(feat_B, aln['aln_none']), {})
        self.assertEqual(fs.choose(feat_B, aln['aln_long']), {"5' Anchored Overlap (-)": {0}})
        self.assertEqual(fs.choose(feat_B, aln['aln_exact']), {"5' Anchored Overlap (-)": {0}})
        self.assertEqual(fs.choose(feat_B, aln['aln_short']), {"5' Anchored Overlap (-)": {0}})

    """Does FeatureSelector.choose() correctly select features with `3' anchored` interval matching rules?"""

    def test_feature_selector_3p_interval(self):
        iv_rule = "3' anchored"
        chrom, start, end = "n/a", 5, 10

        """
                No match  5 |--------> 9 |     aln_none
                  Match 4 --|----------->| 10  aln_long
                  Match   5 |----------->| 10  aln_exact
                  Match     | 6 -------->| 10  aln_short
        (+) 5' -------------|==feat_A===>|-----------> 3'
                  start = 5                end = 10
        (-) 3' <------------|<===feat_B==|------------ 5'
                          5 |<-------- 9 |       Match
                          5 |<-----------| 10    Match
                          5 |<-----------|-- 11  Match
                            | 6 <--------| 10    No match
        """

        # Test feat_A on (+) strand
        feat_A, fs = self.make_feature_for_interval_test(iv_rule, "3' Anchored Overlap (+)", chrom, '+', start, end)
        aln_base = {'Start': start, 'Chrom': chrom, 'Strand': True}
        aln = {
            'aln_none': make_parsed_sam_record(**dict(aln_base, Seq="ATGC")),
            'aln_long': make_parsed_sam_record(**dict(aln_base, Start=start - 1, Seq="ATGCNN")),
            'aln_exact': make_parsed_sam_record(**dict(aln_base, Seq="ATGCN")),
            'aln_short': make_parsed_sam_record(**dict(aln_base, Start=start + 1, Seq="ATGC")),
        }

        self.assertEqual(fs.choose(feat_A, aln['aln_none']), {})
        self.assertEqual(fs.choose(feat_A, aln['aln_long']), {"3' Anchored Overlap (+)": {0}})
        self.assertEqual(fs.choose(feat_A, aln['aln_exact']), {"3' Anchored Overlap (+)": {0}})
        self.assertEqual(fs.choose(feat_A, aln['aln_short']), {"3' Anchored Overlap (+)": {0}})

        # Test feat_B on (-) strand
        feat_B, fs = self.make_feature_for_interval_test(iv_rule, "3' Anchored Overlap (-)", chrom, '-', start, end)
        aln['aln_short'].update({'Start': 5, 'End': 9, 'Strand': False})
        aln['aln_exact'].update({'Start': 5, 'End': 10, 'Strand': False})
        aln['aln_long'].update({'Start': 5, 'End': 11, 'Strand': False})
        aln['aln_none'].update({'Start': 6, 'End': 10, 'Strand': False})

        self.assertEqual(fs.choose(feat_B, aln['aln_none']), {})
        self.assertEqual(fs.choose(feat_B, aln['aln_long']), {"3' Anchored Overlap (-)": {0}})
        self.assertEqual(fs.choose(feat_B, aln['aln_exact']), {"3' Anchored Overlap (-)": {0}})
        self.assertEqual(fs.choose(feat_B, aln['aln_short']), {"3' Anchored Overlap (-)": {0}})

    """Are wildcard keywords in identity selectors properly converted to a Wildcard object copy?"""

    def test_wildcard_identities(self):
        wildcard, non_wild = "all", "Non-wildcard"
        one = [dict(deepcopy(rules_template[0]), Identity=(wildcard, non_wild)),
               dict(deepcopy(rules_template[0]), Identity=(non_wild, wildcard))]
        two = [dict(deepcopy(rules_template[0]), Identity=(wildcard, wildcard))]
        non = [dict(deepcopy(rules_template[0]), Identity=(non_wild, non_wild))]
        dup = deepcopy(two)

        rules = [*one, *two, *non, *dup]

        actual = FeatureSelector(rules, LibraryStats(normalize_by_hits=True)).inv_ident
        expected = {
            (Wildcard(), non_wild):   [0],
            (non_wild,   Wildcard()): [1],
            (Wildcard(), Wildcard()): [2, 4],
            (non_wild,   non_wild):   [3]
        }

        self.assertDictEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
