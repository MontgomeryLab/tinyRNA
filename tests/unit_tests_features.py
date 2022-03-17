import unittest

from unittest.mock import patch, call
from tiny.rna.counter.features import *
from unit_test_helpers import read, make_parsed_sam_record, rules_template

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
        sam_aln = make_parsed_sam_record(**{'start': 2, 'chrom': chrom, 'strand': strand})

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
        sam_aln = make_parsed_sam_record(**{'chrom': chrom, 'strand': strand, 'start': 0, 'seq': 'AT'})

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

        # Note: mock_alignment is wrapped in a double list to imitate multi-alignment bundling
        # It imitates a list of multi-alignment bundles with one bundle having one alignment
        instance.sam_reader.bundle_multi_alignments.return_value = [["mock_alignment"]]
        bundle = instance.stats.count_bundle.return_value
        instance.assign_features.return_value = ({"mock_feat"}, 1)

        expected_calls_to_stats = [
            call.assign_library(library),
            call.count_bundle(["mock_alignment"]),
            call.count_bundle_assignments(bundle, "mock_alignment", {'mock_feat'}, len({'mock_feat'})),
            call.finalize_bundle(bundle),
            call.diags.write_intermediate_file(library["Name"])
        ]

        # CALL FUNCTION
        FeatureCounter.count_reads(instance, library)

        instance.stats.assign_library.assert_called_once()
        instance.assign_features.assert_called_once_with("mock_alignment")
        instance.stats.assert_has_calls(expected_calls_to_stats)

    """Does FeatureSelector build the proper interval selectors?"""

    def test_feature_selector_interval_build(self):
        fs = FeatureSelector(rules_template, LibraryStats())
        iv = HTSeq.GenomicInterval('I', 0, 10, '+')

        # Match tuples formed during GFF parsing
        match_tuples = [('n/a', 'n/a', 'partial'),
                        ('n/a', 'n/a', 'full'),
                        ('n/a', 'n/a', 'exact'),
                        ('n/a', 'n/a', "5' anchored"),
                        ('n/a', 'n/a', "3' anchored")]

        result = fs.build_interval_selectors(iv, match_tuples)

        self.assertIsInstance(result[0][2], IntervalPartialMatch)
        self.assertIsInstance(result[1][2], IntervalFullMatch)
        self.assertIsInstance(result[2][2], IntervalExactMatch)
        self.assertIsInstance(result[3][2], Interval5pMatch)
        self.assertIsInstance(result[4][2], Interval3pMatch)

    """Helper function."""

    def make_feature_for_interval_test(self, iv_rule, feat_id, chrom, strand, start, stop):
        feat_iv = HTSeq.GenomicInterval(chrom, start, stop, strand)
        rules = [dict(rules_template[0], Strict=iv_rule, Identity=('N/A', 'N/A'), nt5end='all', Length='all')]
        fs = FeatureSelector(rules, LibraryStats())

        # Feature with specified coordinates, matching rule 0 with hierarchy 0 and the appropriate selector for iv_rule
        match_tuple = tuple(fs.build_interval_selectors(feat_iv, [(0, 0, iv_rule)]))
        feat = {(feat_id, strand, match_tuple)}

        return feat, fs

    """Does FeatureSelector.choose() correctly select features defining `full` interval matching rules?"""

    def test_feature_selector_full_interval(self):
        iv_rule = 'full'
        chrom, strand, start, stop = "n/a", ".", 5, 10
        feat, fs = self.make_feature_for_interval_test(iv_rule, "Strict Overlap", chrom, strand, start, stop)

        aln_base = {'seq': 'ATGC', 'chrom': chrom, 'strand': strand}
        aln_spill_lo = make_parsed_sam_record(**dict(aln_base, start=start - 1, name="spill"))
        aln_spill_hi = make_parsed_sam_record(**dict(aln_base, start=start + 2, name="spill"))
        aln_contained = make_parsed_sam_record(**dict(aln_base, seq="N", start=7, name="contained"))
        aln_contained_lo = make_parsed_sam_record(**dict(aln_base, start=start, name="contained"))
        aln_contained_hi = make_parsed_sam_record(**dict(aln_base, start=start + 1, name="contained"))

        """
        aln:                  |ATGC|
        feat:               5 |-----| 10
        aln_spill_lo:      4 |ATGC| 8
        aln_spill_hi:         7 |ATGC| 11
        aln_contained:        7 |N| 8       # Fully contained
        aln_contained_lo:   5 |ATGC| 9      # Shared start position
        aln_contained_hi:    6 |ATGC| 10    # Shared end position
        """

        self.assertEqual(fs.choose(feat, aln_spill_lo), set())
        self.assertEqual(fs.choose(feat, aln_spill_hi), set())
        self.assertEqual(fs.choose(feat, aln_contained), {"Strict Overlap"})
        self.assertEqual(fs.choose(feat, aln_contained_lo), {"Strict Overlap"})
        self.assertEqual(fs.choose(feat, aln_contained_hi), {"Strict Overlap"})

    """Does FeatureSelector.choose() correctly select features with `partial` interval matching rules?"""

    def test_feature_selector_partial_interval(self):
        iv_rule = "partial"
        chrom, strand, start, stop = "n/a", ".", 5, 10
        feat, fs = self.make_feature_for_interval_test(iv_rule, "Partial Overlap", chrom, strand, start, stop)

        aln_base = {'seq': 'ATGC', 'chrom': chrom, 'strand': strand}
        aln_spill_lo = make_parsed_sam_record(**dict(aln_base, start=start - 1, name="spill"))
        aln_spill_hi = make_parsed_sam_record(**dict(aln_base, start=start + 2, name="spill"))
        aln_contained_lo = make_parsed_sam_record(**dict(aln_base, start=start, name="contained"))
        aln_contained_hi = make_parsed_sam_record(**dict(aln_base, start=start + 1, name="contained"))

        """
        aln:                  |ATGC|
        feat:               5 |-----| 10
        aln_spill_lo:      4 |ATGC| 8
        aln_spill_hi:         7 |ATGC| 11
        aln_contained_lo:   5 |ATGC| 9      Shared start position
        aln_contained_hi:    6 |ATGC| 10    Shared end position
        """

        self.assertEqual(fs.choose(feat, aln_spill_lo), {"Partial Overlap"})
        self.assertEqual(fs.choose(feat, aln_spill_hi), {"Partial Overlap"})
        self.assertEqual(fs.choose(feat, aln_contained_lo), {"Partial Overlap"})
        self.assertEqual(fs.choose(feat, aln_contained_hi), {"Partial Overlap"})

    """Does FeatureSelector.choose() correctly select features with `exact` interval matching rules?"""

    def test_feature_selector_exact_interval(self):
        iv_rule = "exact"
        chrom, strand, start, stop = "n/a", ".", 5, 9
        feat, fs = self.make_feature_for_interval_test(iv_rule, "Exact Overlap", chrom, strand, start, stop)

        aln_match = {'seq': 'ATGC', 'chrom': chrom, 'strand': strand}
        aln_short = {'seq': 'NNN', 'chrom': chrom, 'strand': strand}
        aln_exact = make_parsed_sam_record(**dict(aln_match, start=start, name="exact"))
        aln_short_lo = make_parsed_sam_record(**dict(aln_short, start=start, name="match lo"))
        aln_short_hi = make_parsed_sam_record(**dict(aln_short, start=start + 1, name="match hi"))
        aln_spill_lo = make_parsed_sam_record(**dict(aln_match, start=start - 1, name="spill lo"))
        aln_spill_hi = make_parsed_sam_record(**dict(aln_match, start=start + 1, name="spill hi"))

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

        self.assertEqual(fs.choose(feat, aln_exact), {"Exact Overlap"})
        self.assertEqual(fs.choose(feat, aln_short_lo), set())
        self.assertEqual(fs.choose(feat, aln_short_hi), set())
        self.assertEqual(fs.choose(feat, aln_spill_lo), set())
        self.assertEqual(fs.choose(feat, aln_spill_hi), set())

    """Does FeatureSelector.choose() correctly select features with `5' anchored` interval matching rules?"""

    def test_feature_selector_5p_interval(self):
        iv_rule = "5' anchored"
        chrom, strand, start, stop = "n/a", ".", 5, 10
        feat, fs = self.make_feature_for_interval_test(iv_rule, "5' Anchored Overlap", chrom, strand, start, stop)

        aln_short = {'seq': 'ATGC', 'chrom': chrom, 'strand': strand}
        aln_long = {'seq': 'ATGCNN', 'chrom': chrom, 'strand': strand}
        aln_match_short = make_parsed_sam_record(**dict(aln_short, start=start, name="match, short"))
        aln_match_long = make_parsed_sam_record(**dict(aln_long, start=start, name="match, long"))
        aln_spill_lo_short = make_parsed_sam_record(**dict(aln_short, start=start - 1, name="spill lo, short"))
        aln_spill_lo_long = make_parsed_sam_record(**dict(aln_long, start=start - 1, name="spill lo, long"))

        """
        aln_short:              |ATGC|
        aln_long:               |ATGCNN|
        feat:                 5 |-----| 10
        aln_match_short:      5 |ATGC| 9
        aln_match_long:       5 |ATGCNN| 10
        aln_spill_lo_short:  4 |ATGC| 8
        aln_spill_lo_long:   4 |ATGCNN| 10
        """

        self.assertEqual(fs.choose(feat, aln_match_short), {"5' Anchored Overlap"})
        self.assertEqual(fs.choose(feat, aln_match_long), {"5' Anchored Overlap"})
        self.assertEqual(fs.choose(feat, aln_spill_lo_short), set())
        self.assertEqual(fs.choose(feat, aln_spill_lo_long), set())

    """Does FeatureSelector.choose() correctly select features with `3' anchored` interval matching rules?"""

    def test_feature_selector_3p_interval(self):
        iv_rule = "3' anchored"
        chrom, strand, start, stop = "n/a", ".", 5, 10
        feat, fs = self.make_feature_for_interval_test(iv_rule, "3' Anchored Overlap", chrom, strand, start, stop)

        aln_short = {'seq': 'ATGC', 'chrom': chrom, 'strand': strand}
        aln_long = {'seq': 'ATGCNN', 'chrom': chrom, 'strand': strand}
        aln_match_short = make_parsed_sam_record(**dict(aln_short, start=start + 1, name="match, short"))
        aln_match_long = make_parsed_sam_record(**dict(aln_long, start=start - 1, name="match, long"))
        aln_spill_lo_short = make_parsed_sam_record(**dict(aln_short, start=start + 2, name="spill hi, short"))
        aln_spill_lo_long = make_parsed_sam_record(**dict(aln_long, start=start, name="spill hi, long"))

        """
        aln_short:               |ATGC|
        aln_long:              |ATGCNN|
        feat:                 5 |-----| 10
        aln_match_short:       6 |ATGC| 10
        aln_match_long:      4 |ATGCNN| 10
        aln_spill_hi_short:     7 |ATGC| 11
        aln_spill_hi_long:    5 |ATGCNN| 11
        """

        self.assertEqual(fs.choose(feat, aln_match_short), {"3' Anchored Overlap"})
        self.assertEqual(fs.choose(feat, aln_match_long), {"3' Anchored Overlap"})
        self.assertEqual(fs.choose(feat, aln_spill_lo_short), set())
        self.assertEqual(fs.choose(feat, aln_spill_lo_long), set())


if __name__ == '__main__':
    unittest.main()
