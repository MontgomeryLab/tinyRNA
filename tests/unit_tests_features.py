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
        htsgas = HTSeq.GenomicArrayOfSets("auto", stranded=True)
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
        htsgas = HTSeq.GenomicArrayOfSets("auto", stranded=True)
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
    @patch("tiny.rna.counter.features.parser")
    def test_count_reads_generic(self, parser, mock):
        instance = mock.return_value
        library = {'File': self.short_sam_file, 'Name': 'short'}

        # Note: mock_alignment is wrapped in a double list to imitate multi-alignment bundling
        # It imitates a list of multi-alignment bundles with one bundle having one alignment
        parser.read_SAM.return_value = [["mock_alignment"]]
        bundle = instance.stats.count_bundle.return_value
        instance.assign_features.return_value = ({"mock_feat"}, 1)

        expected_calls_to_stats = [
            call.assign_library(library),
            call.count_bundle(["mock_alignment"]),
            call.count_bundle_alignments(bundle, "mock_alignment", {'mock_feat'}, len({'mock_feat'})),
            call.finalize_bundle(bundle),
            call.diags.write_intermediate_file(library["Name"])
        ]

        # CALL FUNCTION
        FeatureCounter.count_reads(instance, library)

        instance.stats.assign_library.assert_called_once()
        instance.assign_features.assert_called_once_with("mock_alignment")
        instance.stats.assert_has_calls(expected_calls_to_stats)

    """Does FeatureSelector.choose() correctly select features with strict interval matching rules?"""

    def test_feature_selector_strict_interval(self):
        strict = True
        chrom, strand, start, stop = "I", "+", 5, 10
        rules = [dict(rules_template[0], Strict=strict, Identity=('N/A', 'N/A'), nt5end='all', Length='all')]

        # Feature with coordinates 5..10, matching rule 0 with hierarchy 0 and strict interval
        feat = {("Strict Overlap", start, stop, strand, ((0, 0, strict),))}

        aln_base = {'seq': 'ATGC', 'chrom': chrom, 'strand': strand}
        aln_spill_lo = make_parsed_sam_record(**dict(aln_base, start=start - 1, name="spill"))
        aln_spill_hi = make_parsed_sam_record(**dict(aln_base, start=start + 2, name="spill"))
        aln_contained = make_parsed_sam_record(**dict(aln_base, seq="N", start=7, name="contained"))
        aln_contained_lo = make_parsed_sam_record(**dict(aln_base, start=start, name="contained"))
        aln_contained_hi = make_parsed_sam_record(**dict(aln_base, start=start + 1, name="contained"))

        """
        feat:               5 |-----| 10
        aln_spill_lo:      4 |----| 8
        aln_spill_hi:         7 |----| 11
        aln_contained:        7 |-| 8       # Fully contained
        aln_contained_lo:   5 |----| 9      # Shared start position
        aln_contained_hi:    6 |----| 10    # Shared end position
        """

        fs = FeatureSelector(rules, LibraryStats())
        self.assertEqual(fs.choose(feat, aln_spill_lo), set())
        self.assertEqual(fs.choose(feat, aln_spill_hi), set())
        self.assertEqual(fs.choose(feat, aln_contained), {"Strict Overlap"})
        self.assertEqual(fs.choose(feat, aln_contained_lo), {"Strict Overlap"})
        self.assertEqual(fs.choose(feat, aln_contained_hi), {"Strict Overlap"})

    """Does FeatureSelector.choose() correctly select features with partial interval matching rules?"""

    def test_feature_selector_partial_interval(self):
        strict = False
        chrom, strand, start, stop = "I", "+", 5, 10
        rules = [dict(rules_template[0], Strict=strict, Identity=('N/A', 'N/A'), nt5end='all', Length='all')]

        # Feature with coordinates 5..10, matching rule 0 with hierarchy 0 and strict interval
        feat = {("Partial Overlap", start, stop, strand, ((0, 0, strict),))}

        aln_base = {'seq': 'ATGC', 'chrom': chrom, 'strand': strand}
        aln_spill_lo = make_parsed_sam_record(**dict(aln_base, start=start - 1, name="spill"))
        aln_spill_hi = make_parsed_sam_record(**dict(aln_base, start=start + 2, name="spill"))
        aln_contained_lo = make_parsed_sam_record(**dict(aln_base, start=start, name="contained"))
        aln_contained_hi = make_parsed_sam_record(**dict(aln_base, start=start + 1, name="contained"))

        """
        feat:               5 |-----| 10
        aln_spill_lo:      4 |----| 8
        aln_spill_hi:         7 |----| 11
        aln_contained_lo:   5 |----| 9      Shared start position
        aln_contained_hi:    6 |----| 10    Shared end position
        """

        fs = FeatureSelector(rules, LibraryStats())
        self.assertEqual(fs.choose(feat, aln_spill_lo), {"Partial Overlap"})
        self.assertEqual(fs.choose(feat, aln_spill_hi), {"Partial Overlap"})
        self.assertEqual(fs.choose(feat, aln_contained_lo), {"Partial Overlap"})
        self.assertEqual(fs.choose(feat, aln_contained_hi), {"Partial Overlap"})

if __name__ == '__main__':
    unittest.main()
