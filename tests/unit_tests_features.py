import HTSeq
import unittest
import tests.unit_test_helpers as helpers

from unittest.mock import patch, call
from tiny.rna.counter.hts_parsing import Alignment, read_SAM
from tiny.rna.counter.features import *

resources = "./testdata/counter"

class FeaturesTests(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.gff_file = f"{resources}/identity_choice_test.gff3"
        self.short_gff_file = f"{resources}/single.gff3"
        self.short_gff = helpers.read(self.short_gff_file)

        self.sam_file = f"{resources}/identity_choice_test.sam"
        self.short_sam_file = f"{resources}/single.sam"
        self.short_sam = helpers.read(self.short_sam_file)

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
        iv_feat = HTSeq.GenomicInterval("I", 0, 2, "+")
        iv_none = HTSeq.GenomicInterval("I", 2, 3, "+")

        htsgas[iv_feat] += "Should not match"
        Features.chrom_vectors = htsgas.chrom_vectors
        none_alignment = Alignment(iv_none, "Non-overlap", b"A")

        with patch("tiny.rna.counter.features.FeatureCounter") as mock:
            instance = mock.return_value
            FeatureCounter.assign_features(instance, none_alignment)

        instance.choose.assert_not_called()
        instance.stats.chrom_misses.assert_not_called()

    """Does assign_features return features that overlap the query interval by a single base?"""

    def test_assign_features_single_base_overlap(self):
        htsgas = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        iv_feat = HTSeq.GenomicInterval("I", 0, 2, "+")  # The "feature"
        iv_olap = HTSeq.GenomicInterval("I", 1, 2, "+")  # A single-base overlapping feature
        iv_none = HTSeq.GenomicInterval("I", 2, 3, "+")  # A non-overlapping interval

        htsgas[iv_feat] += 'The "feature"'
        htsgas[iv_none] += "Non-overlapping feature"

        Features.chrom_vectors = htsgas.chrom_vectors
        olap_alignment = Alignment(iv_olap, "Single base overlap", b"A")

        with patch("tiny.rna.counter.features.FeatureCounter") as mock:
            instance = mock.return_value
            FeatureCounter.assign_features(instance, olap_alignment)

        expected_match_list = [(1, 2, {'The "feature"'})]
        instance.selector.choose.assert_called_once_with(expected_match_list, olap_alignment)
        instance.stats.chrom_misses.assert_not_called()

    """Does count_reads call the right functions when handling a single record library?"""

    @patch("tiny.rna.counter.features.FeatureCounter")
    @patch("tiny.rna.counter.features.parser")
    def test_count_reads_generic(self, parser, mock):
        instance = mock.return_value
        library = {'File': self.short_sam_file, 'Name': 'short'}
        alignment = next(read_SAM(library['File']))

        parser.read_SAM.return_value = iter([alignment])  # Need to feed HTSeq's bundler an actual alignment
        bundle = instance.stats.count_bundle.return_value
        instance.assign_features.return_value = ({"mock_feat"}, 1)

        expected_calls_to_stats = [
            call.assign_library(library),
            call.count_bundle([alignment]),
            call.count_bundle_alignments(bundle, alignment, {'mock_feat'}, len({'mock_feat'})),
            call.finalize_bundle(bundle),
            call.diags.write_intermediate_file(library["Name"])
        ]

        # CALL FUNCTION
        FeatureCounter.count_reads(instance, library)

        instance.stats.assign_library.assert_called_once()
        instance.assign_features.assert_called_once_with(alignment)
        instance.stats.assert_has_calls(expected_calls_to_stats)


if __name__ == '__main__':
    unittest.main()
