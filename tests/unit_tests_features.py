import HTSeq
import unittest

from unittest.mock import patch, call, mock_open
from tiny.rna.counter.hts_parsing import read_SAM, parse_alignments
from tiny.rna.counter.features import *
from unit_test_helpers import read, complement

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

    def make_parsed_sam_record(self, name="0_count=1", seq="CAAGACAGAGCTTCACCGTTC", nt5='C',
                        chrom='I', start=15064570, strand='+'):

        if (strand == '+' and nt5 != seq[0]) or (strand == '-' and nt5 != complement[seq[-1]]):
            raise ValueError("Invalid nt5 value for strand")

        return {
            "name": name,
            "len": len(seq),
            "seq": seq,
            "nt5": nt5,
            "chrom": chrom,
            "start": start,
            "end": start + len(seq),
            "strand": strand
        }

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
        sam_aln = self.make_parsed_sam_record(**{'start': 2, 'chrom': chrom, 'strand': strand})

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
        sam_aln = self.make_parsed_sam_record(**{'chrom': chrom, 'strand': strand, 'start': 0, 'seq': 'AT', 'nt5': 'A'})

        """
        iv_none:    2 |-| 3
        iv_olap:  1 |-| 2
        sam_aln: 0 |--| 2
        """

        with patch("tiny.rna.counter.features.FeatureCounter") as mock:
            instance = mock.return_value
            FeatureCounter.assign_features(instance, sam_aln)

        expected_match_list = [(1, 2, {'Overlaps alignment by one base'})]
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


if __name__ == '__main__':
    unittest.main()
