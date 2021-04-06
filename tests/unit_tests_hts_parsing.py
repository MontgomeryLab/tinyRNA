import unittest
import HTSeq

from aquatx.srna.hts_parsing import *


class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, False)
        self.short_sam_file = "./testdata/counter/short.sam"

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

if __name__ == '__main__':
    unittest.main()
