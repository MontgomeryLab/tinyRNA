import unittest

from aquatx.srna.hts_parsing import *


class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, False)
        self.short_sam_file = "testdata/counter/single.sam"

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

    """Does our custom SAM parser produce the same pertinent info as HTSeq's BAM_reader?
    
    A note on SAM files: reads are always stored 5' to 3', so antisense reads are actually
    recorded in reverse complement. HTSeq automatically performs this conversion, but we
    are only really concerned about a sequence's 5' end NT, so we perform this conversion
    more surgically via get_nt_5end() for efficiency.
    """

    def test_sam_parser_comparison(self):
        file = "./run_directory/Lib304_test_aligned_seqs.sam"
        ours = read_SAM(file)
        theirs = HTSeq.BAM_Reader(file)

        for our, their in zip(ours, theirs):
            self.assertEqual(our.iv, their.iv)
            self.assertEqual(our.iv.strand, their.iv.strand)
            self.assertEqual(get_nt_5end(our), chr(their.read.seq[0]))  # See note above
            self.assertEqual(our.read.name, their.read.name)
            self.assertEqual(len(our.read), len(their.read))

if __name__ == '__main__':
    unittest.main()
