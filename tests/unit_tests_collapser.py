import unittest
import json
import os

import aquatx.srna.collapser as collapser


class MyTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        # Change CWD to test folder if test was invoked from project root (ex: by Travis)
        if os.path.basename(os.getcwd()) == 'aquatx-srna':
            os.chdir(f".{os.sep}tests")

        self.fastq_file = './testdata/cel_ws279/Lib303_test.fastq'

    def test_seq_counter(self):
        with open("testdata/collapser/Lib303_counts_reference.json", 'r') as f:
            reference = json.loads(f.read())

        seq_count_dict = collapser.seq_counter(self.fastq_file)
        self.assertDictEqual(seq_count_dict, reference)

        # Test correctness with min filesizes

    def test_seq2fasta_nolowcount(self):
        seq_count_dict = collapser.seq_counter(self.fastq_file)
        collapser.seq2fasta(seq_count_dict, "Lib303.fa")

        # Verify relative order (header index)
        # Verify counts (header count, i.e. _xN)
        # Test correctness with min filesizes (empty and 4 lines (1 record))

    def test_seq2fasta_lowcount(self):
        collapser.seq2fasta()

        # test correctness of above/below rosters
        # test threshold without defined outfile
        # test threshold < 0
        # test threshold = inf


if __name__ == '__main__':
    unittest.main()
