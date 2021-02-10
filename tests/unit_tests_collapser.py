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

        self.fastq_file = './testdata/gonad_seq_full_set/KB1_raw.fastq'

    def test_seq_counter(self):
        seq_count_dict = collapser.seq_counter(self.fastq_file)
        collapser.seq2fasta(seq_count_dict, "out.fa")

    def test_verify_counts_file(self):
        ref = {}
        with open("testdata/collapser/KB1_counts_reference.txt", 'r') as f:
            content = f.read().splitlines(keepends=False)
            for line in content:
                rec = line.split(" ")
                if ref.get(rec[0], None) is not None:
                    print("Duplicate sequence: " + rec[0])
                ref[rec[0]] = int(rec[1])

        seq_count_dict = collapser.seq_counter(self.fastq_file)
        self.assertDictEqual(seq_count_dict, ref)

if __name__ == '__main__':
    unittest.main()
