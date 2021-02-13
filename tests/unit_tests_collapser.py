import sys
import unittest
import json
import os
from io import StringIO

import aquatx.srna.collapser as collapser

from unittest.mock import patch, MagicMock, call, mock_open


class MyTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        # Change CWD to test folder if test was invoked from project root (ex: by Travis)
        if os.path.basename(os.getcwd()) == 'aquatx-srna':
            os.chdir(f".{os.sep}tests")

        # 40k line fastq file for testing
        self.fastq_file = './testdata/cel_ws279/Lib303_test.fastq'
        self.fasta_file = './testdata/collapser/Lib303.fa'
        self.count_file = './testdata/collapser/Lib303_counts_reference.json'
        with open(self.count_file, 'r') as f:
            self.fastq_count_reference = json.loads(f.read())

        # Mock single-record fastq.
        # Binary strings because seq_counter opens files in binary mode and only decodes the sequence for performance
        self.min_fastq_sample = b'\n'.join([
            b"@NS500697:103:HGVJCBGX5:1:11110:9504:2793 1:N:0:ATCACG",
            b"GTTTTGTTGGGCTTTCGCGAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCT",
            b"+",
            b"AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAE<EEAEEEEEEEEEE"
        ])
        self.min_count_reference = {
            'GTTTTGTTGGGCTTTCGCGAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCT': 1
        }
        self.min_fasta_sample = ">seq_0_x1\n" \
                                "GTTTTGTTGGGCTTTCGCGAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCT"

    """
    Testing seq_counter() with minimum inputs: single-record fastq and zero length files.
    File operations are mocked.
    """
    def test_seq_counter_min(self):
        # Single record fastq test
        with patch('aquatx.srna.collapser.open', mock_open(read_data=self.min_fastq_sample)):
            min_count = collapser.seq_counter("mockFileName")
            self.assertDictEqual(min_count, self.min_count_reference)

        # Zero length file test
        with patch('aquatx.srna.collapser.open', mock_open(read_data='')):
            zero_count = collapser.seq_counter("mockFileName")
            self.assertDictEqual(zero_count, {})

    """
    Testing seq_counter() with test-length library files
    """
    def test_seq_counter_full(self):
        #
        seq_count_dict = collapser.seq_counter(self.fastq_file)
        self.assertDictEqual(seq_count_dict, self.fastq_count_reference)

    """
    Testing seq2fasta()
    """
    @patch('sys.stdout', new_callable=StringIO)
    @patch('aquatx.srna.collapser.os')
    @patch('aquatx.srna.collapser.open', new_callable=mock_open())
    def test_seq2fasta_min(self, mock_open_f, mock_os, mock_stdout):
        # Patch os.path.exists to only return True if filename is "mockFileExists"
        mock_os.path.exists.configure_mock(side_effect=lambda x: x == "mockFileExists")

        test_map = {
            'seqs': [{}, self.min_count_reference],
            'out_file': ["mockFileExists", "mockFileDoesntExist"],
            'thresh': [0,1],
            'low_count_file': [None, "mockFileExists", "mockFileDoesntExist"]
        }

        for seqs in test_map['seqs']:
            for out_file in test_map['out_file']:
                for thresh in test_map['thresh']:
                    for low_count_file in test_map['low_count_file']:
                        print(f"Testing with: seqs: {seqs}, out_file: {out_file}, thresh: {thresh}, low_count_file: {low_count_file}", file=sys.stderr)
                        collapser.seq2fasta(seqs, out_file, thresh, low_count_file)

                        if out_file == "mockFileExists":
                            self.assertEqual("Error: mockFileExists already exists.\n", mock_stdout.getvalue())
                            mock_os.path.exists.assert_called_once_with(out_file)
                            mock_open_f.assert_not_called()
                            continue
                        # out_file does NOT exist from here on
                        elif low_count_file == None:
                            mock_os.path.exists.assert_called_once_with(out_file)
                            mock_open_f.assert_called_once_with(out_file, "w")
                            mock_open_f.return_value.__enter__().write.assert_called_once_with(self.min_fasta_sample)
                        elif low_count_file != None:
                            mock_os.path.exists.assert_called_with(low_count_file)

                        mock_open_f.reset_mock()
                        mock_stdout = patch('sys.stdout', new_callable=StringIO)
                        mock_os.reset_mock()


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
