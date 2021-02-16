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

        # Test-length files and file-related messages
        self.fastq_file = './testdata/cel_ws279/Lib303_test.fastq'
        self.fasta_file = 'testdata/collapser/Lib303.fa'
        self.count_file = './testdata/collapser/Lib303_counts_reference.json'
        self.file_exists_fn = lambda x,y: y == "mockFileExists"
        self.file_exists_msg = "Collapser critical error: mockFileExists already exists.\n"
        self.file_required_msg = "Collapser critical error: an output file must be specified.\n"

        with open(self.fasta_file, 'r') as f:
            self.fasta_correct = f.read()
        with open(self.count_file, 'r') as f:
            self.fastq_count_correct = json.loads(f.read())

        # Min-length fastq/fasta (single record)
        self.min_seq = "GTTTTGTTGGGCTTTCGCGAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCT"
        # Binary strings because seq_counter opens files in binary mode and only decodes the sequence for performance
        self.min_fastq_sample = b'\n'.join([
            b"@NS500697:103:HGVJCBGX5:1:11110:9504:2793 1:N:0:ATCACG",
            self.min_seq.encode('ascii'),
            b"+",
            b"AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAE<EEAEEEEEEEEEE"
        ])
        self.min_count_reference = {self.min_seq: 1}
        self.min_fasta_correct = ">seq_0_x1\n" + self.min_seq

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
        self.assertDictEqual(seq_count_dict, self.fastq_count_correct)

    """
    Testing seq2fasta() with all permutations of the "min" parameters defined in test_map
    """
    @patch('sys.stdout', new_callable=StringIO)
    @patch('aquatx.srna.collapser.os', autospec=True)
    @patch('aquatx.srna.collapser.open', new_callable=mock_open())
    def test_seq2fasta_min(self, mock_open_f, mock_os, mock_stdout):
        # Patch os.path.exists to only return True if filename is "mockFileExists"
        mock_os.path.isfile.configure_mock(side_effect=self.file_exists_fn)

        # The parameter sets to permute
        test_map = {
            'seqs': [{}, self.min_count_reference],
            'out_file': [None, "mockFileExists", "mockFileDoesntExist"],
            'thresh': [0, 1],
            'low_count_file': [None, "mockFileExists", "mockFileDoesntExist"]
        }

        def reset_mocks():
            mock_open_f.reset_mock()
            mock_os.reset_mock()
            mock_stdout.truncate(0)
            mock_stdout.seek(0)  # Avoids prepending a null string of previous buffer size

        seqs: dict; out_file: str; thresh: int; low_count_file: str
        for seqs in test_map['seqs']:
            for out_file in test_map['out_file']:
                for thresh in test_map['thresh']:
                    for low_count_file in test_map['low_count_file']:

                        print(f"Testing with: seqs: {seqs}, out_file: {out_file}, thresh: {thresh}, low_count_file: {low_count_file}", file=sys.stderr)
                        collapser.seq2fasta(seqs, out_file, thresh, low_count_file)

                        # Neither above nor below threshold filtering
                        if out_file in ["mockFileExists", None]:
                            mock_open_f.assert_not_called()
                            if out_file == "mockFileExists":
                                self.assertEqual(mock_stdout.getvalue(), self.file_exists_msg)
                                mock_os.path.isfile.assert_called_once_with(out_file)
                            else:
                                self.assertEqual(mock_stdout.getvalue(), self.file_required_msg)
                                mock_os.path.isfile.assert_not_called()

                            reset_mocks()
                            continue
                        # out_file does NOT exist in the following checks due to the above continue statement

                        # Only above threshold filtering
                        elif low_count_file is None:
                            self.assertEqual(mock_stdout.getvalue(), "")
                            mock_os.path.isfile.assert_called_once_with(out_file)
                            mock_open_f.assert_called_once_with(out_file, "w")
                            if seqs == {} or thresh == 1:
                                mock_open_f.return_value.__enter__().write.assert_called_once_with('')
                            else:
                                mock_open_f.return_value.__enter__().write.assert_called_once_with(
                                    self.min_fasta_correct)

                        # Both above and below threshold filtering
                        elif low_count_file is not None:
                            mock_os.path.isfile.assert_has_calls([call(out_file), call(low_count_file)])
                            if low_count_file == "mockFileExists":
                                self.assertEqual(mock_stdout.getvalue(), self.file_exists_msg)
                                mock_open_f.assert_not_called()
                            else:
                                self.assertEqual(mock_stdout.getvalue(), "")
                                mock_open_f.assert_called_with(low_count_file, "w")
                                if seqs == {}:
                                    # Should write nothing for above or below threshold
                                    mock_open_f.return_value.__enter__().write.assert_has_calls([
                                        call(''), call('')
                                    ])
                                elif seqs != {} and thresh == 0:
                                    # Should write something for above but nothing for below threshold
                                    mock_open_f.return_value.__enter__().write.assert_has_calls([
                                        call(self.min_fasta_correct), call('')
                                    ])
                                elif seqs != {} and thresh == 1:
                                    # Should write nothing for above but something for below threshold
                                    mock_open_f.return_value.__enter__().write.assert_has_calls([
                                        call(''), call(self.min_fasta_correct)
                                    ])

                        reset_mocks()


    """
    Testing seq2fasta() with all permutations of the "test-length" parameters defined in test_map
    """
    @patch('sys.stdout', new_callable=StringIO)
    @patch('aquatx.srna.collapser.os', autospec=True)
    @patch('aquatx.srna.collapser.open', new_callable=mock_open())
    def test_seq2fasta_full(self, mock_open_f, mock_os, mock_stdout):
        # Patch os.path.exists to only return True if filename is "mockFileExists"
        mock_os.path.isfile.configure_mock(side_effect=self.file_exists_fn)

        # The parameter sets to permute
        test_map = {
            'seqs': [self.fastq_count_correct],
            'out_file': ["mockFileDoesntExisp"],
            'thresh': [0, 4],
            'low_count_file': [None, "mockFileDoesntExisp"]
        }

        def reset_mocks():
            mock_open_f.reset_mock()
            mock_os.reset_mock()
            mock_stdout.truncate(0)
            mock_stdout.seek(0)  # Avoids prepending a null string of previous buffer size

        seqs: dict; out_file: str; thresh: int; low_count_file: str
        for seqs in test_map['seqs']:
            for out_file in test_map['out_file']:
                for thresh in test_map['thresh']:
                    for low_count_file in test_map['low_count_file']:
                        print(
                            f"Testing with: seqs: Lib303, out_file: {out_file}, thresh: {thresh}, low_count_file: {low_count_file}",
                            file=sys.stderr)
                        collapser.seq2fasta(seqs, out_file, thresh, low_count_file)

                        if thresh == 0 and low_count_file is None:
                            self.assertEqual(mock_stdout.getvalue(), '')
                            mock_open_f.return_value.__enter__().write.assert_called_once_with(self.fasta_correct)
                        if thresh == 0 and low_count_file is not None:
                            self.assertEqual(mock_stdout.getvalue(), '')
                            mock_open_f.return_value.__enter__().write.assert_has_calls([
                                call(self.fasta_correct), call('')
                            ])



        # test correctness of above/below rosters
        # test threshold without defined outfile
        # test threshold < 0
        # test threshold = inf
        # Verify relative order (header index)
        # Verify counts (header count, i.e. _xN)
        # Test correctness with min filesizes (empty and 4 lines (1 record))


if __name__ == '__main__':
    unittest.main()
