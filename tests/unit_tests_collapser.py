import sys
import unittest
import json
import os
from collections import OrderedDict
from io import StringIO

import aquatx.srna.collapser as collapser

from unittest.mock import patch, MagicMock, call, mock_open
from tests.unit_test_helpers import reset_mocks


class MyTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        # Change CWD to test folder if test was invoked from project root (ex: by Travis)
        if os.path.basename(os.getcwd()) == 'aquatx-srna':
            os.chdir(f".{os.sep}tests")

        # Simply for convenience for loading files during setup
        def read(file):
            with open(file, 'r') as f:
                return f.read()

        # Test-length files
        self.fastq_file = './testdata/cel_ws279/Lib303_test.fastq'
        self.fastq_counts_dict = json.loads(read('./testdata/collapser/Lib303_counts_reference.json'))
        self.fasta = {
            "thresh=0": read('testdata/collapser/Lib303_thresh_0.fa'),
            "thresh=4": read('testdata/collapser/Lib303_thresh_4.fa'),
            "thresh=4,low_count": read('testdata/collapser/Lib303_thresh_4_lowcount.fa')
        }

        # File-related messages
        self.prefix_exists_fn = lambda x, y: y in self.output["file"]["exists"].values()
        self.output = {"file": {
            "exists": {"out": "mockPrefixExists_collapsed.fa", "low": "mockPrefixExists_collapsed_lowcounts.fa"},
            "dne": {"out": "mockPrefixDNE_collapsed.fa", "low": "mockPrefixDNE_collapsed_lowcounts.fa"}}}
        self.output["msg"] = {k: "Collapser critical error: "+v+" already exists.\n" for k,v in self.output["file"]["exists"].items()}
        self.prefix_required_msg = "Collapser critical error: an output file must be specified.\n"

        # Min-length fastq/fasta (single record)
        self.min_seq = "GTTTTGTTGGGCTTTCGCGAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCT"
        # Binary strings because seq_counter opens files in binary mode and only decodes the sequence for performance
        self.min_fastq = b'\n'.join([
            b"@NS500697:103:HGVJCBGX5:1:11110:9504:2793 1:N:0:ATCACG",
            self.min_seq.encode('ascii'),
            b"+",
            b"AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAE<EEAEEEEEEEEEE"
        ])
        self.min_counts_dict = {self.min_seq: 1}
        self.min_fasta = ">seq_0_x1\n" + self.min_seq

    """
    Testing seq_counter() with minimum inputs: single-record fastq and zero length files.
    File operations are mocked.
    """
    def test_seq_counter_min(self):
        # Single record fastq test
        with patch('aquatx.srna.collapser.open', mock_open(read_data=self.min_fastq)):
            min_count = collapser.seq_counter("mockFileName")
            self.assertDictEqual(min_count, self.min_counts_dict)

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
        self.assertDictEqual(seq_count_dict, self.fastq_counts_dict)

    """
    Testing that the correct usage messages are produced when improperly calling seq2fasta(),
    or when the specified prefix conflicts with files that already exist.
    """
    @patch('sys.stdout', new_callable=StringIO)
    @patch('aquatx.srna.collapser.os', autospec=True)
    @patch('aquatx.srna.collapser.open', new_callable=mock_open())
    def test_seq2fasta_usage(self, mock_open_f, mock_os, mock_stdout):
        # Simulate that prefix "mockPrefixExists" exists for both output files
        mock_os.path.isfile.configure_mock(side_effect=self.prefix_exists_fn)

        # Test minimum requirements assertions
        with self.assertRaises(AssertionError) as cm:
            collapser.seq2fasta({}, None)
            collapser.seq2fasta({}, "mockPrefixDoesntExist", -1)
        mock_os.path.isfile.assert_not_called()
        mock_open_f.assert_not_called()
        reset_mocks(mock_open_f, mock_os, mock_stdout)

        # Output file exists
        collapser.seq2fasta({}, "mockPrefixExists")
        mock_os.path.isfile.assert_called_once_with(self.output["file"]["exists"]["out"])
        mock_open_f.assert_not_called()
        self.assertEqual(self.output["msg"]["out"], mock_stdout.getvalue())
        reset_mocks(mock_open_f, mock_os, mock_stdout)

        # Low counts file exists
        mock_os.path.isfile.configure_mock(side_effect=lambda x: x == self.output["file"]["exists"]["low"])
        collapser.seq2fasta({}, "mockPrefixExists")
        mock_os.path.isfile.assert_has_calls([
            call(self.output["file"]["exists"]["out"]),
            call(self.output["file"]["exists"]["low"])
        ])
        mock_open_f.assert_not_called()
        self.assertEqual(self.output["msg"]["low"], mock_stdout.getvalue())


    """
    Testing seq2fasta() with all permutations of the "min" parameters defined in test_map.
    The purpose of this test is to see how the script handles minimum conditions: empty
    sequence count dictionary, single record fastq, and thresholds of 0 and 1.
    """
    @patch('sys.stdout', new_callable=StringIO)
    @patch('aquatx.srna.collapser.os', autospec=True)
    @patch('aquatx.srna.collapser.open', new_callable=mock_open())
    def test_seq2fasta_min(self, mock_open_f, mock_os, mock_stdout):
        # Simulate that prefix "mockPrefixExists" exists for both output files
        mock_os.path.isfile.configure_mock(side_effect=self.prefix_exists_fn)

        # The parameter sets to permute
        test_map = {
            'seqs': [{}, self.min_counts_dict],
            'out_prefix': ["mockPrefixDNE"],
            'thresh': [0, 1]
        }

        seqs: dict; out_prefix: str; thresh: int
        for seqs in test_map['seqs']:
            for out_prefix in test_map['out_prefix']:
                for thresh in test_map['thresh']:

                    print(f"Test case: seqs={seqs}, out_prefix={out_prefix}, thresh={thresh}", file=sys.stderr)
                    collapser.seq2fasta(seqs, out_prefix, thresh)

                    # No output messages should have been produced
                    self.assertEqual(mock_stdout.getvalue(), "")

                    if thresh == 0:
                        # Only the outfile should have been opened for writing. No low count file.
                        mock_open_f.assert_called_once_with(self.output["file"]["dne"]["out"], "w")
                        if seqs == {}:
                            # Empty input sequences should result in empty out file
                            mock_open_f.return_value.__enter__().write.assert_called_once_with('')
                        elif seqs == self.min_counts_dict:
                            mock_open_f.return_value.__enter__().write.assert_called_once_with(self.min_fasta)

                    elif thresh == 1:
                        # Both the outfile and low count file should have been written to
                        if seqs == {}:
                            # Empty input sequences should result in empty out and low count file.
                            mock_open_f.return_value.__enter__().write.assert_has_calls([call(''), call('')])
                        elif seqs == self.min_counts_dict:
                            # An empty outfile and a populated low count file should have been written
                            mock_open_f.return_value.__enter__().write.assert_has_calls([
                                call(''), call(self.min_fasta)
                            ])

                    reset_mocks(mock_open_f, mock_os, mock_stdout)


    """
    Testing seq2fasta() with all permutations of the "test-length" parameters defined in test_map.
    Since test_seq2fasta_min() only tested with a single record fastq at most, this test 
    """
    @patch('sys.stdout', new_callable=StringIO)
    @patch('aquatx.srna.collapser.os', autospec=True)
    @patch('aquatx.srna.collapser.open', new_callable=mock_open())
    def test_seq2fasta_full(self, mock_open_f, mock_os, mock_stdout):
        # Simulate that prefix "mockPrefixExists" exists for both output files
        mock_os.path.isfile.configure_mock(side_effect=self.prefix_exists_fn)

        # The parameter sets to permute
        test_map = {
            'seqs': [self.fastq_counts_dict],
            'out_file': ["mockPrefixDNE"],
            'thresh': [0, 4]
        }

        seqs: dict; out_file: str; thresh: int
        for seqs in test_map['seqs']:
            for out_file in test_map['out_file']:
                for thresh in test_map['thresh']:

                    print(f"Test case: seqs=Lib303, out_file={out_file}, thresh={thresh}", file=sys.stderr)
                    collapser.seq2fasta(seqs, out_file, thresh)

                    # No output messages should have been produced
                    self.assertEqual(mock_stdout.getvalue(), "")

                    if thresh == 0:
                        # Only the outfile should have been opened for writing. No low count file.
                        mock_open_f.assert_called_once_with(self.output["file"]["dne"]["out"], "w")
                        mock_open_f.return_value.__enter__().write.assert_called_once_with(self.fasta["thresh=0"])
                    elif thresh == 4:
                        # Both the outfile and low counts file should have been opened for writing. No blank files.
                        mock_open_f.return_value.__enter__().write.assert_has_calls([
                            call(self.fasta["thresh=4"]), call(self.fasta["thresh=4,low_count"])
                        ])

                    reset_mocks(mock_open_f, mock_os, mock_stdout)


    """
    Testing fasta headers for correctness.
    
    The reported index should indicate which UNIQUE sequence it is, where the first unique 
    sequence is index 0, and unique sequence index 1 may have n repetitions of index 0 before 
    it in the fastq file, but its header should report index 1, not index n+1. 
    
    The reported count should indicate how many times the unique sequence occurred in the
    input fastq file. Reported counts are verified for correctness and are verified to be
    in the correct above/below threshold=4 file.
    """
    def test_fasta_headers(self):
        # Build lists of tuples of records for both low and high counts
        # Each record tuple is the pair (fasta sequence header, sequence)
        four_lines = self.fasta["thresh=4"].splitlines()
        lo_lines = self.fasta["thresh=4,low_count"].splitlines()
        f4 = [(header,seq) for header,seq in zip(four_lines[0::2], four_lines[1::2])]
        f4_lo = [(header,seq) for header,seq in zip(lo_lines[0::2], lo_lines[1::2])]

        # Record all unique sequences and their relative position (ctr indicates the nth unique sequence)
        # This will help us verify the index reported in the fasta headers
        seq_to_line, ctr = {}, 0
        with open(self.fastq_file, 'rb') as f:
            while f.readline():
                seq = f.readline()[:-1].decode("utf-8")
                if not seq_to_line.get(seq, None):
                    seq_to_line[seq] = ctr
                    ctr += 1
                f.readline()
                f.readline()

        # Verify the index reported in the high-counts sequence headers
        # Verify the sequence headers indicate counts > threshold=4
        # Verify there are no duplicate sequences
        header_counts_reconstruct = OrderedDict()
        for header,seq in f4:
            seq_index = int(header.split("_")[1])
            seq_count = int(header.split("_")[2][1:])  # Omit the "x" in >seq_INDEX_xCOUNT
            self.assertGreater(int(seq_count), 4, f"Record below threshold was written to out_file: {seq}: {seq_count}")
            self.assertEqual(seq_index, seq_to_line.get(seq, None), f"Unique index was not properly recorded: {seq}: {seq_index}")
            self.assertNotIn(seq, header_counts_reconstruct, f"Duplicate sequence encountered: {seq}")
            header_counts_reconstruct[seq] = seq_count

        # Verify the index reported in the low-counts sequence headers
        # Verify the sequence headers indicate counts <= threshold=4
        # Verify there are no duplicate sequences
        for header,seq in f4_lo:
            seq_index = int(header.split("_")[1])
            seq_count = int(header.split("_")[2][1:])  # Omit the "x" in >seq_INDEX_xCOUNT
            self.assertLessEqual(int(seq_count), 4, f"Record above threshold was written to low_counts_file: {seq}: {seq_count}")
            self.assertEqual(seq_index, seq_to_line.get(seq, None), f"Unique index was not properly recorded: {seq}: {seq_index}")
            self.assertNotIn(seq, header_counts_reconstruct, f"Duplicate sequence encountered: {seq}")
            header_counts_reconstruct[seq] = seq_count

        # Finally, verify all counts
        for seq,count in header_counts_reconstruct.items():
            self.assertEqual(count, self.fastq_counts_dict.get(seq, None),
                             f"Count discrepancy with sequence {seq}")


if __name__ == '__main__':
    unittest.main()
