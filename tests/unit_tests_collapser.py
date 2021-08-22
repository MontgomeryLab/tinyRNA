import unittest
import json
import sys
import os

from tests.unit_test_helpers import read, reset_mocks, ShellCapture, reassemble_gz_w
from unittest.mock import patch, MagicMock, call, mock_open, Mock
from collections import OrderedDict
from io import StringIO

import tiny.srna.collapser as collapser

class MyTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        # Change CWD to test folder if test was invoked from project root (ex: by Travis)
        if os.path.basename(os.getcwd()) == 'tiny':
            os.chdir(f".{os.sep}tests")

        # Test-length files
        self.fastq_file = '../START_HERE/sample_data/Lib303_test.fastq'
        self.fastq_gzip = 'testdata/collapser/Lib303_test.fastq.gz'
        self.fastq_counts_dict = json.loads(read('./testdata/collapser/Lib303_counts_reference.json'))
        self.fasta = {
            "thresh=0": read('testdata/collapser/Lib303_thresh_0_collapsed.fa'),
            "thresh=4": read('testdata/collapser/Lib303_thresh_4_collapsed.fa'),
            "thresh=4,low_count": read('testdata/collapser/Lib303_thresh_4_collapsed_lowcounts.fa'),
            "thresh=0,gz": read('testdata/collapser/Lib303_thresh_0_collapsed.fa.gz', mode='rb'),
            "min,gz": read('testdata/collapser/min_collapsed.fa.gz', mode='rb')
        }

        # File-related messages
        self.prefix_exists_fn = lambda x, y: y in [self.output['file'][f]["exists"] for f in ["out","low"]]
        self.output = {"file": {
            "out": {"exists": "mockPrefixExists_collapsed.fa", "dne": "mockPrefixDNE_collapsed.fa"},
            "low": {"exists": "mockPrefixExists_collapsed_lowcounts.fa", "dne": "mockPrefixDNE_collapsed_lowcounts.fa"}}}
        self.output["msg"] = {k: "Collapser critical error: "+v['exists']+" already exists.\n"
                              for k,v in self.output['file'].items()}
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
        # Binary encoded gzip compressed min_fastq (compression level 6, mtime 0)
        self.min_fastq_gz = b'\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xff\x9d\x8b1\x0e\xc2@\x0c\x04\xfb{\x05=\x8d' \
                            b'\x8fpDgQ`V\xd1"\x8a4X\x88\xff\xbf\x84\xbd \xe5\x01L\xb1\xb6\xe5\xd9\xdb\xfajf\x97>{' \
                            b'\xb5\xc9\x1f|?q\xe7\xa7y\xf5*\xcc{\xb3\xb3\x9f\xe6>\x1d\xaa\xafn\x1e\x89\x00\x0bS(HB' \
                            b'\x1b\x08FPOn\x93\x92\xa4%R\'\x12\x08\x8e\xde\xaf\xac\x94\x97\x91\xc4\xa6\x94c\x89\xc1' \
                            b'\xf2\x1f\xb1\\G\xec|\x01\'\xbdlc\xd2\x00\x00\x00'
        self.min_counts_dict = {self.min_seq: 1}
        self.min_fasta = ">0_count=1\n" + self.min_seq

    """
    Testing seq_counter() with minimum inputs: single-record fastq and zero length files.
    File operations are mocked. Patching looks a little different in this one since
    seq_counter() now has file_reader in its default args, so we have to instead patch
    its __defaults__ attribute.
    """
    def test_seq_counter_min(self):
        # Single record fastq test
        with patch.object(collapser.seq_counter, '__defaults__', new=(mock_open(read_data=self.min_fastq),)) as mo:
            min_count = collapser.seq_counter("mockPrefixDNE")
            mo[0].assert_called_once_with("mockPrefixDNE", "rb")
            self.assertDictEqual(min_count, self.min_counts_dict)
            self.assertEqual([call('mockPrefixDNE', 'rb')], mo[0].call_args_list)
        print("seq_counter: passed single record test.", file=sys.stderr)

        # Zero length file test
        with patch.object(collapser.seq_counter, '__defaults__', new=(mock_open(read_data=''),)) as mo:
            zero_count = collapser.seq_counter("mockPrefixDNE")
            mo[0].assert_called_once_with("mockPrefixDNE", "rb")
            self.assertDictEqual(zero_count, {})
            self.assertEqual([call('mockPrefixDNE', 'rb')], mo[0].call_args_list)
        print("seq_counter: passed zero length input test.", file=sys.stderr)

    """
    Testing seq_counter() with test-length library files
    """
    def test_seq_counter_full(self):
        seq_count_dict = collapser.seq_counter(self.fastq_file)
        self.assertDictEqual(seq_count_dict, self.fastq_counts_dict)
        print("seq_counter: counts verified.", file=sys.stderr)

    """
    Testing gzip reading in seq_counter() 
    """
    def test_seq_counter_gzip(self):
        # MIN TEST
        # Need to patch 1) builtins.open in gzip, 2) file_reader in seq_counter's __defaults__
        with patch('tiny.srna.collapser.gzip.builtins.open', new=mock_open(read_data=self.min_fastq_gz)):
            with patch.object(collapser.seq_counter, '__defaults__', new=(mock_open(read_data=self.min_fastq_gz),)) as mo:
                # Read the mock gzipped single record fastq file
                gz_min_result = collapser.seq_counter("mockPrefixDNE")
                self.assertEqual(self.min_counts_dict, gz_min_result)

        # FULL TEST
        # Verify that full-length fastq.gz files can be read and properly counted
        gz_full_result = collapser.seq_counter(self.fastq_gzip)
        self.assertDictEqual(self.fastq_counts_dict, gz_full_result)

    """
    Testing gzip writing in seq2fasta()
    """
    @patch('tiny.srna.collapser.open', new_callable=mock_open())
    def test_seq2fasta_gzip(self, mock_open_f):
        with patch('tiny.srna.collapser.gzip.builtins.open', new_callable=mock_open) as gz_open:
            # MIN TEST
            collapser.seq2fasta(self.min_counts_dict, "min_gz", gz=True)
            output = reassemble_gz_w(gz_open.mock_calls)
            self.assertEqual(self.fasta['min,gz'], output)
            # Only the gzip.GzipFile() interface should have been used, not builtins.open()
            mock_open_f.assert_not_called()
            # Only a binary write to the output file should have been called
            self.assertEqual([call('min_gz_collapsed.fa.gz', 'wb')], gz_open.call_args_list)
            gz_open.reset_mock()

            # FULL TEST
            collapser.seq2fasta(self.fastq_counts_dict, "Lib303_thresh_0", gz=True)
            output = reassemble_gz_w(gz_open.mock_calls)
            self.assertEqual(self.fasta["thresh=0,gz"], output)
            # Only the gzip.GzipFile() interface should have been used, not builtins.open()
            mock_open_f.assert_not_called()
            # Only a binary write to the output file should have been called
            self.assertEqual([call('Lib303_thresh_0_collapsed.fa.gz', 'wb')], gz_open.call_args_list)

    """
    Testing that the correct usage messages are produced when improperly calling seq2fasta(),
    or when the specified prefix conflicts with files that already exist.
    """
    @patch('sys.stdout', new_callable=StringIO)
    @patch('tiny.srna.collapser.os', autospec=True)
    @patch('tiny.srna.collapser.open', new_callable=mock_open())
    def test_seq2fasta_usage(self, mock_open_f, mock_os, mock_stdout):
        # Simulate that prefix "mockPrefixExists" exists for both output files
        mock_os.path.isfile.configure_mock(side_effect=self.prefix_exists_fn)

        # Test minimum requirements assertions
        with self.assertRaises(AssertionError) as cm:
            collapser.seq2fasta({}, None)
            collapser.seq2fasta({}, "mockPrefixDNE", -1)
        mock_os.path.isfile.assert_not_called()
        mock_open_f.assert_not_called()
        print("seq2fasta: minimum requirements checked.", file=sys.stderr)
        reset_mocks(mock_open_f, mock_os, mock_stdout)

        # Output file exists
        with self.assertRaises(FileExistsError) as cm:
            collapser.seq2fasta({}, "mockPrefixExists")
        mock_os.path.isfile.assert_called_once_with(self.output["file"]["out"]["exists"])
        mock_open_f.assert_not_called()
        reset_mocks(mock_open_f, mock_os, mock_stdout)

        # Low-counts file exists
        with self.assertRaises(FileExistsError) as cm:
            mock_os.path.isfile.configure_mock(side_effect=lambda x: x == self.output["file"]["low"]["exists"])
            collapser.seq2fasta({}, "mockPrefixExists")
        mock_os.path.isfile.assert_has_calls([
            call(self.output["file"]["out"]["exists"]),  # First the output file was checked, and doesn't exist...
            call(self.output["file"]["low"]["exists"])   # but the low count file does
        ])
        mock_open_f.assert_not_called()
        print("seq2fasta: output namespace collision avoidance checked.", file=sys.stderr)

    """
    Testing seq2fasta() with all permutations of the "min" parameters defined in test_map.
    The purpose of this test is to see how the script handles minimum conditions: empty
    sequence count dictionary, single record fastq, and thresholds of 0 and 1.
    """
    @patch('sys.stdout', new_callable=StringIO)
    @patch('tiny.srna.collapser.os', autospec=True)
    @patch('tiny.srna.collapser.open', new_callable=mock_open())
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

                    print(f"seq2fasta: seqs={seqs}, out_prefix={out_prefix}, thresh={thresh}", file=sys.stderr)
                    collapser.seq2fasta(seqs, out_prefix, thresh)

                    # No output messages should have been produced
                    self.assertEqual(mock_stdout.getvalue(), "")

                    if thresh == 0:
                        # Only the outfile should have been opened for writing. No low-count file.
                        mock_open_f.assert_called_once_with(self.output["file"]["out"]["dne"], "w")
                        if seqs == {}:
                            # Empty input sequences should result in empty out file
                            mock_open_f.return_value.__enter__().write.assert_called_once_with('')
                        elif seqs == self.min_counts_dict:
                            mock_open_f.return_value.__enter__().write.assert_called_once_with(self.min_fasta)

                    elif thresh == 1:
                        # Both the outfile and low-count file should have been written to
                        if seqs == {}:
                            # Empty input sequences should result in empty out and low-count file.
                            mock_open_f.return_value.__enter__().write.assert_has_calls([call(''), call('')])
                        elif seqs == self.min_counts_dict:
                            # An empty outfile and a populated low-count file should have been written
                            mock_open_f.return_value.__enter__().write.assert_has_calls([
                                call(''), call(self.min_fasta)
                            ])

                    reset_mocks(mock_open_f, mock_os, mock_stdout)


    """
    Testing seq2fasta() with all permutations of the "test-length" parameters defined in test_map.
    Since test_seq2fasta_min() only tested with a single record fastq at most, this test 
    """
    @patch('sys.stdout', new_callable=StringIO)
    @patch('tiny.srna.collapser.os', autospec=True)
    @patch('tiny.srna.collapser.open', new_callable=mock_open())
    def test_seq2fasta_thresh_0_4(self, mock_open_f, mock_os, mock_stdout):
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

                    print(f"seq2fasta: seqs=Lib303, out_file={out_file}, thresh={thresh}", file=sys.stderr)
                    collapser.seq2fasta(seqs, out_file, thresh)

                    # No output messages should have been produced
                    self.assertEqual(mock_stdout.getvalue(), "")

                    if thresh == 0:
                        # Only the outfile should have been opened for writing. No low-count file.
                        mock_open_f.assert_called_once_with(self.output["file"]["out"]["dne"], "w")
                        mock_open_f.return_value.__enter__().write.assert_called_once_with(self.fasta["thresh=0"])
                    elif thresh == 4:
                        # Both the outfile and low-counts file should have been opened for writing. No blank files.
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
        f4 = [record for record in zip(four_lines[0::2], four_lines[1::2])]
        f4_lo = [record for record in zip(lo_lines[0::2], lo_lines[1::2])]

        # Record all unique sequences and their relative position (first_seen_index indicates the nth unique sequence)
        # This will help us verify the index reported in the fasta headers
        seq_to_index, first_seen_index = {}, 0
        with open(self.fastq_file, 'rb') as f:
            while f.readline():
                seq = f.readline()[:-1].decode("utf-8")
                # If this is the first time this sequence is encountered, record its UNIQUE record number
                if not seq in seq_to_index:
                    seq_to_index[seq] = first_seen_index
                    # Only increment first_seen_index on unique sequences
                    first_seen_index += 1
                f.readline()
                f.readline()

        # Verify the sequence headers in f4 indicate counts > threshold=4
        # Verify the ID/index reported in the high-counts sequence headers
        # Verify there are no duplicate sequences
        header_counts_reconstruct = OrderedDict()
        for header,seq in f4:
            header_seq_index = int(header.split("_")[0][1:])  # Omit the ">" in >ID_count=COUNT
            header_seq_count = int(header.split("=")[1])
            self.assertGreater(header_seq_count, 4, f"Record below threshold was written to out_file: {seq}: {header_seq_count}")
            self.assertEqual(header_seq_index, seq_to_index.get(seq, None), f"Unique ID was not properly recorded: {seq}: {header_seq_index}")
            self.assertNotIn(seq, header_counts_reconstruct, f"Duplicate sequence encountered: {seq}")
            header_counts_reconstruct[seq] = header_seq_count
        print("seq2fasta: fasta headers and content for high counts verified.", file=sys.stderr)

        # Verify the sequence headers in f4_lo indicate counts <= threshold=4
        # Verify the ID/index reported in the low-counts sequence headers
        # Verify there are no duplicate sequences
        for header,seq in f4_lo:
            header_seq_index = int(header.split("_")[0][1:])  # Omit the ">" in >ID_count=COUNT
            header_seq_count = int(header.split("=")[1])
            self.assertLessEqual(header_seq_count, 4, f"Record above threshold was written to low_counts_file: {seq}: {header_seq_count}")
            self.assertEqual(header_seq_index, seq_to_index.get(seq, None), f"Unique ID was not properly recorded: {seq}: {header_seq_index}")
            self.assertNotIn(seq, header_counts_reconstruct, f"Duplicate sequence encountered: {seq}")
            header_counts_reconstruct[seq] = header_seq_count
        print("seq2fasta: fasta headers and content for low counts verified.", file=sys.stderr)

        # Finally, verify all counts by sequence, as reported in the fasta header
        for seq,count in header_counts_reconstruct.items():
            self.assertEqual(count, self.fastq_counts_dict.get(seq, None),
                             f"Count discrepancy with sequence {seq}")
        print("seq2fasta: all counts reported in headers verified.", file=sys.stderr)

    """
    Testing basic command line usage. We're not testing functionality of the script here.
    We just want to know that the installed script is correctly finding input files and
    that execution is halted before writing outputs if outputs already exist.
    """
    def test_collapser_command(self):
        prefix = 'test'
        expected_out_file = prefix + '_collapsed.fa'
        expected_low_file = prefix + '_collapsed_lowcounts.fa'

        # Standard usage test
        with ShellCapture(f'tiny-collapse -i {self.fastq_file} -o {prefix} -t 4') as test:
            test()
            # No expected console output for a non-problematic run
            self.assertEqual(test.get_stdout(), '')
            self.assertEqual(test.get_stderr(), '')
            self.assertTrue(os.path.isfile(expected_out_file))
            self.assertTrue(os.path.isfile(expected_low_file))

        # Namespace collision test
        try:
            test_collapsed_fa_size = os.path.getsize(expected_out_file)
            with ShellCapture(f'tiny-collapse -i /dev/null -o {prefix}') as test:
                test()
                self.assertEqual(test.get_stdout(), '')
                self.assertIn(f"Collapser critical error: {expected_out_file} already exists.\n", test.get_stderr())
                # (Very) roughly tests that the output file of the last test (same prefix) was not modified by this call
                self.assertEqual(test_collapsed_fa_size, os.path.getsize(expected_out_file))
        finally:
            os.remove(expected_out_file)
            os.remove(expected_low_file)

    """
    Testing argparse requirements.
    """
    @patch('tiny.srna.collapser.os', autospec=True)
    @patch('tiny.srna.collapser.gzip.os', autospec=True)
    @patch('sys.stdout', new_callable=StringIO)
    @patch('sys.stderr',  new_callable=StringIO)
    def test_collapser_args(self, mock_stderr, mock_stdout, os_gz, os_aq):
        def collapser_main():
            try:
                collapser.main()
            except SystemExit:
                pass

        def reset_stderr():
            mock_stderr.truncate(0)
            mock_stderr.seek(0)  # Reset stderr capture

        # Negative threshold test
        negative_threshold_args = f"tiny-collapse -i /dev/null -o test -t -1".split(" ")
        with patch('sys.argv', negative_threshold_args) as cm:
            collapser_main()
            self.assertIn("Threshold must be >= 0", mock_stderr.getvalue())
            reset_stderr()

        # Omit prefix test
        no_prefix_args = f"tiny-collapse -i /dev/null -t 1".split(" ")
        with patch('sys.argv', no_prefix_args):
            collapser_main()
            self.assertIn("the following arguments are required: -o/--out-prefix", mock_stderr.getvalue())
            reset_stderr()

        # Omit input test
        no_input_args = f"tiny-collapse -o N/A -t 1".split(" ")
        with patch('sys.argv', no_input_args):
            collapser_main()
            self.assertIn("error: the following arguments are required: -i/--input-file", mock_stderr.getvalue())
            reset_stderr()

        # Ensure helpstring matches the expected
        with patch('sys.argv', ["tiny-collapse", "-h"]):
            collapser_main()
            with open('./testdata/collapser/helpstring.txt', 'r') as f:
                expected_helpstring = f.read()
            # Helpstring is written to stdout, not stderr
            self.assertEqual(expected_helpstring, mock_stdout.getvalue())
            mock_stdout.truncate(0)
            mock_stdout.seek(0)

        # Compression test
        # Lots of patching since we need to mock the file interface for seq_counter default argument,
        # gzip, and seq2fasta. Also need to mock file existence.
        no_compression_args = f"tiny-collapse -i min_gz -o min_gz".split(" ")
        [os.path.isfile.configure_mock(side_effect=self.prefix_exists_fn) for os in [os_aq, os_gz]]
        with patch('tiny.srna.collapser.gzip.builtins.open', new=mock_open(read_data=self.min_fastq_gz)) as gzopen:
            with patch.object(collapser.seq_counter, '__defaults__',
                              new=(mock_open(read_data=self.min_fastq_gz),)) as seq_counter_open:
                with patch('tiny.srna.collapser.open', new_callable=mock_open) as seq2fasta_open:

                    # Without compression (no -c flag)
                    with patch('sys.argv', no_compression_args):
                        collapser_main()
                        # Both interfaces should read at seq_counter, non-gzip write in seq2fasta
                        self.assertEqual([call('min_gz', 'rb')], seq_counter_open[0].call_args_list)
                        self.assertEqual([call('min_gz', 'rb')], gzopen.call_args_list)
                        self.assertEqual([call('min_gz_collapsed.fa', 'w')], seq2fasta_open.call_args_list)
                        gzopen.reset_mock(), seq_counter_open[0].reset_mock(), seq2fasta_open.reset_mock()

                    # With compression (with -c flag)
                    compression_args = f"tiny-collapse -i min_gz -o min_gz -c".split(" ")
                    with patch('sys.argv', compression_args):
                        collapser_main()
                        # Both interfaces should read at seq_counter, but only gzip should write
                        self.assertEqual([call('min_gz', 'rb')], seq_counter_open[0].call_args_list)
                        self.assertEqual([call('min_gz', 'rb'), call('min_gz_collapsed.fa.gz', 'wb')],
                                         gzopen.call_args_list)
                        self.assertEqual([], seq2fasta_open.call_args_list)


if __name__ == '__main__':
    unittest.main()
