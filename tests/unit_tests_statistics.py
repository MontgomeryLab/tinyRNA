import sys
import unittest
from unittest.mock import patch, mock_open, create_autospec

from tiny.rna.counter.statistics import LibraryStats, SummaryStats


class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, False)

    def create_mock_libstats(self, lib=None):
        library = {
            'Name': 'mock_name',
            'File': 'mock_file',
            'collapsed': 'mock_file',
            'fastp_log': 'mock_file'
        }

        if lib: library.update(lib)
        libstat_obj = LibraryStats()
        libstat_obj.assign_library(library)
        return libstat_obj


    def test_collapser_stat(self):
        # Collapser's FASTA headers are structured as (zero-based) INDEX_count=COUNT.
        # Headers are sorted by index, so the index in the last header is the (count - 1)
        # of unique sequences in the quality-filtered sample

        fasta = './testdata/collapser/Lib303_thresh_0_collapsed.fa'
        mock_libstats = self.create_mock_libstats({'collapsed': fasta})
        expected_count = 9631

        with patch('tiny.rna.counter.statistics.SummaryStats') as mock:
            instance = mock.return_value
            uniq_seq_count = SummaryStats.get_collapser_stats(instance, mock_libstats)

        with open(fasta) as f: sanity_check = len(f.readlines())/2

        self.assertEqual(uniq_seq_count, expected_count)
        self.assertEqual(sanity_check, expected_count)


if __name__ == '__main__':
    unittest.main()
