import unittest
import os

import aquatx.srna.collapser as collapser


class MyTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        # Change CWD to test folder if test was invoked from project root (ex: by Travis)
        if os.path.basename(os.getcwd()) == 'aquatx-srna':
            os.chdir(f".{os.sep}tests")

        self.fastq_file = './testdata/gonad_seq_full_set/KB1_test.fastq'

    def test_seq_counter(self):
        collapser.seq_counter()


if __name__ == '__main__':
    unittest.main()
