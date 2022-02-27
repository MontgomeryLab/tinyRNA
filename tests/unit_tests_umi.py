import unittest

from tiny.rna.bitpack.umi import *

class MyTestCase(unittest.TestCase):
    def test_construct(self):
        u = UMI()
        u5 = UMI5p()
        u3 = UMI3p()
        ub = UMIboth()

    def test_factory_construct(self):
        f_5p = UMIFactory(len_5p=1)
        f_3p = UMIFactory(len_3p=1)
        f_bo = UMIFactory(len_5p=1, len_3p=1)

        self.assertIsInstance(f_5p.from_bytes(b"ATGC"), UMI5p)
        self.assertIsInstance(f_3p.from_bytes(b"ATGC"), UMI3p)
        self.assertIsInstance(f_bo.from_bytes(b"ATGC"), UMIboth)

    def test_seq_basic(self):
        UMIFactory(len_5p=5).from_bytes(b"ATGC").printit()


if __name__ == '__main__':
    unittest.main()
