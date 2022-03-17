import unittest
import HTSeq
from copy import deepcopy

from tiny.rna.counter.features import FeatureSelector
from tiny.rna.counter.matching import *
from tiny.rna.counter.statistics import LibraryStats
from unit_test_helpers import rules_template


class MyTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        pass

    """Does FeatureSelector build the proper interval selectors?"""

    def test_feature_selector_interval_build(self):
        fs = FeatureSelector(deepcopy(rules_template), LibraryStats())
        iv = HTSeq.GenomicInterval('I', 0, 10, '+')

        # Match tuples formed during GFF parsing
        match_tuples = [('n/a', 'n/a', 'partial'),
                        ('n/a', 'n/a', 'full'),
                        ('n/a', 'n/a', 'exact'),
                        ('n/a', 'n/a', "5' anchored"),
                        ('n/a', 'n/a', "3' anchored")]

        result = fs.build_interval_selectors(iv, match_tuples)

        self.assertIsInstance(result[0][2], IntervalPartialMatch)
        self.assertIsInstance(result[1][2], IntervalFullMatch)
        self.assertIsInstance(result[2][2], IntervalExactMatch)
        self.assertIsInstance(result[3][2], Interval5pMatch)
        self.assertIsInstance(result[4][2], Interval3pMatch)

    """Are interval selectors hashable and properly compare for equality?
    This is important for storing feature records in GenomicArraysOfSets"""

    def test_interval_selector_equality(self):
        # Test equality check by instance class name
        # All selectors are assigned the same interval
        shared_iv = HTSeq.GenomicInterval('I', 0, 10)
        class_name_test = {
            IntervalPartialMatch(shared_iv),
            IntervalFullMatch(shared_iv),
            IntervalExactMatch(shared_iv),
            Interval5pMatch(shared_iv),
            Interval3pMatch(shared_iv)
        }

        # Test equality check by start/end coordinates
        ivma = IntervalSelector(HTSeq.GenomicInterval('I', 0, 10))
        ivmb = IntervalSelector(HTSeq.GenomicInterval('I', 1, 11))
        ivmc = IntervalSelector(HTSeq.GenomicInterval('I', 0, 10))

        self.assertEqual(len(class_name_test), 5)
        self.assertEqual(ivma, ivmc)
        self.assertNotEqual(ivma, ivmb)



if __name__ == '__main__':
    unittest.main()
