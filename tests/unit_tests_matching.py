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
                        ('n/a', 'n/a', "3' anchored"),
                        ('n/a', 'n/a', "5'anchored"),   # spaces should be optional
                        ('n/a', 'n/a', "3'anchored")]

        result = fs.build_interval_selectors(iv, match_tuples)

        self.assertEqual(len(result), 1)
        self.assertIsInstance(result[iv][0][2], IntervalPartialMatch)
        self.assertIsInstance(result[iv][1][2], IntervalFullMatch)
        self.assertIsInstance(result[iv][2][2], IntervalExactMatch)
        self.assertIsInstance(result[iv][3][2], Interval5pMatch)
        self.assertIsInstance(result[iv][4][2], Interval3pMatch)

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

    """Does StrandMatch properly validate the rule definition?"""

    def test_StrandMatch_rule_validation(self):
        good_defs = ["sense", "antisense", "sense  ", "  antisense"]
        bad_defs = ["sens", "anstisense", "anti sense"]

        for defn in good_defs:
            StrandMatch(defn)  # No error

        for defn in bad_defs:
            errmsg = f'Invalid strand selector: "{defn}"'
            with self.assertRaisesRegex(AssertionError, errmsg):
                StrandMatch(defn)

    """Does NtMatch properly validate the rule definition?"""

    def test_NtMatch_rule_validation(self):
        good_defs = ["A", "T", "G", "C", "N",
                    "A,T", "A,   T", "G,G",
                    "A,T,G,C,N"]
        bad_defs = ["A,", ",", "B", "A,B", " "]
        errmsgs = [
            'Invalid nucleotide selector: ""',
            'Invalid nucleotide selector: ""',
            'Invalid nucleotide selector: "B"',
            'Invalid nucleotide selector: "B"',
            'Invalid nucleotide selector: ""',
        ]

        for defn in good_defs:
            NtMatch(defn)

        for defn, errmsg in zip(bad_defs, errmsgs):
            with self.assertRaisesRegex(AssertionError, errmsg):
                NtMatch(defn)

    """Does NumericalMatch properly validate the rule definition?"""

    def test_NumericalMatch_rule_validation(self):
        good_defs = ["10", "10,11", "10-12", "12-14,   16"]
        bad_defs = [",5", "5 5", " "]

        for defn in good_defs:
            NumericalMatch(defn)

        for defn in bad_defs:
            errmsg = f'Invalid length selector: "{defn}"'
            with self.assertRaisesRegex(AssertionError, errmsg):
                NumericalMatch(defn)


if __name__ == '__main__':
    unittest.main()
