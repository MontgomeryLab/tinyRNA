import unittest
import HTSeq
from copy import deepcopy

from tiny.rna.counter.features import FeatureSelector
from tiny.rna.counter.matching import *
from unit_test_helpers import rules_template


class MatchingTests(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        pass

    """Does FeatureSelector build the proper interval selectors?"""

    def test_feature_selector_interval_build(self):
        fs = FeatureSelector(deepcopy(rules_template))
        iv = HTSeq.GenomicInterval('I', 0, 10, '+')

        # Match tuples formed during GFF parsing
        match_tuples = [('n/a', 'n/a', 'partial', Wildcard()),
                        ('n/a', 'n/a', 'nested', Wildcard()),
                        ('n/a', 'n/a', 'exact', Wildcard()),
                        ('n/a', 'n/a', 'anchored', Wildcard()),
                        ('n/a', 'n/a', "5' anchored", Wildcard()),
                        ('n/a', 'n/a', "3' anchored", Wildcard()),
                        ('n/a', 'n/a', "5'anchored", Wildcard()),   # spaces should be optional
                        ('n/a', 'n/a', "3'anchored", Wildcard())]
        match_tuples += [('n/a', 'n/a', kwd, Wildcard()) for kwd in Wildcard.kwds]

        result = fs.build_interval_selectors(iv, match_tuples)

        self.assertEqual(len(result), 1)
        self.assertIsInstance(result[iv][0][2], IntervalPartialMatch)
        self.assertIsInstance(result[iv][1][2], IntervalNestedMatch)
        self.assertIsInstance(result[iv][2][2], IntervalExactMatch)
        self.assertIsInstance(result[iv][3][2], IntervalAnchorMatch)
        self.assertIsInstance(result[iv][4][2], Interval5pMatch)
        self.assertIsInstance(result[iv][5][2], Interval3pMatch)
        self.assertTrue(isinstance(o, Wildcard) for o in result[iv][6:])

    """Are interval selectors hashable and properly compare for equality?
    This is important for storing feature records in GenomicArraysOfSets"""

    def test_interval_selector_equality(self):
        # Test equality check by instance class name
        # All selectors are assigned the same interval
        shared_iv = HTSeq.GenomicInterval('I', 0, 10)
        class_name_test = {
            IntervalPartialMatch(shared_iv),
            IntervalNestedMatch(shared_iv),
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

    """Does IntervalSelector properly validate the shift parameter?"""

    def test_IntervalSelector_shift_validation(self):
        good_defs = ["1,1", "2, 2", " -3 ,-3", "0,0 ", " -0 , -0 "]
        bad_defs = ["1,", ",2", "3", "4 4", " ", ",", ", "]

        for defn in good_defs:
            IntervalSelector.validate_shift_params(defn)

        for defn in bad_defs:
            errmsg = f'Invalid overlap shift parameters: "{defn}"'
            with self.assertRaisesRegex(AssertionError, errmsg):
                IntervalSelector.validate_shift_params(defn)

    """Does IntervalSelector properly raise IllegalShiftErrors?"""

    def test_IntervalSelector_illegal_shift(self):
        iv = HTSeq.GenomicInterval('n/a', 0, 5, '+')

        with self.assertRaisesRegex(IllegalShiftError, "null interval"):
            IntervalSelector.get_shifted_interval("0,-5", iv)

        with self.assertRaisesRegex(IllegalShiftError, "inverted interval"):
            IntervalSelector.get_shifted_interval("10, 0", iv)

        with self.assertRaisesRegex(IllegalShiftError, "negative start interval"):
            IntervalSelector.get_shifted_interval("-1, 0", iv)


if __name__ == '__main__':
    unittest.main()
