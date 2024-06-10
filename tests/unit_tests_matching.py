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
        good_defs = ["10", "10,11", "1-1", "10-12", "1 -2 ", "12 -14,  16"]
        bad_defs = [",5", "5 5", " ", "1-2-3"]

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

    """Does AdarEditMatch work as expected?"""

    def test_AdarEditMatch(self):
        good_alns = [                                                  # Reference Sequence:
            {'Seq': "G",    'Strand': True,  'MD': "0A0",     'NM': 1},  # A (+)
            {'Seq': "C",    'Strand': False, 'MD': "0T0",     'NM': 1},  # A (-)
            {'Seq': "GGG",  'Strand': True,  'MD': "0A0A0A0", 'NM': 3},  # AAA (+)
            {'Seq': "CCC",  'Strand': False, 'MD': "0T0T0T0", 'NM': 3},  # AAA (-)
            {'Seq': "GTGC", 'Strand': True,  'MD': "0A3",     'NM': 1},  # ATGC (+)
            {'Seq': "GCAC", 'Strand': False, 'MD': "3T0",     'NM': 1},  # ATGC (-)
        ]

        bad_alns = [                                                   # Problem:
            {'Seq': "G",    'Strand': True,  'MD': "1",       'NM': 0},  # No mismatch (+)
            {'Seq': "C",    'Strand': False, 'MD': "1",       'NM': 0},  # No mismatch (-)
            {'Seq': "G",    'Strand': True,  'MD': "0T0",     'NM': 1},  # Mismatch, but not from A (+)
            {'Seq': "C",    'Strand': False, 'MD': "0A0",     'NM': 1},  # Mismatch, but not from A (-)
            {'Seq': "C",    'Strand': True,  'MD': "0T0",     'NM': 1},  # A -> G but wrong strand (+)
            {'Seq': "G",    'Strand': False, 'MD': "0A0",     'NM': 1},  # A -> G but wrong strand (-)
            {'Seq': "GTGC", 'Strand': True,  'MD': "4",       'NM': 0},  # No mismatch (+)
            {'Seq': "GCAC", 'Strand': False, 'MD': "4",       'NM': 0},  # No mismatch (-)
            {'Seq': "GTGC", 'Strand': True,  'MD': "0A2T0",   'NM': 2},  # Mismatch from A and other (+)
            {'Seq': "GCAC", 'Strand': False, 'MD': "0A2T0",   'NM': 2},  # Mismatch from A and other (-)
            {'Seq': "G",    'Strand': True,  'MD': "0A0^T0",  'NM': 2},  # Mismatch from A and deletion (+)
            {'Seq': "C",    'Strand': False, 'MD': "0T0^C0",  'NM': 2},  # Mismatch from A and deletion (-)
        ]

        # Test range
        try:
            range_match = AdarEditMatch("1-3")
            for aln in good_alns:
                self.assertTrue(aln in range_match)
            for aln in bad_alns:
                self.assertFalse(aln in range_match)
        except Exception as e:
            print("Failing alignment: \n" + str(aln))
            raise e

        # Test single value
        val_match_lo = AdarEditMatch("1")
        val_match_hi = AdarEditMatch("3")
        aln_lo, aln_hi = good_alns[0], good_alns[2]
        self.assertTrue(aln_lo['NM'] == 1 and aln_hi['NM'] == 3)  # sanity check
        self.assertTrue(aln_lo in val_match_lo and aln_lo not in val_match_hi)
        self.assertTrue(aln_hi in val_match_hi and aln_hi not in val_match_lo)

    """Does TutEditMatch work as expected?"""

    def test_TutEditMatch(self):
        good_alns = [                                                  # Reference Sequence:
            {'Seq': "T",    'Strand': True,  'MD': "0G0",     'NM': 1},  # G (+)
            {'Seq': "A",    'Strand': False, 'MD': "0C0",     'NM': 1},  # C (-)
            {'Seq': "TTT",  'Strand': True,  'MD': "0C0C0A0", 'NM': 3},  # CCA (+)
            {'Seq': "AAA",  'Strand': False, 'MD': "0G0C0T0", 'NM': 3},  # AGC (-)
            {'Seq': "ATTT", 'Strand': True,  'MD': "2C0C0",   'NM': 2},  # ATCC (+)
            {'Seq': "AAGC", 'Strand': False, 'MD': "0G3",     'NM': 1},  # GCTC (-)
        ]

        bad_alns = [                                                   # Problem:
            {'Seq': "T",    'Strand': True,  'MD': "1",       'NM': 0},  # No mismatch (+)
            {'Seq': "A",    'Strand': False, 'MD': "1",       'NM': 0},  # No mismatch (-)
            {'Seq': "G",    'Strand': True,  'MD': "0T0",     'NM': 1},  # Mismatch, but not to U (+)
            {'Seq': "G",    'Strand': False, 'MD': "0A0",     'NM': 1},  # Mismatch, but not to U (-)
            {'Seq': "A",    'Strand': True,  'MD': "0G0",     'NM': 1},  # N -> U but wrong strand (+)
            {'Seq': "T",    'Strand': False, 'MD': "0G0",     'NM': 1},  # N -> U but wrong strand (-)
            {'Seq': "GTGC", 'Strand': True,  'MD': "4",       'NM': 0},  # No mismatch (+)
            {'Seq': "GCAC", 'Strand': False, 'MD': "4",       'NM': 0},  # No mismatch (-)
            {'Seq': "GTTT", 'Strand': True,  'MD': "0T1G0A0", 'NM': 3},  # Mismatch to U and other (+)
            {'Seq': "AAAG", 'Strand': False, 'MD': "0C2T0",   'NM': 2},  # Mismatch to U and other (-)
            {'Seq': "T",    'Strand': True,  'MD': "0C0^T0",  'NM': 2},  # Mismatch to U and deletion (+)
            {'Seq': "A",    'Strand': False, 'MD': "0G0^A0",  'NM': 2},  # Mismatch to U and deletion (-)
        ]

        # Test range
        try:
            range_match = TutEditMatch("1-3")
            for aln in good_alns:
                self.assertTrue(aln in range_match)
            for aln in bad_alns:
                self.assertFalse(aln in range_match)
        except Exception as e:
            print("Failing alignment: \n" + str(aln))
            raise e

        # Test single value
        val_match_lo = TutEditMatch("1")
        val_match_hi = TutEditMatch("3")
        aln_lo, aln_hi = good_alns[0], good_alns[2]
        self.assertTrue(aln_lo['NM'] == 1 and aln_hi['NM'] == 3)
        self.assertTrue(aln_lo in val_match_lo and aln_lo not in val_match_hi)
        self.assertTrue(aln_hi in val_match_hi and aln_hi not in val_match_lo)

        # Test terminal run of desired length, but below mismatch threshold
        tricky = TutEditMatch("3")
        self.assertTrue(good_alns[4]['Seq'] == "ATTT" and good_alns[4]['NM'] == 2)
        self.assertFalse(good_alns[4] in tricky)


if __name__ == '__main__':
    unittest.main()
