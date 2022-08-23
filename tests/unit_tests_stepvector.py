import unittest
import HTSeq

from tiny.rna.counter.stepvector._stepvector import StepVector


class StepVectorTests(unittest.TestCase):

    """Does our Cython StepVector accept and return steps as expected?"""

    def test_genomicarray_with_cython_stepvec(self):
        # Patch the StepVector reference in the HTSeq module and use a GenomicArray
        # instead of a GenomicArrayOfSets, just as we do in ReferenceTables
        setattr(HTSeq.StepVector, 'StepVector', StepVector)
        gas = HTSeq.GenomicArray('auto', stranded=False)

        iv1 = HTSeq.GenomicInterval("chr", 0, 5)
        iv2 = HTSeq.GenomicInterval("chr", 5, 7)
        iv3 = HTSeq.GenomicInterval("chr", 4, 6)
        iv4 = HTSeq.GenomicInterval("chr", 7, 9)
        ivs = [iv1, iv2, iv3, iv4]

        f1 = {'featA'}
        f2 = {'featB', 'featC'}
        f3 = {'featD'}
        f4 = {'featE', 'featF', 'featG'}
        fs = [f1, f2, f3, f4]

        """
        iv1:  0 |-----| 5       featA
        iv2:       5 |--| 7     featB, featC
        iv3:      4 |--| 6      featD
        iv4:         7 |--| 9   featE, featF, featG
        """

        for iv, feat in zip(ivs, fs):
            gas[iv] += feat

        exp_overlaps = [
            [(0, 4, f1), (4, 5, f1 | f3)],  # f1 overlaps f3
            [(5, 6, f2 | f3), (6, 7, f2)],  # f2 overlaps f3
            [(4, 5, f1 | f3), (5, 6, f3 | f2)],  # f3 overlaps f1 and f2
            [(7, 9, f4)],  # f4 does not overlap anything else
        ]

        for iv, expected in zip(ivs, exp_overlaps):
            actual = list(gas['chr']['.'].array[iv.start:iv.end].get_steps())
            self.assertEqual(actual, expected)

    def test_compare_to_htseq_steps(self):
        """Test currently doesn't work (and it won't be easy to get it to work)

        This is because HTSeq's StepVector with typecode 'O' handles .add_value() using the
        apply() function and a local function which copies the set OUT of the StepVector,
        adds to it, then puts it back. This is exactly the kind of thing that our
        specialized StepVector was intended to avoid.
        """

        from HTSeq.StepVector import StepVector as hStepVector
        def htseq_add(x):
            y = x.copy()
            y.add(value)
            return y

        sv_len = 200
        ours = StepVector.create(length=sv_len)
        theirs = hStepVector.create(length=sv_len, typecode='O')
        
        iv1 = slice(0, 99)
        iv2 = slice(99, 101)
        iv3 = slice(98, 100)
        iv4 = slice(101, 150)
        ivs = [iv1, iv2, iv3, iv4]

        f1 = {'featA'}
        f2 = {'featB', 'featC'}
        f3 = {'featD'}
        f4 = {'featE', 'featF', 'featG'}
        fs = [f1, f2, f3, f4]

        for iv, feat in zip(ivs, fs):
            value = feat
            theirs[iv].apply(htseq_add)
            ours[iv] += feat


if __name__ == '__main__':
    unittest.main()
