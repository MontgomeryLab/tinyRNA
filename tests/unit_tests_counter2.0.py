import unittest
from unittest.mock import patch, mock_open

import HTSeq

import aquatx.srna.counter as counter

class MyTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.sam = "11_count=1	16	I	1043730	255	23M	*	0	0	AATCTACTCGGAAGATCATAATC	IIIIIIIIIIIIIIIIIIIIIII	XA:i:0	MD:Z:23	NM:i:0	XM:i:2"

    def test_sam_reader(self):
        samfile = HTSeq.BAM_Reader("na")
        #with patch('HTSeq.pysam.libcalignmentfile.AlignmentFile._open', new_callable=mock_open(read_data=self.sam)):
        for alignment in samfile:
            print(alignment)


if __name__ == '__main__':
    unittest.main()
