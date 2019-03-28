#!/usr/bin/env python

""" unit tests for functions in counter.py """

import unittest
import pandas as pd
import numpy as np
import smrna.counter as smrna

class test_get_sam_flags(unittest.TestCase):
    """ 
    Testing the get_sam_flags function for returning
    the correct values, raising errors when not given ints
    and incorrect combinations of flags
    """
    def test_no_flag(self):
        self.assertEqual(smrna.get_sam_flags(0),[0])
    def test_one_flag(self):
        self.assertEqual(smrna.get_sam_flags(128), [128])
    def test_two_flags(self):
        self.assertEqual(smrna.get_sam_flags(2180),[2048,128,4])
    # add raise error when value isn't an appropriate flag value

class test_find_strand(unittest.TestCase):
    """
    Testing whether or not the find_strand function
    returns the correct value and if it raises an error
    when incorrect values are given
    """
    def test_return_plus(self):
        self.assertEqual(smrna.find_strand(2180), '+')
    def test_return_minus(self):
        self.assertEqual(smrna.find_strand(16), '-')

    # add raise error if value not a list
    # add raise error if value is empty list

class test_sam_to_df(unittest.TestCase):
    """
    Testing converting a sam file to a dataframe using a
    sam file created with certain values
    """
    def setUp(self):
        self.df = pd.DataFrame([{'qname': 'seq_1_x3', 'flag': 0, 'chr': 'CHROMOSOME_I',
                                'pos': 1738690,'seq': 'CTTCCTCATGTGCTCTGACGT'},
                                {'qname': 'seq_2_x15', 'flag': 0, 'chr': 'CHROMOSOME_V',
                                'pos': 3891144, 'seq': 'GTTTTCAGCTCCCTGACGTTTGGA'},
                                {'qname': 'seq_4_x2', 'flag': 16, 'chr': 'CHROMOSOME_IV',
                                'pos': 162, 'seq': 'AGGAATCCCCTCCATCCACACC'},
                                {'qname': 'seq_5_x25431', 'flag': 16,'chr': 'CHROMOSOME_II',
                                'pos': 7701848, 'seq': 'TCATTGAAGTTCGATCCAAC'},
                                {'qname': 'seq_5_x25431', 'flag': 0, 'chr': 'CHROMOSOME_II',
                                'pos': 770184, 'seq': 'TCATTGAAGTTCGATCCAAC'},
                                {'qname': 'seq_7_x1', 'flag': 2048, 'chr': 'CHROMOSOME_III',
                                'pos': 5017450, 'seq': 'GTTATGTGGTGGATGTGCCTT'}], 
                                columns=['qname','flag','chr','pos','seq'])
        self.df['flag'] = self.df['flag'].astype(np.int32)
        self.df['chr'] = self.df['chr'].astype('category')
        self.df['pos'] = self.df['pos'].astype(np.int32)

    def tearDown(self):
        del self.df
        
    def test_sam_df_content(self):
        pd.testing.assert_frame_equal(self.df, smrna.sam_to_df('tests/unit_test_data/good_test.sam'))
    
class test_seq_counter(unittest.TestCase):
    """
    Testing initial processing of the sam dataframe
    """
    def setUp(self):
        self.df = pd.DataFrame([{'qname': 'seq_1_x3', 'flag': 0, 'chr': 'CHROMOSOME_I',
                                'pos': 1738690, 'seq': 'CTTCCTCATGTGCTCTGACGT', 'length': 21,
                                'counts': 3, '5p_nt': 'C', 'strand': '+', 'num_hits': 1, 
                                'cor_counts': 3},
                                {'qname': 'seq_2_x15', 'flag': 0, 'chr': 'CHROMOSOME_V',
                                'pos': 3891144, 'seq': 'GTTTTCAGCTCCCTGACGTTTGGA','length': 24,
                                'counts': 15, '5p_nt': 'G', 'strand': '+', 'num_hits': 1,
                                'cor_counts': 15},
                                {'qname': 'seq_4_x2', 'flag': 16, 'chr': 'CHROMOSOME_IV',
                                'pos': 162, 'seq': 'AGGAATCCCCTCCATCCACACC','length': 22,
                                'counts': 2,'5p_nt': 'A','strand': '-','num_hits': 1,
                                'cor_counts': 2},
                                {'qname': 'seq_5_x25431', 'flag': 16,'chr': 'CHROMOSOME_II',
                                'pos': 7701848,'seq': 'TCATTGAAGTTCGATCCAAC','length': 20,
                                'counts': 25431,'5p_nt': 'T','strand': '-','num_hits': 2,
                                'cor_counts': 12715.5},
                                {'qname': 'seq_5_x25431', 'flag': 0,'chr': 'CHROMOSOME_II',
                                'pos': 770184,'seq': 'TCATTGAAGTTCGATCCAAC','length': 20,
                                'counts': 25431,'5p_nt': 'T','strand': '+','num_hits': 2,
                                'cor_counts': 12715.5},
                                {'qname': 'seq_7_x1', 'flag': 2048,'chr': 'CHROMOSOME_III',
                                'pos': 5017450,'seq': 'GTTATGTGGTGGATGTGCCTT', 'length': 21,
                                'counts': 1,'5p_nt': 'G','strand': '+','num_hits': 1,
                                'cor_counts': 1} ],
                                columns=['qname','flag','chr','pos','seq','length','counts',
                                        '5p_nt','strand','num_hits','cor_counts'])
                      
        self.df['flag'] = self.df['flag'].astype(np.int32)
        self.df['chr'] = self.df['chr'].astype('category')
        self.df['pos'] = self.df['pos'].astype(np.int32)
        self.df['length'] = self.df['length'].astype(np.int32)
        self.df['counts'] = self.df['counts'].astype(np.int32)
        self.df['5p_nt'] = self.df['5p_nt'].astype('category')
        self.df['strand'] = self.df['strand'].astype('category')
        self.df['cor_counts'] = self.df['cor_counts'].astype(np.float32)
        self.df['num_hits'] = self.df['num_hits'].astype(np.int32)
        self.testdf = smrna.sam_to_df('tests/unit_test_data/good_test.sam')
    
    def tearDown(self):
        del self.df, self.testdf

    def test_sam_additions_to_df(self):
        pd.testing.assert_frame_equal(self.df, smrna.seq_counter(self.testdf))

class test_nt_counter(unittest.TestCase):
    """
    Testing if the nucleotide/length counter returns
    expected values based on a generated input
    """
    def setUp(self):
        self.df = pd.DataFrame([{'length': 20, 'A': 0, 'C': 0, 'G': 0, 'T': 25431}, 
                                {'length': 21, 'A': 0, 'C': 3, 'G': 1, 'T':0},
                                {'length': 22, 'A': 2, 'C': 0, 'G': 0, 'T': 0},
                                {'length': 24, 'A': 0, 'C': 0, 'G': 15, 'T': 0}],
                                 columns=['length','A','C','G','T'], dtype=np.int32).set_index('length')
        self.df.columns = pd.CategoricalIndex(['A', 'C', 'G', 'T'], categories=['A', 'C', 'G', 'T'], ordered=False, name='5p_nt', dtype='category')
        self.testdf = smrna.sam_to_df('tests/unit_test_data/good_test.sam')
        self.testdf = smrna.seq_counter(self.testdf)
    def tearDown(self):
        del self.df, self.testdf
    def test_nt_count_df(self):
        pd.testing.assert_frame_equal(self.df, smrna.nt_counter(self.testdf, save=False))

class test_feature_reader(unittest.TestCase):
    """ 
    Testing if a gff file is properly read in
    with expected headers
    """
    def setUp(self):
        self.df = pd.DataFrame([{'chr': 'CHROMOSOME_I', 'source':'.', 'feature': 'miRNA',
                                'start': 1738648, 'end': 1738679, 'score':'.', 'strand':'+', 
                                'frame': '.', 'attr':'cel-miR-50-5p'},
                                {'chr': 'CHROMOSOME_I', 'source':'.', 'feature': 'miRNA',
                                'start': 1738690, 'end': 1738719, 'score':'.', 'strand':'+',
                                'frame': '.', 'attr':'cel-miR-50-3p'},
                                {'chr': 'CHROMOSOME_I', 'source':'.', 'feature': 'miRNA',
                                'start': 2888510, 'end': 2888540, 'score':'.', 'strand':'-',
                                'frame': '.', 'attr':'cel-miR-5546-5p'},
                                {'chr': 'CHROMOSOME_I', 'source':'.', 'feature': 'miRNA',
                                'start': 2888468, 'end': 2888497, 'score':'.', 'strand':'-',
                                'frame': '.', 'attr':'cel-miR-5546-3p'}],
                                columns=['chr','source','feature','start','end','score','strand','frame','attr'])

    def tearDown(self):
        del self.df

    def test_read_features(self):
        pd.testing.assert_frame_equal(self.df, smrna.feature_reader('tests/unit_test_data/good_features.gff'))

class test_feature_splitter(unittest.TestCase):
    """
    Testing if a feature data frame or sequence data frame
    is split properly
    """
    def setUp(self):
        self.df = pd.DataFrame([{'chr': 'CHROMOSOME_I', 'source':'.', 'feature': 'miRNA',
                                'start': 1738648, 'end': 1738679, 'score':'.', 'strand':'+',
                                'frame': '.', 'attr':'cel-miR-50-5p'},
                                {'chr': 'CHROMOSOME_I', 'source':'.', 'feature': 'miRNA',
                                'start': 1738690, 'end': 1738719, 'score':'.', 'strand':'+',
                                'frame': '.', 'attr':'cel-miR-50-3p'},
                                {'chr': 'CHROMOSOME_I', 'source':'.', 'feature': 'miRNA',
                                'start': 2888510, 'end': 2888540, 'score':'.', 'strand':'-',
                                'frame': '.', 'attr':'cel-miR-5546-5p'},
                                {'chr': 'CHROMOSOME_I', 'source':'.', 'feature': 'miRNA',
                                'start': 2888468, 'end': 2888497, 'score':'.', 'strand':'-',
                                'frame': '.', 'attr':'cel-miR-5546-3p'}],
                                columns=['chr','source','feature','start','end','score','strand','frame','attr'])
        self.dfdict = {'+': self.df.iloc[0:2], '-': self.df.iloc[2:]}

        self.testdf = smrna.feature_reader('tests/unit_test_data/good_features.gff')

    def tearDown(self):
        del self.df, self.testdf

    def test_strand_split(self):
        pd.testing.assert_frame_equal(self.dfdict['+'], smrna.feature_splitter(self.testdf)['+'])
        pd.testing.assert_frame_equal(self.dfdict['-'], smrna.feature_splitter(self.testdf)['-'])

class test_assign_counts(unittest.TestCase):
    """
    Testing if counts are assigned properly to a feature
    """
    def setUp(self):
        self.feat = pd.DataFrame({'chr': 'CHROMOSOME_I', 'feature': 'miRNA',
                    'start': 1738648, 'end': 1738679, 'strand':'+',
                    'frame': '.', 'attr':'cel-miR-50-5p'}, index=[0])
        self.df = pd.DataFrame([{'qname': 'seq_1_x3', 'flag': 0, 'chr': 'CHROMOSOME_I',
                                'pos': 1738690, 'seq': 'CTTCCTCATGTGCTCTGACGT', 'length': 21,
                                'counts': 3, '5p_nt': 'C', 'strand': '+', 'num_hits': 1, 
                                'cor_counts': 3},
                                {'qname': 'seq_2_x15', 'flag': 0, 'chr': 'CHROMOSOME_V',
                                'pos': 3891144, 'seq': 'GTTTTCAGCTCCCTGACGTTTGGA','length': 24,
                                'counts': 15, '5p_nt': 'G', 'strand': '+', 'num_hits': 1,
                                'cor_counts': 15},
                                {'qname': 'seq_4_x2', 'flag': 16, 'chr': 'CHROMOSOME_IV',
                                'pos': 162, 'seq': 'AGGAATCCCCTCCATCCACACC','length': 22,
                                'counts': 2,'5p_nt': 'A','strand': '-','num_hits': 1,
                                'cor_counts': 2},
                                {'qname': 'seq_5_x25431', 'flag': 16,'chr': 'CHROMOSOME_II',
                                'pos': 7701848,'seq': 'TCATTGAAGTTCGATCCAAC','length': 20,
                                'counts': 25431,'5p_nt': 'T','strand': '-','num_hits': 2,
                                'cor_counts': 12715.5},
                                {'qname': 'seq_5_x25431', 'flag': 0,'chr': 'CHROMOSOME_II',
                                'pos': 770184,'seq': 'TCATTGAAGTTCGATCCAAC','length': 20,
                                'counts': 25431,'5p_nt': 'T','strand': '+','num_hits': 2,
                                'cor_counts': 12715.5},
                                {'qname': 'seq_7_x1', 'flag': 2048,'chr': 'CHROMOSOME_III',
                                'pos': 5017450,'seq': 'GTTATGTGGTGGATGTGCCTT', 'length': 21,
                                'counts': 1,'5p_nt': 'G','strand': '+','num_hits': 1,
                                'cor_counts': 1} ],
                                columns=['qname','flag','chr','pos','seq','length','counts',
                                        '5p_nt','strand','num_hits','cor_counts'])                     
    def tearDown(self):
        del self.feat, self.df
    def test_assigned_good_count(self):
        # not sure the best way to test/set this up so I may end up not testing it
        # or changing the function because apparently I don't understand fully what
        # I am doing :)
        print('unimplemented test. will finish this later')

class test_feature_counter(unittest.TestCase):
    """
    Test whether or not the feature counter returns expected
    counts based on a made up file
    """
    def setUp(self):
        self.testdf = pd.DataFrame([{'qname': 'seq_1_x3', 'flag': 0, 'chr': 'CHROMOSOME_I',
                                'pos': 1738690, 'seq': 'CTTCCTCATGTGCTCTGACGT', 'length': 21,
                                'counts': 3, '5p_nt': 'C', 'strand': '+', 'num_hits': 1, 
                                'cor_counts': 3},
                                {'qname': 'seq_2_x15', 'flag': 0, 'chr': 'CHROMOSOME_V',
                                'pos': 3891144, 'seq': 'GTTTTCAGCTCCCTGACGTTTGGA','length': 24,
                                'counts': 15, '5p_nt': 'G', 'strand': '+', 'num_hits': 1,
                                'cor_counts': 15},
                                {'qname': 'seq_4_x2', 'flag': 16, 'chr': 'CHROMOSOME_IV',
                                'pos': 162, 'seq': 'AGGAATCCCCTCCATCCACACC','length': 22,
                                'counts': 2,'5p_nt': 'A','strand': '-','num_hits': 1,
                                'cor_counts': 2},
                                {'qname': 'seq_5_x25431', 'flag': 16,'chr': 'CHROMOSOME_II',
                                'pos': 7701848,'seq': 'TCATTGAAGTTCGATCCAAC','length': 20,
                                'counts': 25431,'5p_nt': 'T','strand': '-','num_hits': 2,
                                'cor_counts': 12715.5},
                                {'qname': 'seq_5_x25431', 'flag': 0,'chr': 'CHROMOSOME_II',
                                'pos': 770184,'seq': 'TCATTGAAGTTCGATCCAAC','length': 20,
                                'counts': 25431,'5p_nt': 'T','strand': '+','num_hits': 2,
                                'cor_counts': 12715.5},
                                {'qname': 'seq_7_x1', 'flag': 2048,'chr': 'CHROMOSOME_III',
                                'pos': 5017450,'seq': 'GTTATGTGGTGGATGTGCCTT', 'length': 21,
                                'counts': 1,'5p_nt': 'G','strand': '+','num_hits': 1,
                                'cor_counts': 1} ],
                                columns=['qname','flag','chr','pos','seq','length','counts',
                                        '5p_nt','strand','num_hits','cor_counts'])                     
        
        self.testfeat = pd.DataFrame([{'chr': 'CHROMOSOME_I', 'source':'.', 'feature': 'miRNA',
                                'start': 1738648, 'end': 1738679, 'score':'.', 'strand':'+',
                                'frame': '.', 'attr':'cel-miR-50-5p'},
                                {'chr': 'CHROMOSOME_I', 'source':'.', 'feature': 'miRNA',
                                'start': 1738690, 'end': 1738719, 'score':'.', 'strand':'+',
                                'frame': '.', 'attr':'cel-miR-50-3p'},
                                {'chr': 'CHROMOSOME_I', 'source':'.', 'feature': 'miRNA',
                                'start': 2888510, 'end': 2888540, 'score':'.', 'strand':'-',
                                'frame': '.', 'attr':'cel-miR-5546-5p'},
                                {'chr': 'CHROMOSOME_I', 'source':'.', 'feature': 'miRNA',
                                'start': 2888468, 'end': 2888497, 'score':'.', 'strand':'-',
                                'frame': '.', 'attr':'cel-miR-5546-3p'}],
                                columns=['chr','source','feature','start','end','score','strand','frame','attr'])
        self.df = pd.DataFrame([{'chr': 'CHROMOSOME_I', 'feature': 'miRNA',
                                'start': 1738648, 'end': 1738679, 'strand':'+',
                                 'attr':'cel-miR-50-5p','counts':0.0},
                                {'chr': 'CHROMOSOME_I', 'feature': 'miRNA',
                                'start': 1738690, 'end': 1738719, 'strand':'+',
                                 'attr':'cel-miR-50-3p','counts':3.0},
                                {'chr': 'CHROMOSOME_I', 'feature': 'miRNA',
                                'start': 2888510, 'end': 2888540, 'strand':'-',
                                 'attr':'cel-miR-5546-5p', 'counts':0.0},
                                {'chr': 'CHROMOSOME_I', 'feature': 'miRNA',
                                'start': 2888468, 'end': 2888497, 'strand':'-',
                                 'attr':'cel-miR-5546-3p', 'counts':0.0}],
                                columns=['chr','feature','start','end','strand','attr','counts'])
    def tearDown(self):
        del self.df, self.testfeat, self.testdf
    def test_good_feature_count(self):
        pd.testing.assert_frame_equal(self.df, smrna.feature_counter(self.testfeat, self.testdf))

if __name__ == '__main__':
    unittest.main()
