""" 
script to time functions used in counter.py 

starting point to figure out where slow parts are
and where I should focus my efforts to optimize.

"""

import timeit

def time_sam_to_df():
    SETUP_CODE = '''from smrna.counter import sam_to_df'''
    TEST_CODE = '''sam_to_df('tests/testdata/KB1_test.sam')'''

    times = timeit.repeat(setup=SETUP_CODE, stmt = TEST_CODE, repeat=2, number=10)
    ### [x/10 for x in times] is just dividing each time by the number of iterations being timed in order
    ### to return an average time for just ONE iteration
    print('Average time to convert sam to a pandas dataframe: {}'.format(sum([x/10 for x in times])/len(times)))
    print('Max time to convert sam to a pandas dataframe: {}'.format(max([x/10 for x in times])))

def time_seq_counter():
    SETUP_CODE = '''
from smrna.counter import sam_to_df, seq_counter
sam_df = sam_to_df('tests/testdata/KB1_test.sam')'''

    TEST_CODE = '''seq_df = seq_counter(sam_df)'''

    times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE, repeat=2, number=10)
    print('Average time to set up sam dataframe after reading: {}'.format(sum([x/10 for x in times])/len(times)))
    print('Max time to set up sam dataframe after reading: {}'.format(max([x/10 for x in times])))

def time_nt_counter():
    SETUP_CODE = '''
from smrna.counter import sam_to_df, seq_counter, nt_counter
sam_df = sam_to_df('tests/testdata/KB1_test.sam')
seq_df = seq_counter(sam_df)
'''
    TEST_CODE = '''nt_counts = nt_counter(seq_df, save=True, out_prefix='KB1_test')'''

    times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE, repeat=2, number=10)
    print('Average time to count 5\'nt and sizes from df: {}'.format(sum([x/10 for x in times])/len(times)))
    print('Max time to count 5\'nt and sizes from df: {}'.format(max([x/10 for x in times])))

def time_feature_reader():
    SETUP_CODE='''from smrna.counter import feature_reader'''
    TEST_CODE='''features = feature_reader('tests/testdata/fixed_miRNAs_4nt.gff')'''
    
    times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE, repeat=2, number=10)
    print('Average time to read in features as dataframe: {}'.format(sum([x/10 for x in times])/len(times)))
    print('Max time to read in features as dataframe: {}'.format(max([x/10 for x in times])))

def time_feature_splitter_annotations():
    SETUP_CODE='''
from smrna.counter import feature_reader, feature_splitter
features = feature_reader('tests/testdata/fixed_miRNAs_4nt.gff')'''
    TEST_CODE='''split_feats = feature_splitter(features)'''
    
    times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE, repeat=2, number=10)
    print('Average time to split read in features: {}'.format(sum([x/10 for x in times])/len(times)))
    print('Max time to split in features: {}'.format(max([x/10 for x in times])))
    
def time_feature_splitter_both():
    SETUP_CODE='''
from smrna.counter import sam_to_df, seq_counter, nt_counter, feature_reader, feature_splitter, feature_counter
sam_df = sam_to_df('tests/testdata/KB1_test.sam')
seq_df = seq_counter(sam_df)
features = feature_reader('tests/testdata/fixed_miRNAs_4nt.gff')'''
    TEST_CODE='''
split_feats = feature_splitter(features)
split_seqs = feature_splitter(seq_df)'''

    times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE, repeat=2, number=10)
    print('Average time to split read in features & seqs: {}'.format(sum([x/10 for x in times])/len(times)))
    print('Max time to split in features & seqs: {}'.format(max([x/10 for x in times])))
    
def time_feature_counter_full():
    SETUP_CODE = '''
from smrna.counter import sam_to_df, seq_counter, nt_counter,feature_reader, feature_counter
sam_df = sam_to_df('tests/testdata/KB1_test.sam')
seq_df = seq_counter(sam_df)
features = feature_reader('tests/testdata/fixed_miRNAs_4nt.gff')
'''
    TEST_CODE = '''feature_counter(features, seq_df)'''

    times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE, repeat=2, number=10)
    print('Average time to assign counts to full features: {}'.format(sum([x/10 for x in times])/len(times)))
    print('Max time to assign counts to full features: {}'.format(max([x/10 for x in times])))

def time_feature_counter_split_both_strand():
    SETUP_CODE = '''
from smrna.counter import sam_to_df, seq_counter, nt_counter, feature_reader, feature_splitter, feature_counter
sam_df = sam_to_df('tests/testdata/KB1_test.sam')
seq_df = seq_counter(sam_df)
features = feature_reader('tests/testdata/fixed_miRNAs_4nt.gff')
split_feats = feature_splitter(features)
split_seqs = feature_splitter(seq_df)
'''
    TEST_CODE = '''
feature_counter(split_feats['+'], split_seqs['+'])
feature_counter(split_feats['-'], split_seqs['-'])'''

    times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE, repeat=2, number=10)
    print('Average time to assign counts to all features by strand: {}'.format(sum([x/10 for x in times])/len(times)))
    print('Max time to assign counts to all features by strand: {}'.format(max([x/10 for x in times])))

def time_feature_counter_split_plus_strand():
    SETUP_CODE = '''
from smrna.counter import sam_to_df, seq_counter, nt_counter, feature_reader, feature_splitter, feature_counter
sam_df = sam_to_df('tests/testdata/KB1_test.sam')
seq_df = seq_counter(sam_df)
features = feature_reader('tests/testdata/fixed_miRNAs_4nt.gff')
split_feats = feature_splitter(features)
split_seqs = feature_splitter(seq_df)
'''
    TEST_CODE = '''feature_counter(split_feats['+'], split_seqs['+'])'''

    times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE, repeat=2, number=10)
    print('Average time to assign counts to + strand features: {}'.format(sum([x/10 for x in times])/len(times)))
    print('Max time to assign counts to + strand features: {}'.format(max([x/10 for x in times])))

def time_main_routine():
    SETUP_CODE = '''
import argparse
import pandas as pd
from smrna.counter import sam_to_df, seq_counter, nt_counter, feature_reader, feature_splitter, feature_counter
args = argparse.Namespace
args.input_file='tests/testdata/KB1_test.sam'
args.ref_annotations='tests/testdata/fixed_miRNAs_4nt.gff'
args.out_prefix='KB1_test'
args.mask_file=None
args.antisense=None
'''
    TEST_CODE = '''
#Read in and process sam file
sam_df = sam_to_df(args.input_file)
seq_df = seq_counter(sam_df)

# Read in and process feature file
features = feature_reader(args.ref_annotations)
if args.mask_file:
    masks = feature_reader(args.mask_file)
    split_mask = feature_splitter(masks)
# Split by strand/chrom
split_feats = feature_splitter(features)
split_seqs = feature_splitter(seq_df)

split_seqs = {fkey: split_seqs[fkey] for fkey in split_feats.keys()}
feat_counts = pd.DataFrame()

count_s = 0
for skey in split_seqs.keys():
    if skey not in ['+', '-']:
        raise ValueError('Strands are not +/-.')
    elif args.antisense:
        if count_s >= 1:
            pass
        count_s += 1
        anti_counts = pd.DataFrame()
        anti_counts = pd.concat([anti_counts, feature_counter(split_feats['+'], split_seqs['-'])])
        anti_counts = pd.concat([anti_counts, feature_counter(split_feats['-'], split_seqs['+'])])
        anti_counts.to_csv(args.out_prefix + '_antisense_feature_counts.csv')

if args.mask_file:
    try:
        feat_counts = pd.concat([feat_counts, feature_counter(split_feats['-'], split_seqs['-'], mask=split_mask['-'])])
    # Might expect a key error since mask file doesn't have to have features on both strands
    except KeyError:
        pass
    try:
        feat_counts = pd.concat([feat_counts, feature_counter(split_feats['+'], split_seqs['+'], mask=split_mask['+'])])
    except KeyError:
        pass
else:
    feat_counts = pd.concat([feat_counts, feature_counter(split_feats['+'], split_seqs['+'])])
    feat_counts = pd.concat([feat_counts, feature_counter(split_feats['-'], split_seqs['-'])])

feat_counts.to_csv(args.out_prefix + '_feature_counts.csv')    
'''

    times = timeit.repeat(setup=SETUP_CODE, stmt = TEST_CODE, repeat=2, number=10)
    print('Average time to run main counter routine: {}'.format(sum([x/10 for x in times])/len(times)))
    print('Max time to run main counter routine: {}'.format(max([x/10 for x in times])))

def main():
    time_main_routine()
    time_feature_counter_split_plus_strand()
    time_feature_counter_split_both_strand()
    time_feature_counter_full()
    time_feature_splitter_both()
    time_feature_splitter_annotations()
    time_feature_reader()
    time_nt_counter()
    time_seq_counter()
    time_sam_to_df()

if __name__ == '__main__':
    main()
