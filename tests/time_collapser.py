#!/usr/bin/env python
""" 
script to time functions used in collapser.py 

starting point to figure out where slow parts are
and where I should focus my efforts to optimize.

"""

import timeit

def time_seq_counter_from_fastq():
    SETUP_CODE = '''from aquatx.srna.collapser import seq_counter'''
    TEST_CODE = '''seq_counter('tests/testdata/KB1_cleaned.fq')'''

    times = timeit.repeat(setup=SETUP_CODE, stmt = TEST_CODE, repeat=2, number=10)
    ### [x/10 for x in times] is just dividing each time by the number of iterations being timed in order
    ### to return an average time for just ONE iteration
    print('Average time to convert fastq sequences to a counter dictionary: {}'.format(sum([x/10 for x in times])/len(times)))
    print('Max time to convert fastq sequences to a counter dictionary: {}'.format(max([x/10 for x in times])))

def time_seq_2_fasta():
    SETUP_CODE = '''
from aquatx.srna.collapser import seq_counter, seq2fasta
seqs = seq_counter('tests/testdata/KB1_cleaned.fq')'''
    TEST_CODE = '''seq2fasta(seqs, 'tests/testdata/KB1_uniq.fa')'''

    times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE, repeat=2, number=10)
    print('Average time to write collapsed fasta: {}'.format(sum([x/10 for x in times])/len(times)))
    print('Max time to write collapsed fasta: {}'.format(max([x/10 for x in times])))

def time_collapse_fastq_file_full():
    SETUP_CODE = '''
from aquatx.srna.collapser import seq_counter, seq2fasta'''
    TEST_CODE = '''
seqs = seq_counter('tests/testdata/KB1_cleaned.fq')
seq2fasta(seqs, 'tests/testdata/KB1_uniq.fa')'''

    times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE, repeat=2, number=10)
    print('Average time to collapse fastq: {}'.format(sum([x/10 for x in times])/len(times)))
    print('Max time to to collapse fastq: {}'.format(max([x/10 for x in times])))

def main():
    time_seq_counter_from_fastq()
    time_seq_2_fasta()
    time_collapse_fastq_file_full()

if __name__ == '__main__':
    main()
