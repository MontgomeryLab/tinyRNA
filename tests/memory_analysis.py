#!/usr/bin/env python

""" A script for analyzing memory usage of code elements """

import argparse
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_args():
    """ Get input arguments. """

    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--report', metavar='REPORT',
                        help='report file to write results to')
    parser.add_argument('-l', '--limit', type=int, help='Limit number of lines')
    parser.add_argument('-i', '--step', type=int, help='Increase lines by step size')
    parser.add_argument('-s', '--start', type=int, help='Number of lines to start with')
  
    args = parser.parse_args()

    return args

def set_up_data(rows):
    """ Create the data frame """
    chrom = ["CHROMOSOME_I", "CHROMOSOME_II", "CHROMOSOME_III", "CHROMOSOME_IV", "CHROMOSOME_V"]
    attr = ["-3p", "-5p"]
    
    # Create a dataframe of random elements that resemble elements in a real dataset
    df = pd.DataFrame({'chr': [chrom[random.choice(range(5))] for _ in range(rows)],
                      'pos': np.random.randint(0,20000, size=(rows, 1)).tolist(),
                      'seq': [(''.join(random.choice("ATCG")) for _ in range(np.random.randint(15,32))) for _ in range(rows)],
                      'length': np.random.randint(15,32, size=(rows,1)).tolist(),
                      'strand': "+",
                      'cor_counts': np.random.uniform(0,1000, [rows,1]).tolist(),
                      'feature': "mirna",
                      'start': np.random.randint(0,20000, size=(rows, 1)).tolist(),
                      'end': np.random.randint(0,20000, size=(rows, 1)).tolist(),
                      'attr': ["mir-"+str(np.random.randint(1,500))+random.choice(attr) for _ in range(rows)] })
    return df

def fetch_mem_usage(start, limit, step):
    """ 
    Analyze the memory usage of data frames increasing in
    size.
    
    Inputs:
    start = minimum number of rows to calculate memory for
    limit = maximum number of rows to calculate memory for
    step = step size of rows to calculate memory for
    """
    
    counter = [start,0,0]
    
    while (counter[0] <= limit):
        df = set_up_data(counter[0])
        mem = df.memory_usage(deep=True).sum()

        yield mem
        
        counter[0] += step
        counter[1] += 1
        counter[2] += 1

def main():
    """ main routine """
    # get inputs from user
    args = get_args()
    # calculate steps for plotting
    steps = np.arange(args.start, args.limit + args.step, args.step)
    # calculate memory usage per step TODO: don't recreate df everytime, concat eat step to go faster
    nums = np.fromiter(fetch_mem_usage(args.start, args.limit, args.step), float)
    print("Number of bytes per line: {}".format(np.mean(nums/steps)))
    # plot
    plt.plot(steps,nums)
    plt.show()

if __name__ == '__main__':
    main()
