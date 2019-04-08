""" Plotting functions for small RNA data. 

This script produces basic static plots for publications. It uses a default style
called smrna-light, which is a style sheet we defined for this tool. Other styles
may be used if desired. This is the master script for all plots from the QC step to
the final DEG plots."""

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings(action='once')

def get_args():
    """
    Get input arguments from the user/command line.

    Requires the user provide the output of counter.py.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-files', metavar='DATAFILE', required=True,
                        help='input files with data from final merged tables')
    parser.add_argument('-o', '--out-prefix', metavar='OUTFILE', required=True,
                        help='prefix to use for output PDF files. If mode=all/by-sample, sample.'\
                        'names will also be appended to the prefix.')
    parser.add_argument('-m', '--data-types', metavar='MODE', required=True,
                        help='List of data types corresponding to input files. Options: '\
                        'raw_counts, norm_counts, degs, len_dist, class_counts')
    parser.add_argument('-s', '--style', metavar='PLOTSTYLE', default='smrna-light',
                        help='plotting style to use. Default: smrna-light.')
    parser.add_argument('-c', '--color', metavar='PLOTCOLOR', 
                        help='color palette to use for data instead of default')

    args = parser.parse_args()

    return args

def size_5p_nt_dist_barplot(size_dist_file, prefix):
    """
    Generate the size and 5'nt distribution plot for a sample.

    Input:
        size_dist_file = the file containing the 5p nt vs length matrix
        outprefix = prefix to use for saving the file

    Output:
        PDF file containing the plot
    """
    # Read the size_dist file
    size_dist = pd.read_csv(size_dist_file, index_col=0)
    total_reads = size_dist.sum().sum()

    size_dist_rpm = size_dist/total_reads

    # Create the plot
    #fig, ax = plt.subplots()
    nt_colors = ['#F78E2D', '#CBDC3F', '#4D8AC8', '#E06EAA'] # Orange, Yellow-green, Blue, Pink
    ax = size_dist_rpm.plot(kind='bar', stacked=True, color=nt_colors,
                            title='Distribution of aligned reads')

    # Label the plot
    plt.xticks(rotation=0)
    plt.ylabel('Proportion of Reads')
    plt.xlabel('Length of Sequence')

    # Save the plot
    plt.savefig(out_name, bbox_inches='tight')

def counts_scatter_plot(infile, outprefix):
    """
    Generates scatter plot of normalized small RNA counts.
    
    Inputs:
        x_sample: counts for x axis
        y_sample: counts for y axis
        out_name: filename to use for output file

    Outputs:
        PDF file containing the plot
    """
    # TODO: create all possible combinations of scatter plots based on columns
    x_data = pd.read_csv(x_sample, index_col=0, sep='\t', usecols=[1,2], names=['sample_1_feats', 'sample_1_counts'])
    y_data = pd.read_csv(y_sample, index_col=0, sep='\t', usecols=[1,2], names=['sample_2_feats', 'sample_2_counts'])
    
    xy = pd.merge(x_data.apply(np.log), y_data.apply(np.log), how='outer', right_index=True, left_index=True).fillna(0).drop('_no_feature')
    print(xy.head())
    ax = xy.plot(kind='scatter', x='sample_1_counts', y='sample_2_counts', color='#333333',
                  title='Reads in Sample 1 vs Sample 2')
    plt.xticks(rotation=0)
    plt.ylabel('Normalized Reads (log)')
    plt.xlabel('Normalized Reads (log)')

    plt.savefig(out_name, bbox_inches='tight')

def class_plots(class_file, out_name):
    """
    Generates pie and bar charts for each library according to counts
    that correspond to different classes.
    """
    # TODO: Create all class plots based on columns (looped)
    class_counts = pd.read_csv(class_file, index_col=0, usecols=[1,2], header=0, names=['class', 'counts']).drop('_no_class')
    class_counts = class_counts/(class_counts.sum())
    print(class_counts)

    fig, axes = plt.subplots(nrows=1, ncols=2)

    class_counts.plot(kind='pie', subplots=True, ax=axes[0], labels=['' for _ in class_counts.index])
    axes[0].set_aspect("equal")
    axes[0].legend(loc='best', bbox_to_anchor= (1,0.5), fontsize=10, labels=class_counts.index)
    axes[0].set_ylabel('')
    axes[0].set_xlabel('')
    
    axes[1].barh(class_counts.index, class_counts['counts'], color=plt.rcParams['axes.prop_cycle'].by_key()['color'])
    axes[1].set_xlabel('Proportion of reads')
    axes[1].set_xlim(0,1)
    
    fig.suptitle("Proportion of classes of small RNAs", fontsize=22)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.subplots_adjust(top=0.85)
    plt.savefig(out_name, bbox_inches='tight')

def main():
    """ 
    Main routine
    """
    args = get_args()

    # Define basic plot parameters
    try:
        plt.style.use(args.style)
    except:
        print("%s theme not installed. Using seaborn-poster instead." % args.style)
        plt.style.use('seaborn-poster')

    # Create all plots by reading in a list of files & types
    for infile, dtype in zip(args.input_files, args.data_types):
        if dtype == 'raw_counts':
            counts_scatter_plot(infile, args.out_prefix + '_raw_count_scatter')
        elif dtype == 'norm_counts':
            counts_scatter(infile, args.out_prefix + '_norm_count_scatter')
        elif dtype == 'degs':
            continue
        elif dtype == 'len_dist':
            size_5p_nt_dist_barplot(infile, args.out_prefix + '_len_dist')
        elif dtype == 'class_counts':
            class_plots(infile, args.out_prefix + '_class_plots')


if __name__ == '__main__':
    main()


