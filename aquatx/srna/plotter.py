""" Plotting script for the sRNA AQuATx workflow.

This script produces basic static plots for publications as part of the AQuATx
workflow. It uses a default style and produces multiple PDFs. This script is
intended to be used immediately following the aquatx-count/counter.py step of
the workflow. It creates a specific set of plots through the mode argument.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import aquatx.srna.plotterlib as aqplt

def get_args():
    """Get input arguments from the user/command line."""

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-files', metavar='DATAFILE', required=True, nargs='+',
                        help='input files with data from final merged tables')
    parser.add_argument('-o', '--out-prefix', metavar='OUTFILE', required=True,
                        help='prefix to use for output PDF files. If mode=all/by-sample, sample.'\
                        'names will also be appended to the prefix.')
    parser.add_argument('-m', '--data-types', metavar='MODE', required=True, nargs='+',
                        help='List of data types corresponding to input files. Options: '\
                        'raw_counts, norm_counts, degs, len_dist, class_counts')
    args = parser.parse_args()

    return args

def size_dist_plot(size_file, pdf_name, **kwargs):
    """
    Create a PDF of size and 5'nt distribution plot for a sample.

    Args:
        size_file: Path to file containing size + 5p-nt counts
        pdf name: String to use as prefix for saving the output
        kwargs: Keyword arguments to pass to matplotlib rc or plot 

    Returns:
        None: Saves a PDF file containing the plot
    """
    # Read the size_dist file
    size_dist = pd.read_csv(size_file, index_col=0)
    
    # Create the plot
    nt_colors = ['#F78E2D', '#CBDC3F', '#4D8AC8', '#E06EAA'] # Orange, Yellow-green, Blue, Pink
    ax = aqplt.size_dist_bar(size_dist_rpm, color=nt_colors, **kwargs)
    
    # Save the plot
    ax.savefig(out_name, bbox_inches='tight')

def class_plots(class_file, pdf_name, **kwargs):
    """Create a PDF of the proportion of counts assigned to
    a feature in a particular class as a pie and bar chart.

    Args:
        class_counts: Path to the file containing a table of class counts
        pdf_name: String to use as a prefix for saving the output
        kwargs: Keyword arguments to pass to matplotlib rc or plot 

    Returns:
        None: Saves a PDF file containing the plot
    """
    class_counts = pd.read_csv(class_file, index_col=0).drop('_no_class')
    fig, ax = aqplt.class_pie_barh(class_counts, pdf_name, **kwargs)
    
    # Save the plot
    fig.savefig(pdf_name, bbox_inches='tight')

def main():
    """ 
    Main routine
    """
    args = get_args()

    # Create all plots by reading in a list of files & types
    for infile, dtype in zip(args.input_files, args.data_types):
        if dtype == 'raw_counts':
            continue
        elif dtype == 'norm_counts':
            continue
        elif dtype == 'degs':
            continue
        elif dtype == 'len_dist':
            size_5p_nt_dist_barplot(infile, args.out_prefix + '_len_dist')
        elif dtype == 'class_counts':
            continue

if __name__ == '__main__':
    main()

