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
import itertools

def get_args():
    """Get input arguments from the user/command line."""

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-files', metavar='DATAFILE', required=True, nargs='+',
                        help='input files with data from final merged tables')
    parser.add_argument('-o', '--out-prefix', metavar='OUTFILE', required=True,
                        help='prefix to use for output PDF files. If mode=all/by-sample, sample.'\
                        'names will also be appended to the prefix.')
    parser.add_argument('-d', '--data-types', metavar='DATATYPE', required=True, nargs='+',
                        help='List of data types corresponding to input files. Options: '\
                        '\nraw_counts: the raw counts per feature data,'\
                        '\nnorm_counts: the normalized counts per feature data,'\
                        '\ndegs: the table of differential gene expression analysis including p-values,'\
                        '\nlen_dist: the size and 5p nucleotide count matrix,'\ 
                        '\nclass_counts: the table of raw counts per class')
    parser.add_argument('-p', '--plots', metavar='PLOTS', required=True, nargs='+',
                        help='List of plots to create. Options: '\
                        '\nlen_dist: A stacked barplot showing size & 5p-nt distribution,'\
                        '\nclass_charts: A pie and barchart showing proportions of counts per class,'\
                        '\nreplicate_scatter: A scatter plot comparing replicates for all count files given,'\
                        '\nsample_avg_scatter: A scatter plot comparing all sample groups,'\
                        'averaged by replicate. Uses the normalized counts for averaging.'\
                        '\nsample_avg_scatter_by_class: A scatter plot comparing all sample groups,'\
                        'with different classes highlighted.'\
                        '\nsample_avg_scatter_by_deg: A scatter plot comparing all sample groups,'\
                        'with significantly different (padj<0.05) genes highlighted.'\
                        '\nsample_avg_scatter_by_both: A scatter plot comparing all sample groups,'\
                        'with classes and significantly different genes highlighted')
    args = parser.parse_args()

    return args

def get_pairs(samples):
    """Returns all pairs of comparisons from a list of samples to compare."""

    return itertools.product(samples, 2)

def size_dist_plot(size_file, pdf_name, **kwargs):
    """Create a PDF of size and 5'nt distribution plot for a sample.

    Args:
        size_file: Path to file containing size + 5p-nt counts
        pdf name: String to use as prefix for saving the output
        kwargs: Keyword arguments to pass to matplotlib rc or plot 
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
    """
    class_counts = pd.read_csv(class_file, index_col=0, header=None).drop('_no_class')
    fig, ax = aqplt.class_pie_barh(class_counts, pdf_name, **kwargs)
    
    # Save the plot
    fig.savefig(pdf_name, bbox_inches='tight')

def get_replicates(df):
    """Determine which replicates to compare in a dataframe based on sample_replicate_N format.
    
    Args:
        df: The dataframe with columns using sample_replicate_N format and features for rows

    Returns:
        sample_dict: A dictionary containing all column names corresponding to unique samples.
    """
    
    sample_dict = dict()

    for col in df.columns:
        sample = col.split("_replicate_")[0]
        sample_dict.setdefault(sample, []).append(col)

    return sample_dict

def get_sample_averages(df):
    """Average the counts per sample across replicates.
    
    Args:
        df: The dataframe with columns using sample_replicate_N format and features for rows
            to be averaged.
    
    Returns:
        new_df: The new dataframe with averaged counts per unique sample.
    """

    samples = get_replicates(df)
    
    new_df = pd.DataFrame(index=df.index, columns=samples.keys())

    for key, val in samples.items():
        new_df.loc[key] = df[val].mean(axis=1)
    
    return new_df

def main():
    """ 
    Main routine
    """
    args = get_args()

    # Sort files & types
    data_dict = dict(zip(args.input_files, args.data_types))

if __name__ == '__main__':
    main()

