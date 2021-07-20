""" Plotting script for the sRNA AQuATx workflow.

This script produces basic static plots for publications as part of the AQuATx
workflow. It uses a default style and produces multiple PDFs. This script is
intended to be used immediately following the aquatx-count/counter.py step of
the workflow. It creates a specific set of plots through the mode argument.
"""

import argparse
import os.path
import itertools
from collections import defaultdict
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import aquatx.srna.plotterlib as aqplt

def get_args():
    """Get input arguments from the user/command line."""

    parser = argparse.ArgumentParser()
    counter_files = parser.add_argument_group("Input files produced by Counter")
    diffexp_files = parser.add_argument_group("Input files produced by DEG")

    # Single file inputs
    counter_files.add_argument('-rc', '--raw-counts', metavar='RAW_COUNTS',
                               help='The ...feature_counts.csv file.')
    diffexp_files.add_argument('-nc', '--norm-counts', metavar='NORM_COUNTS',
                               help='The ...norm_counts.csv file.')

    # Multi-file inputs
    diffexp_files.add_argument('-deg', '--deg-tables', metavar='DEG_COMPARISONS,...', nargs='+',
                               help='The ...cond1...cond2...deseq.csv files')
    counter_files.add_argument('-len', '--len-dist', metavar='5P_LEN_DISTS,...', nargs='+',
                               help='The ...nt_len_dist.csv files')

    # Outputs options
    parser.add_argument('-o', '--out-prefix', metavar='OUTFILE', default='',
                        help='Optional prefix to use for output PDF files.')

    parser.add_argument('-p', '--plots', metavar='PLOTS', required=True, nargs='+',
                        help='List of plots to create. Options: \n'
                             'len_dist: A stacked barplot showing size & 5p-nt distribution,\n'
                             'class_charts: A pie and barchart showing proportions of counts per class,\n'
                             'replicate_scatter: A scatter plot comparing replicates for all count files given,\n'
                             'sample_avg_scatter: A scatter plot comparing all sample groups,'
                                'averaged by replicate. Uses the normalized counts for averaging.\n'
                             'sample_avg_scatter_by_class: A scatter plot comparing all sample groups,'
                                'with different classes highlighted.\n'
                             'sample_avg_scatter_by_deg: A scatter plot comparing all sample groups,'
                                'with significantly different (padj<0.05) genes highlighted.\n'
                             'sample_avg_scatter_by_both: A scatter plot comparing all sample groups,'
                                'with classes and significantly different genes highlighted')
    args = parser.parse_args()

    return args

def get_pairs(samples):
    """Get pairs of comparisons from a list of samples to compare.
    
    Args:
        samples: A list of samples to compare
    
    Returns:
        pair: An generator of all possible pairs
    """
    
    for pair in itertools.combinations(samples, 2):
        yield pair

def size_dist_plot(size_file, pdf_name, **kwargs):
    """Create a PDF of size and 5'nt distribution plot for a sample.

    Args:
        size_file: Path to file containing size + 5p-nt counts
        pdf_name: Filename for output PDF plot
        kwargs: Keyword arguments to pass to matplotlib rc or plot 
    """
    # Read the size_dist file
    size_dist = pd.read_csv(size_file, index_col=0)
    
    # Create the plot
    nt_colors = ['#F78E2D', '#CBDC3F', '#4D8AC8', '#E06EAA'] # Orange, Yellow-green, Blue, Pink
    ax = aqplt.size_dist_bar(size_dist, color=nt_colors, **kwargs)
    
    # Save the plot
    plt.savefig(pdf_name, bbox_inches='tight')

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

def get_sample_rep_dict(df):
    """Determine which replicates to compare in a dataframe based on sample_rep_N format.
    
    Args:
        df: The dataframe with columns using sample_rep_N format and features for rows

    Returns:
        sample_dict: A dictionary containing all column names corresponding to unique samples.
    """
    
    sample_dict = defaultdict(list)

    for col in df.columns:
        if col == "Feature.Class": continue
        sample = col.split("_rep_")[0]
        sample_dict[sample].append(col)

    return sample_dict

def get_sample_averages(df, samples):
    """Average the counts per sample across replicates.
    
    Args:
        df: The dataframe with columns using sample_rep_N format and features for rows
            to be averaged.
        samples: A dictionary containing sample names and their associated "sample_rep_N" replicates

    Returns:
        new_df: The new dataframe with averaged counts per unique sample.
    """
    
    new_df = pd.DataFrame(index=df.index, columns=samples.keys())

    for key, val in samples.items():
        new_df.loc[:, key] = df[val].mean(axis=1)
    
    return new_df


def scatter_replicates(count_df, output_prefix, samples, norm=False):
    """Creates PDFs of all pairwise comparison scatter plots from a count table.
    
    Args:
        count_df: A dataframe of counts per features
        output_prefix: A string to use as a prefix for saving files
        samples: A dictionary containing sample names and their associated "sample_rep_N" replicates
        norm: Boolean indicating if data was normalized
    """

    for samp, reps in samples.items():
        for pair in get_pairs(reps):
            rscat = aqplt.scatter_simple(count_df.loc[:,pair[0]], count_df.loc[:,pair[1]], 
                                        color='#888888', marker='s', alpha=0.5, s=50, 
                                        edgecolors='none', log_norm=True)
            rscat.set_title(samp)
            rscat.set_xlabel('Replicate ' + pair[0].split('_rep_')[1])
            rscat.set_ylabel('Replicate ' + pair[1].split('_rep_')[1])
            pdf_name = '_'.join([output_prefix, samp, 'replicates', pair[0].split('_rep_')[1],
                                pair[1].split('_rep_')[1], 'scatter.pdf'])
            rscat.figure.savefig(pdf_name, bbox_inches='tight')

def get_degs(comparisons):
    """Create a new dataframe containing all features and pvalues for each comparison
    from a set of DEG tables. 

    Args:
        comparisons: Files (DEG tables) to process

    Returns:
        de_table: A single table of p-values per feature/comparison
    """
    count = 0
    samples = {}
    
    for degfile in comparisons:
        degs = pd.read_csv(degfile, index_col=0)
        name_split = os.path.basename(degfile).split('_')
        # Negative indexes are used since prefixes are allowed to contain underscores
        samples[degfile] = name_split[-5] + "vs" + name_split[-3]

        if count == 0:
            de_table = pd.DataFrame(index=degs.index, columns=comparisons)

        de_table.loc[:, degfile] = degs.loc[:, 'padj'] # Are you sure? My hunch is this is supposed to be p-value
        count += 1
        
    de_table.rename(samples, axis=1, inplace=True)
    
    return de_table

def scatter_samples(count_df, output_prefix, classes=None, degs=None, unknown=False, pval=0.05):
    """Creates PDFs of all pairwise comparison scatter plots from a count table.
    Can highlight classes and/or differentially expressed genes as different colors.

    Args:
        count_df: A dataframe of counts per feature
        output_prefix: A string to use as a prefix for saving files
        classes: A dataframe containing class per feature
        degs: A dataframe of differential gene table output to highlight
        unknown: Boolean indicating if ambiguous classes should be highlighted
    """
    samples = get_pairs(list(count_df.columns))

    if (classes is not None) and (degs is not None):
        uniq_classes = list(pd.unique(classes.loc[:,'class'].dropna()))
        
        if not unknown:
            uniq_classes.remove('unknown')
        
        for pair in samples:
            grp_args = []
            try:
                samp_comp = '%s_vs_%s' % (pair[0], pair[1])
                deg_list = list(degs.loc[degs[samp_comp] < pval].index)
            except KeyError:
                try: 
                    samp_comp = '%s_vs_%s' % (pair[1], pair[0])
                    deg_list = list(degs.loc[degs[samp_comp] < pval].index)
                except KeyError as ke:
                    print('Sample names in count data frame do not match sample comparisons')
                    print('in DEG table. Make sure formatting is correct in your tables for')
                    print('AQuATx. Error occurred with %s, %s in count table.' % (pair[0], pair[1]))
                    raise ke
            
            class_degs = classes.loc[deg_list, :]

            for cls in uniq_classes:
                grp_args.append(list(class_degs.loc[class_degs['class'] == cls].index))

            labels = ['p > %6.2f' % pval] + uniq_classes
            sscat = aqplt.scatter_grouped(count_df.loc[:,pair[0]], count_df.loc[:,pair[1]],
                                          *grp_args, log_norm=True, labels=labels) 
            sscat.set_title('%s vs %s' % (pair[0], pair[1]))
            sscat.set_xlabel(pair[0])
            sscat.set_ylabel(pair[1])
            pdf_name = '_'.join([output_prefix, pair[0], 'vs', pair[1], 'scatter_by_deg_class.pdf'])
            sscat.figure.savefig(pdf_name, bbox_inches='tight')

    elif classes is not None:
        uniq_classes = list(pd.unique(classes.loc[:, 'Feature.Class'].dropna()))
        
        if not unknown:
            uniq_classes.remove('unknown')
        
        grp_args = []
        for cls in uniq_classes:
            grp_args.append(list(classes.loc[classes['Feature.Class'] == cls].index))
    
        for pair in samples:
            sscat = aqplt.scatter_grouped(count_df.loc[:,pair[0]], count_df.loc[:,pair[1]],
                                          *grp_args, log_norm=True, labels=uniq_classes) 
            sscat.set_title('%s vs %s' % (pair[0], pair[1]))
            sscat.set_xlabel(pair[0])
            sscat.set_ylabel(pair[1])
            pdf_name = '_'.join([output_prefix, pair[0], 'vs', pair[1], 'scatter_by_class.pdf'])
            sscat.figure.savefig(pdf_name, bbox_inches='tight')
    
    elif degs is not None:
        
        for pair in samples:    
            grp_args = []
            try:
                samp_comp = '%s_vs_%s' % (pair[0], pair[1])
                grp_args.append(list(degs.loc[degs[samp_comp] < pval].index))
            except KeyError:
                try: 
                    samp_comp = '%s_vs_%s' % (pair[1], pair[0])
                    grp_args.append(list(degs.loc[degs[samp_comp] < pval].index))
                except KeyError as ke:
                    print('Sample names in count data frame do not match sample comparisons')
                    print('in DEG table. Make sure formatting is correct in your tables for')
                    print('AQuATx. Error occurred with %s, %s in count table.' % (pair[0], pair[1]))
                    raise ke
            labels = ['p < %d' % pval]
            sscat = aqplt.scatter_grouped(count_df.loc[:,pair[0]], count_df.loc[:,pair[1]],
                                          *grp_args, log_norm=True, labels=labels) 
            sscat.set_title('%s vs %s' % (pair[0], pair[1]))
            sscat.set_xlabel(pair[0])
            sscat.set_ylabel(pair[1])
            pdf_name = '_'.join([output_prefix, pair[0], 'vs', pair[1], 'scatter_by_deg.pdf'])
            sscat.figure.savefig(pdf_name, bbox_inches='tight')
    
    else:
        for pair in samples:
            sscat = aqplt.scatter_simple(count_df.loc[:,pair[0]], count_df.loc[:,pair[1]], 
                                        color='#888888', marker='s', alpha=0.5, s=50, 
                                        edgecolors='none', log_norm=True)
            sscat.set_title('%s vs %s' % (pair[0], pair[1]))
            sscat.set_xlabel(pair[0])
            sscat.set_ylabel(pair[1])
            pdf_name = '_'.join([output_prefix, pair[0], 'vs', pair[1], 'scatter.pdf'])
            sscat.figure.savefig(pdf_name, bbox_inches='tight')


def load_raw_counts(args) -> (pd.DataFrame, Optional[dict]):
    if args.raw_counts is not None:
        raw_count_df = pd.read_csv(args.raw_counts, index_col=0)
        sample_rep_dict = get_sample_rep_dict(raw_count_df)
    else:
        raw_count_df, sample_rep_dict = None, None

    return raw_count_df, sample_rep_dict


def load_norm_counts(args) -> (pd.DataFrame, pd.DataFrame, Optional[dict]):
    if args.norm_counts is not None:
        norm_count_df = pd.read_csv(args.norm_counts, index_col=0)
        sample_rep_dict = get_sample_rep_dict(norm_count_df)  # Refactor
        norm_count_avg_df = get_sample_averages(norm_count_df, sample_rep_dict)

        classes = pd.DataFrame(norm_count_df["Feature.Class"])
    else:
        norm_count_df, norm_count_avg_df, sample_rep_dict, classes = None, None, None, None
        scatters = ['sample_avg_scatter', 'sample_avg_scatter_by_class',
                    'sample_avg_scatter_by_deg', 'sample_avg_scatter_by_both']
        if any([dependency in args.plots for dependency in scatters]):
            raise Exception('The normalized count file is required for scatter plots, but it was not provided.')

    return norm_count_df, norm_count_avg_df, sample_rep_dict, classes

def load_deg_tables(args) -> pd.DataFrame:
    if args.deg_tables is not None:
        de_table = get_degs(args.deg_tables)
    else:
        de_table = None
        de_dependencies = ['sample_avg_scatter_by_deg', 'sample_avg_scatter_by_both']
        if any([dependency in args.plots for dependency in de_dependencies]):
            raise Exception('Differential expression tables are required for deg scatter plots, but none were provided.')

    return de_table

def main():
    """
    Main routine
    """
    args = get_args()

    # Create dataframes for raw/norm counts and deg tables if given
    raw_count_df, sample_rep_dict = load_raw_counts(args)
    norm_count_df, norm_count_avg_df, sample_rep_dict, classes = load_norm_counts(args)
    de_table = load_deg_tables(args)

    if all([counts is None for counts in [raw_count_df, norm_count_df]]) and 'replicate_scatter' in args.plots:
        raise Exception('Count files (norm or raw) are required for replicate scatter plots, '
                        'but neither were provided.')

    # generate plots requested
    for plot in args.plots:
        if plot == 'len_dist':
            if args.len_dist is not None:
                for infile in args.len_dist:
                    outfile = args.out_prefix + os.path.splitext(os.path.basename(infile))[0] + "_len_dist.pdf"
                    size_dist_plot(infile, outfile)
            else:
                print("5' nucleotide/length tables are required for length distribution plots, but none were provided.")

        elif plot == 'class_charts':
            pass
            # try:
            #     for infile in data_dict['class_counts']:
            #         outprefix = args.out_prefix + os.path.splitext(os.path.basename(infile[0]))[0]
            #         class_plots(infile, outprefix)
            # except KeyError as ke:
            #     print('No class_counts files found in args, but class charts requested.')
            #     print('Please check that -d/--dtypes has class_counts.')
            #     raise ke

        elif plot == 'replicate_scatter':
            # Refactor. I believe this will lead to overwritten outputs if both raw and norm counts are provided
            if raw_count_df is not None:
                scatter_replicates(raw_count_df, args.out_prefix, sample_rep_dict)
            
            if norm_count_df is not None:
                scatter_replicates(norm_count_df, args.out_prefix, sample_rep_dict)
        
        elif plot == 'sample_avg_scatter':
            if norm_count_df is not None:
                scatter_samples(norm_count_avg_df, args.out_prefix + "_avg")
        elif plot == 'sample_avg_scatter_by_class':
            if (norm_count_df is not None):
                scatter_samples(norm_count_avg_df, args.out_prefix + "_avg_by_class", classes=classes)
        elif plot == 'sample_avg_scatter_by_deg':
            if (norm_count_df is not None) and (de_table is not None):
                scatter_samples(norm_count_avg_df, args.out_prefix + "_avg_by_deg", degs=de_table)
        elif plot == 'sample_avg_scatter_by_both':
            if (norm_count_df is not None) and (de_table is not None):
                scatter_samples(norm_count_avg_df, args.out_prefix + "_avg_by_both", classes=classes, degs=de_table)
        else:
            print('Plot type %s not recognized, please check the -p/--plot arguments' % plot)


if __name__ == '__main__':
    main()

