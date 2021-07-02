""" Plotting script for the sRNA AQuATx workflow.

This script produces basic static plots for publications as part of the AQuATx
workflow. It uses a default style and produces multiple PDFs. This script is
intended to be used immediately following the aquatx-count/counter.py step of
the workflow. It creates a specific set of plots through the mode argument.
"""

import argparse
import os.path
import itertools
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import aquatx.srna.plotterlib as aqplt
from HTSeq import GFF_Reader

def get_args():
    """Get input arguments from the user/command line."""

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-files', metavar='DATAFILE', required=True, nargs='+',
                        help='input files with data from final merged tables')
    parser.add_argument('-o', '--out-prefix', metavar='OUTFILE', default='',
                        help='Optional prefix to use for output PDF files.')
    parser.add_argument('-d', '--data-types', metavar='DATATYPE', required=True, nargs='+',
                        help='List of data types corresponding to input files. Options: '\
                        '\nraw_counts: the raw counts per feature data,'\
                        '\nnorm_counts: the normalized counts per feature data,'\
                        '\ndegs: the table of differential gene expression analysis including p-values,'\
                        '\nlen_dist: the size and 5p nucleotide count matrix,'\
                        '\nclass_counts: the table of raw counts per class')
    parser.add_argument('-r', '--references', metavar='GFF3_FILES', nargs='+',
                        help='The reference annotation files containing class information to'\
                        'highlight on plots.')
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
        new_df.loc[:, key] = df[val].mean(axis=1)
    
    return new_df

def get_annotations(count_df, *args):
    """Get feature + class information from reference annotation files

    Args:
        args: A list of all files to process

    Returns:
        classes: A dataframe containing features in index + class information.
    """
    classes = pd.DataFrame(index=count_df.index, columns=['class'])

    for afile in args:
        feat_gff = GFF_Reader(afile)
        for feat in feat_gff:
            if pd.isnull(classes.loc[feat.attr["ID"]]).any():
                classes.loc[feat.attr["ID"]] = feat.type
            else:
                classes.loc[feat.attr["ID"]] = 'ambiguous'

    return classes

def scatter_replicates(count_df, output_prefix, norm=False):
    """Creates PDFs of all pairwise comparison scatter plots from a count table.
    
    Args:
        count_df: A dataframe of counts per features
        output_prefix: A string to use as a prefix for saving files
        norm: Boolean indicating if data was normalized
    """

    samples = get_replicates(count_df)

    for samp, reps in samples.items():
        for pair in get_pairs(reps):
            rscat = aqplt.scatter_simple(count_df.loc[:,pair[0]], count_df.loc[:,pair[1]], 
                                        color='#888888', marker='s', alpha=0.5, s=50, 
                                        edgecolors='none', log_norm=True)
            rscat.set_title(samp)
            rscat.set_xlabel('Replicate ' + pair[0].split('_replicate_')[1])
            rscat.set_ylabel('Replicate ' + pair[1].split('_replicate_')[1])
            pdf_name = '_'.join([output_prefix, samp, 'replicates', pair[0].split('_replicate_')[1],
                                pair[1].split('_replicate_')[1], 'scatter.pdf'])
            rscat.figure.savefig(pdf_name, bbox_inches='tight')

def get_degs(*args):
    """Create a new dataframe containing all features and pvalues for each comparison
    from a set of DEG tables. 

    Args:
        args: Files (DEG tables) to process

    Returns:
        de_table: A single table of p-values per feature/comparison
    """
    count = 0
    samples = {}
    
    for degfile in args:
        degs = pd.read_csv(degfile, index_col=0)
        samples[degfile] = os.path.basename(degfile).split('.')[:-1][0].split('_cond1_')[1].split('_deseq_')[0].replace('_cond2_','_vs_')

        if count == 0:
            de_table = pd.DataFrame(index=degs.index, columns=args)
        
        de_table.loc[:, degfile] = degs.loc[:, 'padj']
        count += 1
        
    de_table.rename(samples, axis=1, inplace=True)
    
    return de_table

def scatter_samples(count_df, output_prefix, classes=None, degs=None, ambig=False, pval=0.05):
    """Creates PDFs of all pairwise comparison scatter plots from a count table.
    Can highlight classes and/or differentially expressed genes as different colors.

    Args:
        count_df: A dataframe of counts per feature
        output_prefix: A string to use as a prefix for saving files
        classes: A dataframe containing class per feature
        degs: A dataframe of differential gene table output to highlight
        ambig: Boolean indicating if ambiguous classes should be highlighted
    """
    samples = get_pairs(list(count_df.columns))

    if (classes is not None) and (degs is not None):
        uniq_classes = list(pd.unique(classes.loc[:,'class'].dropna()))
        
        if not ambig:
            uniq_classes.remove('ambiguous')
        
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
        uniq_classes = list(pd.unique(classes.loc[:, 'class'].dropna()))
        
        if not ambig:
            uniq_classes.remove('ambiguous')
        
        grp_args = []
        for cls in uniq_classes:
            grp_args.append(list(classes.loc[classes['class'] == cls].index))
    
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

def main():
    """ 
    Main routine
    """
    args = get_args()

    # Sort files & types
    data_dict=dict()
    for infile, dtype in zip(args.input_files, args.data_types):
        data_dict.setdefault(dtype, []).append(infile)
    
    # Create count dataframes if given
    try: 
        raw_count_df = pd.read_csv(data_dict['raw_counts'][0], index_col=0)
    except KeyError:
        raw_count_df = None
            
    try:
        norm_count_df = pd.read_csv(data_dict['norm_counts'][0], index_col=0)
        norm_count_avg_df = get_sample_averages(norm_count_df)
    
        # get class information if given
        try:
            classes = get_annotations(norm_count_df, args.references)
        except NameError:
            classes = None
            if bool(set('sample_avg_scatter_by_class', 'sample_avg_scatter_by_both') & set(args.plots)):
                raise Exception('No reference annotation file given, but class highlighting on scatter'\
                                'plots requested. Please check -r/--references was given.')

    except KeyError:
        norm_count_df = None
        scatters = ['sample_avg_scatter', 'sample_avg_scatter_by_class', 
                     'sample_avg_scatter_by_deg', 'sample_avg_scatter_by_both']
        if bool(set(scatters) & set(args.plots)):
            raise Exception('No normalized count files specified in -d/--dtypes, but sample scatter'\
                            'plots were requested. Please check your arguments for norm_counts.')
    
    if ('replicate_scatter' in args.plots) and (raw_count_df is None) and (norm_count_df is None):
        raise Exception('No count files specified in -d/--dtypes, but replicate scatter plots requested.'\
                        'Please check your arguments.')

    # create merged de table for highlighting
    try:
        de_table = get_degs(data_dict['degs'])

    except KeyError:
        de_table = None
        if bool(set(['sample_avg_scatter_by_deg', 'sample_avg_scatter_by_both']) & set(args.plots)):
            raise Exception('No differential expression tables specified in -d/--dtypes, but deg scatter'\
                            'plots were requested. Please check your arguments for degs.')
    
    # generate plots requested
    for plot in args.plots:
        if plot is 'len_dist':
            try:
                for infile in data_dict['len_dist']:
                    outprefix = args.out_prefix + ''.join(os.path.basename(infile[0]).split('.')[:-1])
                    size_dist_plot(infile, outprefix)
            except KeyError as ke:
                print('No len_dist files found in args, but length distribution plots requested.')
                print('Please check that -d/--dtypes has len_dist.')
                raise ke

        elif plot is 'class_charts':
            try:
                for infile in data_dict['class_counts']:
                    outprefix = args.out_prefix + ''.join(os.path.basename(infile[0]).split('.')[:-1])
                    class_plots(infile, outprefix)
            except KeyError as ke:
                print('No class_counts files found in args, but class charts requested.')
                print('Please check that -d/--dtypes has class_counts.')
                raise ke

        elif plot is 'replicate_scatter':
            if raw_count_df is not None:
                scatter_replicates(raw_count_df, args.out_prefix)
            
            if norm_count_df is not None:
                scatter_replicates(norm_count_df, args.out_prefix)
        
        elif plot is 'sample_avg_scatter':
            if norm_count_df is not None:
                scatter_samples(norm_count_avg_df, args.out_prefix)
        elif plot is 'sample_avg_scatter_by_class':
            if (norm_count_df is not None) and (classes is not None):
                scatter_samples(norm_count_avg_df, args.out_prefix, classes=classes)
        elif plot is 'sample_avg_scatter_by_deg':
            if (norm_count_df is not None) and (de_table is not None):
                scatter_samples(norm_count_avg_df, args.out_prefix, degs=de_table)
        elif plot is 'sample_avg_scatter_by_both':
            if (norm_count_df is not None) and (de_table is not None) and  (classes is not None):
                scatter_samples(norm_count_avg_df, args.out_prefix, classes=classes, degs=de_table)
        else:
            print('Plot type %s not recognized, please check the -p/--plot arguments' % plot)


if __name__ == '__main__':
    main()

