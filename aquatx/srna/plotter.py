""" Plotting script for the sRNA AQuATx workflow.

This script produces basic static plots for publications as part of the AQuATx
workflow. It uses a default style and produces multiple PDFs. This script is
intended to be used immediately following the aquatx-count/counter.py step of
the workflow. It creates a specific set of plots through the mode argument.
"""

import argparse
import os.path
import itertools
import re
from collections import defaultdict
from typing import Optional, Dict, Union

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from aquatx.srna.plotterlib import plotterlib as lib
from aquatx.srna.util import report_execution_time

mpl.use("PDF")


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


def len_dist_plots(files_list, out_prefix, **kwargs):
    """Create a PDF of size and 5'nt distribution plot for a sample.

    Args:
        files_list: A list of files containing size + 5p-nt counts per library
        out_prefix: The prefix to use when naming output PDF plots
        kwargs: Keyword arguments to pass to matplotlib rc or plot 
    """
    for size_file in files_list:
        # Read the size_dist file
        size_dist = pd.read_csv(size_file, index_col=0)

        # Create the plot
        nt_colors = ['#F78E2D', '#CBDC3F', '#4D8AC8', '#E06EAA'] # Orange, Yellow-green, Blue, Pink
        ax = aqplt.size_dist_bar(size_dist, color=nt_colors, **kwargs)

        # Save the plot
        pdf_name = '_'.join([out_prefix, os.path.splitext(os.path.basename(size_file))[0], "len_dist.pdf"])
        plt.savefig(pdf_name, bbox_inches='tight')


def class_plots(class_counts, out_prefix, **kwargs):
    """Create a PDF of the proportion of counts assigned to
    a feature in a particular class as a pie and bar chart.

    Args:
        class_counts: A dataframe containing class counts per library
        out_prefix: The prefix to use when naming output PDF plots
        kwargs: Keyword arguments to pass to matplotlib rc or plot 
    """
    for library in class_counts:
        fig, ax = aqplt.class_pie_barh(class_counts[library], **kwargs)

        # Save the plot
        pdf_name = f"{out_prefix}_{library}_class_chart.pdf"
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


@report_execution_time("Scatter replicates")
def scatter_replicates(count_df, output_prefix, output_postfix, samples, norm=False):
    """Creates PDFs of all pairwise comparison scatter plots from a count table.
    
    Args:
        count_df: A dataframe of counts per features
        output_prefix: A string to use as a prefix for saving files
        output_postfix: A string to use as postfix for differentiating between raw and norm counts
        samples: A dictionary containing sample names and their associated "sample_rep_N" replicates
        norm: Boolean indicating if data was normalized
    """

    for samp, reps in samples.items():
        for pair in get_pairs(reps):
            rscat = aqplt.scatter_simple(count_df.loc[:,pair[0]], count_df.loc[:,pair[1]], 
                                        color='#888888', marker='s', alpha=0.5, s=50, 
                                        edgecolors='none', log_norm=True)
            rscat.set_title(samp if not norm else samp + " (normalized)")
            rep1, rep2 = pair[0].split('_rep_')[1], pair[1].split('_rep_')[1]
            rscat.set_xlabel('Replicate ' + rep1)
            rscat.set_ylabel('Replicate ' + rep2)
            pdf_name = '_'.join([output_prefix, samp, 'replicates', rep1, rep2,
                                 output_postfix, 'scatter.pdf'])
            rscat.figure.savefig(pdf_name, bbox_inches='tight')


def load_deg_tables(comparisons):
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
        samples[degfile] = name_split[-5] + "_vs_" + name_split[-3]

        if count == 0:
            de_table = pd.DataFrame(index=degs.index, columns=comparisons)

        de_table.loc[:, degfile] = degs.loc[:, 'padj']
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
        uniq_classes = list(pd.unique(classes))
        
        if not unknown:
            uniq_classes.remove('unknown')
        
        for pair in samples:
            grp_args = []
            try:
                samp_comp = '%s_vs_%s' % (pair[0], pair[1])
                deg_list = list(degs.index[degs[samp_comp] < pval])
            except KeyError:
                try: 
                    samp_comp = '%s_vs_%s' % (pair[1], pair[0])
                    deg_list = list(degs.index[degs[samp_comp] < pval])
                except KeyError as ke:
                    print('Sample names in count data frame do not match sample comparisons')
                    print('in DEG table. Make sure formatting is correct in your tables for')
                    print('AQuATx. Error occurred with %s, %s in count table.' % (pair[0], pair[1]))
                    raise ke
            
            class_degs = classes.loc[deg_list]

            for cls in uniq_classes:
                grp_args.append(list(class_degs.index[class_degs == cls]))

            labels = ['p > %6.2f' % pval] + uniq_classes
            sscat = aqplt.scatter_grouped(count_df.loc[:,pair[0]], count_df.loc[:,pair[1]],
                                          *grp_args, log_norm=True, labels=labels) 
            sscat.set_title('%s vs %s' % (pair[0], pair[1]))
            sscat.set_xlabel(pair[0])
            sscat.set_ylabel(pair[1])
            pdf_name = '_'.join([output_prefix, pair[0], 'vs', pair[1], 'scatter_by_deg_class.pdf'])
            sscat.figure.savefig(pdf_name, bbox_inches='tight')

    elif classes is not None:
        uniq_classes = list(pd.unique(classes))
        
        if not unknown:
            uniq_classes.remove('unknown')
        
        grp_args = []
        for cls in uniq_classes:
            grp_args.append(list(classes.index[classes == cls]))
    
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


def get_r_safename(name):
    leading_char = lambda x: re.sub(r"^(?=[^a-zA-Z.]+|\.\d)", "X", x)
    special_char = lambda x: re.sub(r"[^a-zA-Z0-9_.]", ".", x)
    return special_char(leading_char(name))


def load_raw_counts(raw_counts_file) -> Optional[pd.DataFrame]:
    if raw_counts_file is None: return None
    raw_count_df = pd.read_csv(raw_counts_file, index_col=0)
    raw_count_df.rename(get_r_safename, axis="columns")

    return raw_count_df


def load_norm_counts(norm_counts_file) -> Optional[pd.DataFrame]:
    if norm_counts_file is None: return None
    return pd.read_csv(norm_counts_file, index_col=0)


def get_flat_classes(count_df) -> pd.Series:
    # Features with multiple associated classes must have these classes flattened
    count_df["Feature.Class"] \
        .apply(lambda x: [cls.strip() for cls in x.split(',')]) \
        .explode()


def get_class_counts(counts_df) -> pd.DataFrame:
    # 1. Group and sum by Feature Class. This makes the next step less expensive.
    grouped = counts_df.groupby("Feature.Class").sum()

    # 2. Groups may also contain list (string) entries which need to be flattened for proper attribution.
    grouped_flat_index = [row._replace(Index=split.strip()) if ',' in row[0] else row
                          for row in grouped.itertuples()
                          for split in row[0].split(',')]

    # 3. Finally, group and sum again now that class lists have been flattened.
    class_counts = pd.DataFrame(grouped_flat_index, columns=counts_df.columns).groupby("Feature.Class").sum()

    return class_counts


@report_execution_time("Input validation")
def validate_inputs(args):
    dependencies_satisfied_for = {
        'len_dist': all([args.len_dist]),
        'class_charts': all([args.norm_counts]),
        'replicate_scatter': any([args.norm_counts, args.raw_counts]),
        'sample_avg_scatter': all([args.norm_counts]),
        'sample_avg_scatter_by_class': all([args.norm_counts]),
        'sample_avg_scatter_by_deg': all([args.norm_counts, args.deg_tables]),
        'sample_avg_scatter_by_both': all([args.norm_counts, args.deg_tables]),
    }

    nc_file = "normalized count file (-nc/--norm-counts)"
    rc_file = "raw count file (-rc/--raw-counts)"
    deg_files = "differential expression tables (-deg/--deg-tables)"

    error_message = {
        'len_dist': "5' nucleotide/length tables (-len/--len-dist) are required for ",
        'class_charts': f"{nc_file} is required for ",
        'replicate_scatter': f"A {nc_file} or a {rc_file} is required for ",
        'sample_avg_scatter': f"A {nc_file} is required for ",
        'sample_avg_scatter_by_class': f"A {nc_file} is required for ",
        'sample_avg_scatter_by_deg': f"A {nc_file} and {deg_files} are required for ",
        'sample_avg_scatter_by_both': f"A {nc_file} and {deg_files} are required for ",
    }

    unsatisfied = []
    for plot_name in args.plots:
        if not dependencies_satisfied_for[plot_name]:
            print(error_message[plot_name] + plot_name)
            unsatisfied.append(plot_name)

    for bad_request in unsatisfied:
        args.plots.remove(bad_request)


@report_execution_time("Setup")
def setup(args):
    required_inputs = {
        'len_dist': [],
        'class_charts': ["norm_count_df", "class_counts"],
        'replicate_scatter': ["raw_count_df", "norm_count_df", "sample_rep_dict"],
        'sample_avg_scatter': ["norm_count_df", "norm_count_avg_df"],
        'sample_avg_scatter_by_class': ["norm_count_df", "norm_count_avg_df", "feat_classes"],
        'sample_avg_scatter_by_deg': ["norm_count_df", "norm_count_avg_df", "de_table"],
        'sample_avg_scatter_by_both': ["norm_count_df", "norm_count_avg_df", "feat_classes", "de_table"],
    }

    relevant_vars: Dict[str, Union[pd.DataFrame, pd.Series, dict, None]] = {}
    input_getters = {
        'raw_count_df': lambda: load_raw_counts(args.raw_counts),
        'norm_count_df': lambda: load_norm_counts(args.norm_counts),
        'de_table': lambda: load_deg_tables(args.deg_tables),
        'sample_rep_dict': lambda: get_sample_rep_dict(relevant_vars["norm_count_df"]),
        'norm_count_avg_df': lambda: get_sample_averages(relevant_vars["norm_count_df"], relevant_vars["sample_rep_dict"]),
        'feat_classes': lambda: get_flat_classes(relevant_vars["norm_count_df"]),
        'class_counts': lambda: get_class_counts(relevant_vars["norm_count_df"])
    }

    for plot in args.plots:
        for req in required_inputs[plot]:
            if req is not None and req not in relevant_vars:
                relevant_vars[req] = input_getters[req]()

    return relevant_vars


@report_execution_time("Main routine")
def main():
    """
    Main routine
    """
    args = get_args()
    validate_inputs(args)
    inputs = setup(args)

    global aqplt
    aqplt = lib()

    # generate plots requested
    for plot in args.plots:
        if plot == 'len_dist':
            len_dist_plots(args.len_dist, args.out_prefix)
        elif plot == 'class_charts':
            class_plots(inputs["class_counts"], args.out_prefix)
        elif plot == 'replicate_scatter':
            if inputs["raw_count_df"] is not None:
                scatter_replicates(inputs["raw_count_df"], args.out_prefix, "raw_count", inputs["sample_rep_dict"])
            
            if inputs["norm_count_df"] is not None:
                scatter_replicates(inputs["norm_count_df"], args.out_prefix, "norm_count", inputs["sample_rep_dict"], norm=True)
        elif plot == 'sample_avg_scatter':
            scatter_samples(inputs["norm_count_avg_df"], args.out_prefix)
        elif plot == 'sample_avg_scatter_by_class':
            scatter_samples(inputs["norm_count_avg_df"], args.out_prefix, classes=inputs["feat_classes"])
        elif plot == 'sample_avg_scatter_by_deg':
            scatter_samples(inputs["norm_count_avg_df"], args.out_prefix, degs=inputs["de_table"])
        elif plot == 'sample_avg_scatter_by_both':
            scatter_samples(inputs["norm_count_avg_df"], args.out_prefix, classes=inputs["feat_classes"], degs=inputs["de_table"])
        else:
            print('Plot type %s not recognized, please check the -p/--plot arguments' % plot)


if __name__ == '__main__':
    main()

