""" Plotting script for the sRNA tinyRNA workflow.

This script produces basic static plots for publications as part of the tinyRNA
workflow. It uses a default style and produces multiple PDFs, but you may provide
your own styles sheet.
"""

import multiprocessing as mp
import pandas as pd
import itertools
import argparse
import os.path
import re

from collections import defaultdict
from typing import Optional, Dict, Union, Tuple
from pkg_resources import resource_filename

from tiny.rna.configuration import timestamp_format
from tiny.rna.plotterlib import plotterlib as lib
from tiny.rna.util import report_execution_time, make_filename


def get_args():
    """Get input arguments from the user/command line."""

    parser = argparse.ArgumentParser()
    counter_files = parser.add_argument_group("Input files produced by Counter")
    diffexp_files = parser.add_argument_group("Input files produced by DGE")

    # Single file inputs
    diffexp_files.add_argument('-nc', '--norm-counts', metavar='NORM_COUNTS',
                               help='The ...norm_counts.csv file.')

    # Multi-file inputs
    diffexp_files.add_argument('-dge', '--dge-tables', metavar='DGE_COMPARISONS', nargs='+',
                               help='The ...cond1...cond2...deseq.csv files, separated by a space.')
    counter_files.add_argument('-len', '--len-dist', metavar='5P_LEN_DISTS', nargs='+',
                               help='The ...nt_len_dist.csv files, separated by a space.')

    # Outputs options
    parser.add_argument('-o', '--out-prefix', metavar='OUTPREFIX',
                        help='Optional prefix to use for output PDF files.')
    parser.add_argument('-pv', '--p-value', metavar='VAL', default=0.05, type=float,
                        help='Optional p-value to use in DGE scatter plots.')
    parser.add_argument('-s', '--style-sheet', metavar='MPLSTYLE',
                        default=resource_filename('tiny', 'templates/tinyrna-light.mplstyle'),
                        help='Optional matplotlib style sheet to use for plots.')
    parser.add_argument('-v', '--vector-scatter', action='store_true',
                        help='Produce scatter plots with vectorized points (slower).\n'
                        'Note: only the points on scatter plots will be raster if '
                        'this option is not provided.')
    parser.add_argument('-p', '--plots', metavar='PLOTS', required=True, nargs='+',
                        help='List of plots to create. Options: \n'
                             'len_dist: A stacked barplot showing size & 5p-nt distribution,\n'
                             'class_charts: A pie and barchart showing proportions of counts per class,\n'
                             'replicate_scatter: A scatter plot comparing replicates for all count files given,\n'
                             'sample_avg_scatter_by_dge: A scatter plot comparing all sample groups,'
                                'with significantly different (padj<0.05) genes highlighted.\n'
                             'sample_avg_scatter_by_dge_class: A scatter plot comparing all sample groups,'
                                'with classes and significantly different genes highlighted')

    return parser.parse_args()


def get_pairs(samples):
    """Get pairs of comparisons from a list of samples to compare.
    
    Args:
        samples: A list of samples to compare
    
    Returns:
        pair: An generator of all possible pairs
    """
    
    for pair in itertools.combinations(samples, 2):
        yield pair


def len_dist_plots(files_list: list, out_prefix:str, **kwargs):
    """Create a PDF of size and 5'nt distribution plot for a sample.

    Args:
        files_list: A list of files containing size + 5p-nt counts per library
        out_prefix: The prefix to use when naming output PDF plots
        kwargs: Additional keyword arguments to pass to pandas.DataFrame.plot()
    """

    for size_file in files_list:
        # Read the size_dist file
        size_dist = pd.read_csv(size_file, index_col=0)

        # Create the plot
        plot = aqplt.len_dist_bar(size_dist, **kwargs)

        # Parse the "sample_rep_N" string from the input filename to avoid duplicate out_prefix's in the basename
        basename = os.path.splitext(os.path.basename(size_file))[0]
        date_prefix_pos = re.search(timestamp_format, basename).span()

        if date_prefix_pos is not None:
            # File is a pipeline product
            begin = date_prefix_pos[1] + 1
            end = basename.rfind("_nt_len_dist")
            condition_and_rep = basename[begin:end]
        else:
            # File does not appear to have been produced by the pipeline
            condition_and_rep = basename

        pdf_name = make_filename([out_prefix, condition_and_rep, "len_dist"], ext='.pdf')
        plot.figure.savefig(pdf_name)


def class_plots(class_counts: pd.DataFrame, out_prefix:str, **kwargs):
    """Create a PDF of the proportion of counts assigned to
    a feature in a particular class as a pie and bar chart.

    Args:
        class_counts: A dataframe containing class counts per library
        out_prefix: The prefix to use when naming output PDF plots
        kwargs: Additional keyword arguments to pass to pandas.DataFrame.plot()
    """

    for library in class_counts:
        fig = aqplt.class_pie_barh(class_counts[library], **kwargs)

        # Save the plot
        pdf_name = make_filename([out_prefix, library, 'class_chart'], ext='.pdf')
        fig.savefig(pdf_name)


def get_sample_rep_dict(df: pd.DataFrame) -> dict:
    """Group sample_rep_N format strings by sample
    
    Args:
        df: The dataframe with columns using sample_rep_N format and features for rows

    Returns:
        sample_dict: A dictionary containing all column names corresponding to unique samples.
    """
    
    sample_dict = defaultdict(list)

    for col in df.columns:
        if col == "Feature Class": continue
        sample = col.split("_rep_")[0]
        sample_dict[sample].append(col)

    return sample_dict


def get_sample_averages(df: pd.DataFrame, samples:dict) -> pd.DataFrame:
    """Average counts across replicates on a per-sample basis
    
    Args:
        df: A dataframe with columns using sample_rep_N format and features for rows
        samples: A dictionary containing sample names and their associated "sample_rep_N" replicates

    Returns:
        new_df: The new dataframe with averaged counts per unique sample.
    """
    
    new_df = pd.DataFrame(index=df.index, columns=samples.keys())

    for key, val in samples.items():
        new_df[key] = df.loc[:, val].mean(axis=1)
    
    return new_df


def scatter_replicates(count_df: pd.DataFrame, output_prefix: str, samples: dict, view_lims: Tuple[float, float] = None):
    """Creates pairwise scatter plots comparing replicates' counts and saves the plot as a PDF
    
    Args:
        count_df: A dataframe of counts per feature
        output_prefix: A string to use as a prefix for saving files
        samples: A dictionary containing sample names and their associated "sample_rep_N" replicate names
        view_lims: Optional plot view limits as tuple(min, max)
    """

    for samp, reps in samples.items():
        for pair in get_pairs(reps):
            rscat = aqplt.scatter_simple(count_df.loc[:,pair[0]], count_df.loc[:,pair[1]],
                                         color='#B3B3B3', alpha=0.5, log_norm=True, rasterized=RASTER)
            aqplt.set_square_scatter_view_lims(rscat, view_lims)
            aqplt.set_scatter_ticks(rscat)
            rscat.set_title(samp)
            rep1, rep2 = pair[0].split('_rep_')[1], pair[1].split('_rep_')[1]
            rscat.set_xlabel('Log$_{2}$ normalized reads in replicate ' + rep1)
            rscat.set_ylabel('Log$_{2}$ normalized reads in replicate ' + rep2)
            pdf_name = make_filename([output_prefix, samp, 'replicates', rep1, rep2, 'scatter'], ext='.pdf')
            rscat.figure.savefig(pdf_name)


def load_dge_tables(comparisons: list) -> pd.DataFrame:
    """Creates a new dataframe containing all features and padj values for each comparison
    from a list of DGE tables.

    Args:
        comparisons: DGE tables (files) produced by DESeq2

    Returns:
        de_table: A single table of padj values per feature/comparison
    """

    de_table = pd.DataFrame()

    for dgefile in comparisons:
        # Negative indexes are used since prefixes are allowed to contain underscores
        name_split = os.path.basename(dgefile).split('_')
        comparison_name = name_split[-5] + "_vs_" + name_split[-3]

        de_table[comparison_name] = pd.read_csv(dgefile, index_col=0)['padj']
    
    return de_table


def scatter_dges(count_df, dges, output_prefix, viewLims, classes=None, show_unknown=False, pval=0.05):
    """Creates PDFs of all pairwise comparison scatter plots from a count table.
    Can highlight classes and/or differentially expressed genes as different colors.

    Args:
        count_df: A dataframe of counts per feature
        dges: A dataframe of differential gene table output to highlight
        output_prefix: A string to use as a prefix for saving files
        classes: A dataframe containing class(es) per feature
        show_unknown: If true, class "unknown" will be included if highlighting by classes
    """

    if classes is not None:
        uniq_classes = list(pd.unique(classes))

        if not show_unknown and 'unknown' in uniq_classes:
            uniq_classes.remove('unknown')

        for pair in dges:
            p1, p2 = pair.split("_vs_")
            dge_list = list(dges.index[dges[pair] < pval])
            class_dges = classes.loc[dge_list]

            grp_args = []
            for cls in uniq_classes:
                grp_args.append(list(class_dges.index[class_dges == cls]))

            labels = ['p ≥ %g' % pval] + uniq_classes
            sscat = aqplt.scatter_grouped(count_df.loc[:,p1], count_df.loc[:,p2], viewLims, *grp_args,
                                          log_norm=True, labels=labels, rasterized=RASTER)
            sscat.set_title('%s vs %s' % (p1, p2))
            sscat.set_xlabel("Log$_{2}$ normalized reads in " + p1)
            sscat.set_ylabel("Log$_{2}$ normalized reads in " + p2)
            pdf_name = make_filename([output_prefix, pair, 'scatter_by_dge_class'], ext='.pdf')
            sscat.figure.savefig(pdf_name)

    else:
        for pair in dges:
            grp_args = list(dges.index[dges[pair] < pval])
            p1, p2 = pair.split("_vs_")

            labels = ['p ≥ %g' % pval, 'p < %g' % pval]
            sscat = aqplt.scatter_grouped(count_df.loc[:,p1], count_df.loc[:,p2], viewLims, grp_args,
                                          log_norm=True, labels=labels, alpha=0.5, rasterized=RASTER)
            sscat.set_title('%s vs %s' % (p1, p2))
            sscat.set_xlabel("Log$_{2}$ normalized reads in " + p1)
            sscat.set_ylabel("Log$_{2}$ normalized reads in " + p2)
            pdf_name = make_filename([output_prefix, pair, 'scatter_by_dge'], ext='.pdf')
            sscat.figure.savefig(pdf_name)


def load_norm_counts(norm_counts_file: str) -> Optional[pd.DataFrame]:
    """Loads a norm_counts CSV as a DataFrame
    Args:
        norm_counts_file: The norm counts CSV produced by DESeq2
    Returns:
        The norm counts DataFrame
    """

    if norm_counts_file is None: return None
    return pd.read_csv(norm_counts_file, index_col=0)


def get_flat_classes(counts_df: pd.DataFrame) -> pd.Series:
    """Features with multiple associated classes must have these classes flattened
    In these cases, the feature will be listed multiple times in the index with one associated class per row.

    Args:
        counts_df: A counts DataFrame
    Returns:
        A Series with an index of features, and entries for each associated class
    """

    return counts_df["Feature Class"] \
        .apply(lambda x: [cls.strip() for cls in x.split(',')]) \
        .explode()


def get_class_counts(counts_df: pd.DataFrame) -> pd.DataFrame:
    """Calculates class counts from a counts DataFrame

    If there are multiple classes associated with a feature, each class will receive
    the feature's full count without normalizing by the number of associated classes.

    Args:
        counts_df: A DataFrame containing features and their associated classes and counts
    Returns:
        class_counts: A DataFrame with an index of classes and count entries, per comparison
    """

    # 1. Group by Feature Class. Class "lists" are strings. Normalize their counts by the number of commas.
    grouped = counts_df.groupby("Feature Class").apply(lambda grp: grp.sum() / (grp.name.count(',') + 1))

    # 2. Convert class "list" strings to true lists since we no longer need a hashable index
    grouped.index = grouped.index.map(lambda x: [cls.strip() for cls in x.split(',')])

    # 3. Finally, flatten class lists and add the normalized counts to their groups
    return grouped.reset_index().explode("Feature Class").groupby("Feature Class").sum()


def validate_inputs(args: argparse.Namespace) -> None:
    """Determines if the necessary input files have been provided for the requested plots
    This is necessary because we allow users to run the tool with only the files necessary
    for the plot types they've requested, rather than all possible inputs.

    Args:
        args: Command line arguments parsed by argparse
    """

    dependencies_satisfied_for = {
        'len_dist': all([args.len_dist]),
        'class_charts': all([args.norm_counts]),
        'replicate_scatter': all([args.norm_counts]),
        'sample_avg_scatter_by_dge': all([args.norm_counts, args.dge_tables]),
        'sample_avg_scatter_by_dge_class': all([args.norm_counts, args.dge_tables]),
    }

    nc_file = "normalized count file (-nc/--norm-counts)"
    dge_files = "differential expression tables (-dge/--dge-tables)"

    error_message = {
        'len_dist': "5' nucleotide/length tables (-len/--len-dist) are required for ",
        'class_charts': f"A {nc_file} is required for ",
        'replicate_scatter': f"A {nc_file} is required for ",
        'sample_avg_scatter_by_dge': f"A {nc_file} and {dge_files} are required for ",
        'sample_avg_scatter_by_dge_class': f"A {nc_file} and {dge_files} are required for ",
    }

    for plot_name in args.plots.copy():
        if plot_name not in dependencies_satisfied_for:
            print(f"Skipping unrecognized plot type: {plot_name}")
            args.plots.remove(plot_name)
            continue
        if not dependencies_satisfied_for[plot_name]:
            print(error_message[plot_name] + plot_name)
            args.plots.remove(plot_name)


def setup(args: argparse.Namespace) -> dict:
    """Dynamically loads input files based on the plots requested
    We don't want to waste time processing files that are unnecessary for the request

    Args:
        args: Command line arguments parsed by argparse
    Returns:
        relevant_vars: A dictionary of input names and their processed data
    """

    required_inputs = {
        'len_dist': [],
        'class_charts': ["norm_count_df", "class_counts"],
        'replicate_scatter': ["norm_count_df", "sample_rep_dict", "norm_view_lims"],
        'sample_avg_scatter_by_dge':
            ["norm_count_df", "sample_rep_dict", "norm_count_avg_df", "de_table", "avg_view_lims"],
        'sample_avg_scatter_by_dge_class':
            ["norm_count_df", "sample_rep_dict", "norm_count_avg_df", "feat_classes", "de_table", "avg_view_lims"],
    }

    relevant_vars: Dict[str, Union[pd.DataFrame, pd.Series, dict, None]] = {}
    input_getters = {
        'norm_count_df': lambda: load_norm_counts(args.norm_counts),
        'de_table': lambda: load_dge_tables(args.dge_tables),
        'sample_rep_dict': lambda: get_sample_rep_dict(relevant_vars["norm_count_df"]),
        'norm_count_avg_df': lambda: get_sample_averages(relevant_vars["norm_count_df"], relevant_vars["sample_rep_dict"]),
        'feat_classes': lambda: get_flat_classes(relevant_vars["norm_count_df"]),
        'class_counts': lambda: get_class_counts(relevant_vars["norm_count_df"]),
        'avg_view_lims': lambda: aqplt.get_scatter_view_lims(relevant_vars["norm_count_avg_df"]),
        'norm_view_lims': lambda: aqplt.get_scatter_view_lims(relevant_vars["norm_count_df"].select_dtypes(['number']))
    }

    for plot in args.plots:
        for req in required_inputs[plot]:
            if req is not None and req not in relevant_vars:
                relevant_vars[req] = input_getters[req]()

    return relevant_vars


@report_execution_time("Plotter runtime")
def main():
    """
    Main routine
    """
    args = get_args()
    validate_inputs(args)

    global aqplt, RASTER
    aqplt = lib(args.style_sheet)
    RASTER = not args.vector_scatter
    inputs = setup(args)

    # Assemble work units for multiprocessing pool
    itinerary = []

    # generate plots requested
    for plot in args.plots:
        if plot == 'len_dist':
            func = len_dist_plots
            arg, kwd = (args.len_dist, args.out_prefix), {}
        elif plot == 'class_charts':
            func = class_plots
            arg, kwd = (inputs["class_counts"], args.out_prefix), {}
        elif plot == 'replicate_scatter':
            func = scatter_replicates
            arg = (inputs["norm_count_df"], args.out_prefix, inputs["sample_rep_dict"], inputs["norm_view_lims"])
            kwd = {}
        elif plot == 'sample_avg_scatter_by_dge':
            func = scatter_dges
            arg = (inputs["norm_count_avg_df"], inputs["de_table"], args.out_prefix, inputs["avg_view_lims"])
            kwd = {"pval": args.p_value}
        elif plot == 'sample_avg_scatter_by_dge_class':
            func = scatter_dges
            arg = (inputs["norm_count_avg_df"], inputs["de_table"], args.out_prefix, inputs["avg_view_lims"])
            kwd = {"classes": inputs["feat_classes"], "pval": args.p_value}
        else:
            print('Plot type %s not recognized, please check the -p/--plot arguments' % plot)
            continue

        itinerary.append((func, arg, kwd))

    if len(itinerary) > 1:
        with mp.Pool(len(itinerary)) as pool:
            results = []
            for task, args, kwds in itinerary:
                results.append(pool.apply_async(task, args, kwds))
            for result in results:
                result.wait()
    elif len(itinerary) == 1:
        # Don't use multiprocessing if only one plot type requested
        task = itinerary.pop()
        task[0](*task[1], **task[2])

if __name__ == '__main__':
    main()

