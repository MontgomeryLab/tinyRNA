"""This script produces basic static plots for publication as part of the tinyRNA workflow.

Input file requirements vary by plot type and you are free to supply only the files necessary
for your plot selections. If you are sourcing all of your input files from the same run
directory, you may find it easier to instead run `tiny replot` within that run directory.
"""
import multiprocessing as mp
import pandas as pd
import itertools
import traceback
import argparse
import os.path
import sys
import csv
import re

from collections import defaultdict
from typing import Optional, Dict, Union, Tuple
from pkg_resources import resource_filename

from tiny.rna.configuration import timestamp_format
from tiny.rna.plotterlib import plotterlib as lib
from tiny.rna.util import report_execution_time, make_filename, SmartFormatter


def get_args():
    """Get input arguments from the user/command line."""

    parser = argparse.ArgumentParser(description=__doc__, add_help=False, formatter_class=SmartFormatter)
    required_args = parser.add_argument_group("Required arguments")
    counter_files = parser.add_argument_group("Input files produced by Counter")
    diffexp_files = parser.add_argument_group("Input files produced by DGE")
    optional_args = parser.add_argument_group("Optional arguments")

    # Single file inputs
    counter_files.add_argument('-rc', '--raw-counts', metavar='RAW_COUNTS',
                               help='The ...feature_counts.csv file')
    diffexp_files.add_argument('-nc', '--norm-counts', metavar='NORM_COUNTS',
                               help='The ...norm_counts.csv file')
    counter_files.add_argument('-ss', '--summary-stats', metavar='STAT',
                               help='The ...summary_stats.csv file')

    # Multi-file inputs
    diffexp_files.add_argument('-dge', '--dge-tables', metavar='COMPARISON', nargs='+',
                               help='The ...cond1...cond2...deseq.csv files')
    counter_files.add_argument('-len', '--len-dist', metavar='5P_LEN', nargs='+',
                               help='The ...nt_len_dist.csv files')

    # Outputs options
    optional_args.add_argument('-h', '--help', action="help", help="show this help message and exit")
    optional_args.add_argument('-o', '--out-prefix', metavar='PREFIX',
                               help='Prefix to use for output filenames.')
    optional_args.add_argument('-pv', '--p-value', metavar='VALUE', default=0.05, type=float,
                               help='P-value to use in DGE scatter plots.')
    optional_args.add_argument('-s', '--style-sheet', metavar='MPLSTYLE',
                               default=resource_filename('tiny', 'templates/tinyrna-light.mplstyle'),
                               help='Optional matplotlib style sheet to use for plots.')
    optional_args.add_argument('-v', '--vector-scatter', action='store_true',
                               help='Produce scatter plots with vectorized points (slower).\n'
                               'Note: only the points on scatter plots will be raster if '
                               'this option is not provided.')

    # Required arguments
    required_args.add_argument('-p', '--plots', metavar='PLOT', required=True, nargs='+',
                               help="R|List of plots to create. Options: \n"
                               "• len_dist: A stacked barplot showing size & 5' nucleotide distribution.\n"
                               "• class_charts: A pie and barchart showing proportions of counts per class.\n"
                               "• replicate_scatter: A scatter plot comparing replicates for all count files given.\n"
                               "• sample_avg_scatter_by_dge: A scatter plot comparing all sample groups, with "
                               "significantly different genes highlighted. P-value can be set using --p-value. "
                               "Default: 0.05.\n"
                               "• sample_avg_scatter_by_dge_class: A scatter plot comparing all sample groups, with "
                               "classes and significantly different genes highlighted.")

    return parser.parse_args()


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


def class_charts(raw_class_counts: pd.DataFrame, mapped_reads: pd.Series, out_prefix: str, **kwargs):
    """Create a PDF of the proportion of raw assigned counts vs total mapped reads for each class

    Args:
        raw_class_counts: A dataframe containing RAW class counts per library
        mapped_reads: A Series containing the total number of mapped reads per library
        out_prefix: The prefix to use when naming output PDF plots
        kwargs: Additional keyword arguments to pass to pandas.DataFrame.plot()
    """

    for library in raw_class_counts:
        fig = aqplt.class_pie_barh(raw_class_counts[library], mapped_reads[library], **kwargs)

        # Save the plot
        pdf_name = make_filename([out_prefix, library, 'class_chart'], ext='.pdf')
        fig.savefig(pdf_name)


def get_mapped_reads(summary_stats_file: str) -> pd.Series:
    """Produces a Series of mapped reads per library for calculating proportions in class charts.

    Args:
        summary_stats_file: the summary stats csv produced by Counter

    Returns: a Series containing mapped reads per sample
    """

    return pd.read_csv(summary_stats_file, index_col=0).loc["Mapped Reads",:]


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
        for pair in itertools.combinations(reps, 2):
            rscat = aqplt.scatter_simple(count_df.loc[:,pair[0]], count_df.loc[:,pair[1]],
                                         color='#A1A1A1', alpha=0.5, log_norm=True, rasterized=RASTER)
            aqplt.set_square_scatter_view_lims(rscat, view_lims)
            aqplt.set_scatter_ticks(rscat)
            rscat.set_title(samp)
            rep1, rep2 = pair[0].split('_rep_')[1], pair[1].split('_rep_')[1]
            rscat.set_xlabel("Log"+u'\u2082'+" normalized reads in replicate " + rep1)
            rscat.set_ylabel("Log"+u'\u2082'+" normalized reads in replicate " + rep2)
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
        comparison = re.findall(r"_cond1_(.*?)_cond2_(.*?)_deseq_table.csv", os.path.basename(dgefile))

        if not comparison:
            raise ValueError("Could not find condition names in DGE filename: " + dgefile)
        if len(comparison) > 1:
            print("Warning: multiple conditions matched in DGE filename. Using first match.")

        comparison_name = "_vs_".join(comparison[0])
        de_table[comparison_name] = pd.read_csv(dgefile, index_col=0)['padj']

    return de_table


def scatter_dges(count_df, dges, output_prefix, viewLims, classes=None, show_unknown=True, pval=0.05):
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
            sscat.set_xlabel(f"Log"+u'\u2082'+" normalized reads in " + p1)
            sscat.set_ylabel(f"Log"+u'\u2082'+" normalized reads in " + p2)
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
            sscat.set_xlabel("Log"+u'\u2082'+" normalized reads in " + p1)
            sscat.set_ylabel("Log"+u'\u2082'+" normalized reads in " + p2)
            pdf_name = make_filename([output_prefix, pair, 'scatter_by_dge'], ext='.pdf')
            sscat.figure.savefig(pdf_name)


def load_raw_counts(raw_counts_file: str) -> Optional[pd.DataFrame]:
    """Loads a raw_counts CSV as a DataFrame
    Args:
        raw_counts_file: The raw counts CSV produced by Counter
    Returns:
        The raw counts DataFrame
    """

    if raw_counts_file is None: return None
    return pd.read_csv(raw_counts_file, index_col=0)


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


def get_class_counts(raw_counts_df: pd.DataFrame) -> pd.DataFrame:
    """Calculates class counts from the Raw Counts DataFrame

    If there are multiple classes associated with a feature, then that feature's
    counts are normalized across all associated classes.

    Args:
        raw_counts_df: A DataFrame containing features and their associated classes and counts
    Returns:
        class_counts: A DataFrame with an index of classes and count entries, per comparison
    """

    # 1. Group by Feature Class. Class "lists" are strings. Normalize their counts by the number of commas.
    grouped = (raw_counts_df.drop("Feature Name", 1, errors='ignore')
               .groupby("Feature Class")
               .apply(lambda grp: grp.sum() / (grp.name.count(',') + 1)))

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
        'class_charts': all([args.raw_counts, args.summary_stats]),
        'replicate_scatter': all([args.norm_counts]),
        'sample_avg_scatter_by_dge': all([args.norm_counts, args.dge_tables]),
        'sample_avg_scatter_by_dge_class': all([args.norm_counts, args.dge_tables]),
    }

    nc_file = "normalized count file (-nc/--norm-counts)"
    rc_file = "raw count file (-rc/--raw-counts)"
    sc_file = "summary stats file (-sc/--summary-stats)"
    dge_files = "differential expression tables (-dge/--dge-tables)"

    error_message = {
        'len_dist': "5' nucleotide/length tables (-len/--len-dist) are required for ",
        'class_charts': f"A {rc_file} and {sc_file} are required for ",
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
        'class_charts': ["raw_count_df", "class_counts", "mapped_reads"],
        'replicate_scatter': ["norm_count_df", "sample_rep_dict", "norm_view_lims"],
        'sample_avg_scatter_by_dge':
            ["norm_count_df", "sample_rep_dict", "norm_count_avg_df", "de_table", "avg_view_lims"],
        'sample_avg_scatter_by_dge_class':
            ["norm_count_df", "sample_rep_dict", "norm_count_avg_df", "feat_classes", "de_table", "avg_view_lims"],
    }

    relevant_vars: Dict[str, Union[pd.DataFrame, pd.Series, dict, None]] = {}
    input_getters = {
        'raw_count_df': lambda: load_raw_counts(args.raw_counts),
        'mapped_reads': lambda: get_mapped_reads(args.summary_stats),
        'norm_count_df': lambda: load_norm_counts(args.norm_counts),
        'de_table': lambda: load_dge_tables(args.dge_tables),
        'sample_rep_dict': lambda: get_sample_rep_dict(relevant_vars["norm_count_df"]),
        'norm_count_avg_df': lambda: get_sample_averages(relevant_vars["norm_count_df"], relevant_vars["sample_rep_dict"]),
        'feat_classes': lambda: get_flat_classes(relevant_vars["norm_count_df"]),
        'class_counts': lambda: get_class_counts(relevant_vars["raw_count_df"]),
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
            func = class_charts
            arg, kwd = (inputs["class_counts"], inputs['mapped_reads'], args.out_prefix), {}
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
                results.append(pool.apply_async(task, args, kwds, error_callback=err))
            for result in results:
                result.wait()
    elif len(itinerary) == 1:
        # Don't use multiprocessing if only one plot type requested
        task = itinerary.pop()
        task[0](*task[1], **task[2])

def err(e):
    """Allows us to print errors from a MP worker without discarding the other results"""
    print("An error occurred. See console log for details.", file=sys.stderr)
    print(''.join(traceback.format_exception(type(e), e, e.__traceback__)))

if __name__ == '__main__':
    main()

