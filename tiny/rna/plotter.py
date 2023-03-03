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
import re

from collections import defaultdict
from typing import Dict, Union, Tuple, DefaultDict, Iterable
from pkg_resources import resource_filename

from tiny.rna.plotterlib import plotterlib
from tiny.rna.util import report_execution_time, make_filename, SmartFormatter, timestamp_format, add_transparent_help

aqplt: plotterlib = None
RASTER: bool = True


def get_args():
    """Get input arguments from the user/command line."""

    parser = argparse.ArgumentParser(description=__doc__, add_help=False, formatter_class=SmartFormatter)
    required_args = parser.add_argument_group("Required arguments")
    counter_files = parser.add_argument_group("Input files produced by tiny-count")
    diffexp_files = parser.add_argument_group("Input files produced by tiny-deseq.r")
    optional_args = parser.add_argument_group("Optional arguments")

    # Single file inputs
    counter_files.add_argument('-rc', '--raw-counts', metavar='RAW_COUNTS',
                               help='The ...feature_counts.csv file')
    diffexp_files.add_argument('-nc', '--norm-counts', metavar='NORM_COUNTS',
                               help='The ...norm_counts.csv file')
    counter_files.add_argument('-uc', '--rule-counts', metavar='RULE_COUNTS',
                               help='The ...counts-by-rule.csv file')
    counter_files.add_argument('-ss', '--summary-stats', metavar='STAT',
                               help='The ...summary_stats.csv file')

    # Multi-file inputs
    diffexp_files.add_argument('-dge', '--dge-tables', metavar='COMPARISON', nargs='+',
                               help='The ...cond1...cond2...deseq.csv files')
    counter_files.add_argument('-len', '--len-dist', metavar='5P_LEN', nargs='+',
                               help='The ...nt_len_dist.csv files')

    # Outputs options
    optional_args.add_argument('-o', '--out-prefix', metavar='PREFIX',
                               help='Prefix to use for output filenames.')
    optional_args.add_argument('-pv', '--p-value', metavar='VALUE', default=0.05, type=float,
                               help='P value to use in DGE scatter plots.')
    optional_args.add_argument('-s', '--style-sheet', metavar='MPLSTYLE',
                               default=resource_filename('tiny', 'templates/tinyrna-light.mplstyle'),
                               help='Optional matplotlib style sheet to use for plots.')
    optional_args.add_argument('-v', '--vector-scatter', action='store_true',
                               help='Produce scatter plots with vectorized points (slower).\n'
                               'Note: only the points on scatter plots will be raster if '
                               'this option is not provided.')
    optional_args.add_argument('-ldi', '--len-dist-min', metavar='VALUE', type=int,
                               help='len_dist plots will start at this value')
    optional_args.add_argument('-lda', '--len-dist-max', metavar='VALUE', type=int,
                               help='len_dist plots will end at this value')
    optional_args.add_argument('-dgi', '--dge-min', metavar='VALUE', type=float,
                               help='scatter_by_dge plots will start at this log2 value')
    optional_args.add_argument('-dga', '--dge-max', metavar='VALUE', type=float,
                               help='scatter_by_dge plots will end at this log2 value')
    optional_args.add_argument('-una', '--unassigned-class', metavar='LABEL', default='_UNASSIGNED_',
                               help='Use this label in class-related plots for unassigned counts'),
    optional_args.add_argument('-unk', '--unknown-class', metavar='LABEL', default='_UNKNOWN_',
                               help='Use this label in class-related plots for counts which were '
                                    'assigned by rules lacking a "Classify as..." value')

    # Class filtering options
    mutex_class_filter = optional_args.add_mutually_exclusive_group()
    mutex_class_filter.add_argument('-ic', '--classes-include', metavar='CLASS', nargs='+', type=str,
                                    help='Only include these classes, if present, in class scatter '
                                         'plots (applies regardless of P value)')
    mutex_class_filter.add_argument('-ec', '--classes-exclude', metavar='CLASS', nargs='+', type=str,
                                    help='Omit these classes, if present, from class scatter plots '
                                         '(applies regardless of P value)')

    # Required arguments
    required_args.add_argument('-p', '--plots', metavar='PLOT', required=True, nargs='+',
                               help="R|List of plots to create. Options: \n"
                               "• len_dist: A stacked barchart showing size & 5' nucleotide distribution.\n"
                               "• rule_charts: A barchart showing percentages of counts by matched rule.\n"
                               "• class_charts: A barchart showing percentages of counts per classification.\n"
                               "• replicate_scatter: A scatter plot comparing replicates for all count files given.\n"
                               "• sample_avg_scatter_by_dge: A scatter plot comparing all sample groups, with "
                                    "differentially expressed small RNAs highlighted based on P value cutoff.\n"
                               "• sample_avg_scatter_by_dge_class: A scatter plot comparing all sample groups, "
                                    "with classes highlighted for differentially expressed small RNAs "
                                    "based on P value cutoff.")

    add_transparent_help(parser)
    return parser.parse_args()


def len_dist_plots(matrices: dict, out_prefix:str, vmin: int = None, vmax: int = None, **kwargs):
    """Create a PDF of size and 5'nt distribution plot for a sample.

    Lengths are plotted as a continuous closed interval from vmin to vmax. If either are
    unspecified, min and/or max are calculated across libraries. When this happens,
    bounds are determined separately per plot subtype (i.e. the min/max of Assigned
    len_dists does not determine the min/max of Mapped len_dists)

    Args:
        matrices: A dictionary of size + 5p-nt counts per library per subtype
        out_prefix: The prefix to use when naming output PDF plots
        vmin: The optional first length to plot in the range
        vmax: The optional last length to plot in the range
        kwargs: Additional keyword arguments to pass to pandas.DataFrame.plot()
    """

    orig_vmin, orig_vmax = vmin, vmax

    for subtype, libraries in matrices.items():
        if orig_vmin is None: vmin = min([l.index.min() for l in libraries.values()])
        if orig_vmax is None: vmax = max([l.index.max() for l in libraries.values()])
        if vmin > vmax:
            print(f'ERROR: len_dist min > max. Skipping "{subtype}" plot subtype', file=sys.stderr)
            continue

        for condition_and_rep, len_dist_df in libraries.items():
            # Convert reads to proportion
            size_prop = len_dist_df / len_dist_df.sum().sum()

            # Ensure that there are no gaps & close range's half-open interval
            size_prop = size_prop.reindex(range(vmin, vmax + 1)).fillna(0)

            # Create the plot
            plot = aqplt.len_dist_bar(size_prop, subtype, **kwargs)

            # Save plot
            pdf_name = make_filename([out_prefix, condition_and_rep, "len_dist"], ext='.pdf')
            save_plot(plot, "len_dist", pdf_name)


def get_len_dist_dict(files_list: list) -> DefaultDict[str, Dict[str, pd.DataFrame]]:
    """Reads 5' nt/len matrices into a dictionary of DataFrames

    Args:
        files_list: a list of 5' nt/len matrices which follow the naming convention
            {sample}_rep_{replicate}_{subtype}_nt_len_dist.csv

    Returns: A dictionary of size + 5p-nt counts per library per subtype
    """

    matrices = defaultdict(dict)

    for file in sorted(files_list):
        # Parse the "sample_rep_N" string from the input filename to avoid duplicate out_prefixes in the basename
        basename = os.path.splitext(os.path.basename(file))[0]
        date_prefix_pos = re.search(timestamp_format, basename)

        if date_prefix_pos is not None:
            # File is a pipeline product
            begin = date_prefix_pos.end() + 1
            end = basename.rfind("_nt_len_dist")
            condition_and_rep = basename[begin:end]
        else:
            # File does not appear to have been produced by the pipeline
            condition_and_rep = basename

        subtype = "assigned" if "assigned" in condition_and_rep else "mapped"
        matrices[subtype][condition_and_rep] = pd.read_csv(file, index_col=0)

    return matrices


def class_charts(raw_class_counts: pd.DataFrame, mapped_reads: pd.Series, out_prefix: str, class_na: str, scale=2, **kwargs):
    """Create a PDF of the proportion of raw assigned counts vs total mapped reads for each class

    Args:
        raw_class_counts: A dataframe containing RAW class counts per library
        mapped_reads: A Series containing the total number of mapped reads per library
        out_prefix: The prefix to use when naming output PDF plots
        class_na: The label to use for representing unassigned reads
        scale: The decimal scale for table percentages, and for determining inclusion of "unassigned"
        kwargs: Additional keyword arguments to pass to pandas.DataFrame.plot()
    """

    class_props = get_proportions_df(raw_class_counts, mapped_reads, class_na, scale).sort_index()
    max_prop = class_props.max().max()

    for library in raw_class_counts:
        chart = aqplt.barh_proportion(class_props[library], max_prop, scale, **kwargs)
        chart.set_title("Percentage of small RNAs by class")
        chart.set_ylabel("Class")

        # Save the plot
        pdf_name = make_filename([out_prefix, library, 'class_chart'], ext='.pdf')
        save_plot(chart, "class_chart", pdf_name)


def rule_charts(rule_counts: pd.DataFrame, out_prefix: str, scale=2, **kwargs):
    """Create a PDF of the proportion of raw assigned counts vs total mapped reads for each matched rule

    Args:
        rule_counts: A dataframe containing raw rule counts per library
        out_prefix: The prefix to use when naming output PDF plots
        scale: The decimal scale for table percentages, and for determining inclusion of "unassigned"
        **kwargs: Additional keyword arguments to pass to pandas.DataFrame.plot()
    """

    rule_strings = rule_counts.loc[:, 'Rule String']
    rule_counts.drop('Rule String', axis=1, inplace=True)
    mapped_reads = rule_counts.loc['Mapped Reads']
    rule_counts.drop('Mapped Reads', inplace=True)

    rule_props = get_proportions_df(rule_counts, mapped_reads, "N", scale)
    max_prop = rule_props.max().max()

    for library, prop_df in rule_props.items():
        chart = aqplt.barh_proportion(prop_df, max_prop, scale, **kwargs)
        chart.set_title("Percentage of small RNAs by matched rule")
        chart.set_ylabel("Rule")

        # Save the plot
        pdf_name = make_filename([out_prefix, library, 'rule_chart'], ext='.pdf')
        save_plot(chart, "rule_chart", pdf_name)


def get_proportions_df(counts_df: pd.DataFrame, mapped_totals: pd.Series, un: str, scale=2):
    """Calculates proportions and unassigned counts, rounded according to scale
    Unassigned counts below scale threshold will be replaced with NaNs

    Args:
        counts_df: A dataframe of counts with library columns
        mapped_totals: A series of mapped totals indexed by library
        un: The name to use for the unassigned counts proportion
        scale: The desired number of *percentage* decimal places"""

    scale += 2  # Convert percentage scale to decimal scale
    props = (counts_df / mapped_totals).round(scale)

    unassigned_threshold = round(10 ** (-1 * scale), scale)
    un_ds = (1 - props.sum()).round(scale).rename(un)
    un_ds = un_ds.mask(un_ds < unassigned_threshold)
    return pd.concat([pd.DataFrame(un_ds).T, props])


def load_mapped_reads(summary_stats_file: str) -> pd.Series:
    """Produces a Series of mapped reads per library for calculating proportions in class charts.

    Args:
        summary_stats_file: the summary stats csv produced by tiny-count

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
    non_numeric_cols = ["Feature Name"]

    for col in df.columns:
        if col in non_numeric_cols: continue
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


def scatter_replicates(count_df: pd.DataFrame, samples: dict, output_prefix: str,
                       view_lims: Tuple[float, float] = None):
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
                                         color='#B3B3B3', alpha=0.5, rasterized=RASTER)
            aqplt.set_square_scatter_view_lims(rscat, view_lims)
            aqplt.set_scatter_ticks(rscat)
            rscat.set_title(samp)
            rep1, rep2 = pair[0].split('_rep_')[1], pair[1].split('_rep_')[1]
            rscat.set_xlabel("Log$_{2}$ normalized reads in replicate " + rep1)
            rscat.set_ylabel("Log$_{2}$ normalized reads in replicate " + rep2)
            pdf_name = make_filename([output_prefix, samp, 'replicates', rep1, rep2, 'scatter'], ext='.pdf')
            save_plot(rscat, "replicate_scatter", pdf_name)


def load_dge_tables(comparisons: list, class_fillna: str) -> pd.DataFrame:
    """Creates a new dataframe containing all features and padj values for each comparison
    from a list of DGE tables.

    Args:
        comparisons: DGE tables (files) produced by DESeq2
        class_fillna: the value to fill for NaN values in multiindex level 1

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
        table = set_counts_table_multiindex(pd.read_csv(dgefile), class_fillna)

        de_table[comparison_name] = table['padj']

    return de_table


def filter_dge_classes(count_df: pd.DataFrame, dges: pd.DataFrame, include: Iterable = None, exclude: Iterable = None):
    """Filters features by classification in counts and DGE tables.
    Arguments `include` and `exclude` are mutually exclusive; providing both
    arguments will result in an error.

    Args:
        count_df: A dataframe with an index comprised of (feature ID, classifier)
        dges: A dataframe with an index comprised of (feature ID, classifier)
        include: An iterable of classifiers to allow
        exclude: An iterable of classifiers to exclude

    Returns: A filtered counts_avg_df and a filtered dge table
    """

    if not (include or exclude):
        return count_df, dges
    elif include and exclude:
        raise ValueError("Include/exclude filters are mutually exclusive.")

    if include is not None:
        include_lc = [cls.lower() for cls in include]
        mask = count_df.index.get_level_values('Classifier').str.lower().isin(include_lc)
    elif exclude is not None:
        exclude_lc = [cls.lower() for cls in exclude]
        mask = ~count_df.index.get_level_values('Classifier').str.lower().isin(exclude_lc)
    else:
        mask = None  # to appease linter

    return count_df[mask], dges[mask]


def scatter_by_dge_class(counts_avg_df, dges, output_prefix, view_lims, include=None, exclude=None, pval=0.05):
    """Creates PDFs of all pairwise comparison scatter plots with differentially
    expressed features colored by class. Counts for features with P value >= `pval`
    will be assigned the color grey.

    Args:
        counts_avg_df: A dataframe of normalized counts per multiindex (feature ID, classification) that
            have been averaged across replicates within each condition.
        dges: A differential gene expression dataframe with multiindex (feature ID, classification),
            and columns containing adjusted P values for each pairwise comparison. Each column should be
            labeled as "conditionA_vs_conditionB" where conditionB is the untreated condition.
        view_lims: A tuple of (min, max) data bounds for the plot view
        include: A list of classes to plot, if present (default: all classes)
        exclude: A list of classes to omit from plots (default: no classes)
        output_prefix: A string to use as a prefix for saving files
        pval: The P value threshold for determining the outgroup
    """

    counts_avg_df, dges = filter_dge_classes(counts_avg_df, dges, include, exclude)
    if counts_avg_df.empty or dges.empty:
        print('ERROR: No classes passed filtering. Skipping scatter_by_dge_class.', file=sys.stderr)
        return

    uniq_classes = pd.unique(counts_avg_df.index.get_level_values(1))
    class_colors = aqplt.assign_class_colors(uniq_classes)
    aqplt.set_dge_class_legend_style()

    for pair in dges:
        tr, ut = pair.split("_vs_")  # treated, untreated
        dge_classes = dges[dges[pair] < pval].groupby(level=1).groups

        labels, grp_args = zip(*dge_classes.items()) if dge_classes else ((), ())
        sscat = aqplt.scatter_grouped(counts_avg_df.loc[:, ut], counts_avg_df.loc[:, tr], *grp_args,
                                      colors=class_colors, pval=pval, view_lims=view_lims, labels=labels,
                                      rasterized=RASTER)

        sscat.set_title('%s vs %s' % (tr, ut))
        sscat.set_xlabel("Log$_{2}$ normalized reads in " + ut)
        sscat.set_ylabel("Log$_{2}$ normalized reads in " + tr)
        sscat.get_legend().set_bbox_to_anchor((1, 1))
        pdf_name = make_filename([output_prefix, pair, 'scatter_by_dge_class'], ext='.pdf')
        save_plot(sscat, "scatter_by_dge_class", pdf_name)


def scatter_by_dge(counts_avg_df, dges, output_prefix, view_lims, pval=0.05):
    """Creates PDFs of all pairwise comparison scatter plots with differentially
    expressed features highlighted. Counts for features with P value >= `pval`
    will be assigned the color grey.

    Args:
        counts_avg_df: A dataframe of normalized counts per multiindex (feature ID, classification) that
            have been averaged across replicates within each condition.
        dges: A differential gene expression dataframe with multiindex (feature ID, classification),
            and columns containing adjusted P values for each pairwise comparison. Each column should be
            labeled as "conditionA_vs_conditionB" where conditionB is the untreated condition.
        view_lims: A tuple of (min, max) data bounds for the plot view
        output_prefix: A string to use as a prefix for saving files
        pval: The pvalue threshold for determining the outgroup
    """

    if counts_avg_df.empty or dges.empty:
        print('ERROR: Received empty counts data. Skipping scatter_by_dge.', file=sys.stderr)
        return

    for pair in dges:
        grp_args = dges.index[dges[pair] < pval]
        tr, ut = pair.split("_vs_")  # treated, untreated

        labels = ['p < %g' % pval] if not grp_args.empty else []
        colors = aqplt.assign_class_colors(labels)
        sscat = aqplt.scatter_grouped(counts_avg_df.loc[:, ut], counts_avg_df.loc[:, tr], grp_args,
                                      colors=colors, alpha=0.5, pval=pval, view_lims=view_lims, labels=labels,
                                      rasterized=RASTER)

        sscat.set_title('%s vs %s' % (tr, ut))
        sscat.set_xlabel("Log$_{2}$ normalized reads in " + ut)
        sscat.set_ylabel("Log$_{2}$ normalized reads in " + tr)
        pdf_name = make_filename([output_prefix, pair, 'scatter_by_dge'], ext='.pdf')
        save_plot(sscat, 'scatter_by_dge', pdf_name)


def load_raw_counts(raw_counts_file: str, class_fillna: str) -> pd.DataFrame:
    """Loads a raw_counts CSV as a DataFrame
    Args:
        raw_counts_file: The raw counts CSV produced by tiny-count
        class_fillna: the value to fill for NaN values in multiindex level 1
    Returns:
        The raw counts DataFrame with a multiindex of (feature id, classifier)
    """

    return set_counts_table_multiindex(pd.read_csv(raw_counts_file), class_fillna)


def load_rule_counts(rule_counts_file: str) -> pd.DataFrame:
    """Loads a rule_counts CSV as a DataFrame
    Args:
        rule_counts_file: The rule counts CSV produced by tiny-count
    Returns:
        The raw counts DataFrame with a multiindex of (feature id, classifier)
    """

    return pd.read_csv(rule_counts_file, index_col=0)


def load_norm_counts(norm_counts_file: str, class_fillna: str) -> pd.DataFrame:
    """Loads a norm_counts CSV as a DataFrame
    Args:
        norm_counts_file: The norm counts CSV produced by DESeq2
        class_fillna: the value to fill for NaN values in multiindex level 1
    Returns:
        The raw counts DataFrame with a multiindex of (feature id, classifier)
    """

    return set_counts_table_multiindex(pd.read_csv(norm_counts_file), class_fillna)


def set_counts_table_multiindex(counts: pd.DataFrame, fillna: str) -> pd.DataFrame:
    """Assigns a multiindex composed of (Feature ID, Classifier) to the counts table,
    and fills NaN values in the classifier column"""

    level0 = "Feature ID"
    level1 = "Classifier"

    counts[level1] = counts[level1].fillna(fillna)
    return counts.set_index([level0, level1])


def get_class_counts(raw_counts_df: pd.DataFrame) -> pd.DataFrame:
    """Calculates class counts from level 1 of the raw counts multiindex

    Args:
        raw_counts_df: A DataFrame with a multiindex of (feature ID, classifier)
    Returns:
        class_counts: A DataFrame with an index of classes and count entries, per comparison
    """

    return (raw_counts_df
            .drop("Feature Name", axis="columns", errors='ignore')
            .groupby(level=1)
            .sum())


def save_plot(plot, dir_name, filename):
    """Complete plots are sent here for output patterns that apply to all plots
    Args:
        - plot: the finished axes object
        - dir_name: all plots are placed in subdirectories by this name
        - filename: the file basename for the plot, including extension
    """

    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)

    final_save_path = os.path.join(dir_name, filename)
    plot.figure.savefig(final_save_path)


def validate_inputs(args: argparse.Namespace) -> None:
    """Determines if the necessary input files have been provided for the requested plots
    This is necessary because we allow users to run the tool with only the files necessary
    for the plot types they've requested, rather than all possible inputs.

    Args:
        args: Command line arguments parsed by argparse
    """

    dependencies_satisfied_for = {
        'len_dist': all([args.len_dist]),
        'rule_charts': all([args.rule_counts]),
        'class_charts': all([args.raw_counts, args.summary_stats]),
        'replicate_scatter': all([args.norm_counts]),
        'sample_avg_scatter_by_dge': all([args.norm_counts, args.dge_tables]),
        'sample_avg_scatter_by_dge_class': all([args.norm_counts, args.dge_tables]),
    }

    nc_file = "normalized count file (-nc/--norm-counts)"
    rc_file = "raw count file (-rc/--raw-counts)"
    ss_file = "summary stats file (-sc/--summary-stats)"
    dge_files = "differential expression tables (-dge/--dge-tables)"

    error_message = {
        'len_dist': "5' nucleotide/length tables (-len/--len-dist) are required for ",
        'rule_charts': "A rule counts file (-uc/--rule-counts) is required for ",
        'class_charts': f"A {rc_file} and {ss_file} are required for ",
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
        'len_dist':
            ["ld_matrices_dict"],
        'rule_charts':
            ["rule_counts_df"],
        'class_charts':
            ["raw_counts_df", "class_counts_df", "mapped_reads_ds"],
        'replicate_scatter':
            ["norm_counts_df", "sample_rep_dict", "norm_view_lims"],
        'sample_avg_scatter_by_dge':
            ["norm_counts_df", "sample_rep_dict", "norm_counts_avg_df",
             "de_table_df", "avg_view_lims"],
        'sample_avg_scatter_by_dge_class':
            ["norm_counts_df", "sample_rep_dict", "norm_counts_avg_df",
             "de_table_df", "avg_view_lims"]
    }

    # These are frozen function pointers; both the function and its
    # parameters aren't called until they are supplied the call operator ()
    fetched: Dict[str, Union[pd.DataFrame, pd.Series, dict, None]] = {}
    input_getters = {
        'ld_matrices_dict': lambda: get_len_dist_dict(args.len_dist),
        'mapped_reads_ds': lambda: load_mapped_reads(args.summary_stats),
        'rule_counts_df': lambda: load_rule_counts(args.rule_counts),
        'norm_counts_df': lambda: load_norm_counts(args.norm_counts, args.unknown_class),
        'raw_counts_df': lambda: load_raw_counts(args.raw_counts, args.unknown_class),
        'de_table_df': lambda: load_dge_tables(args.dge_tables, args.unknown_class),
        'sample_rep_dict': lambda: get_sample_rep_dict(fetched["norm_counts_df"]),
        'norm_counts_avg_df': lambda: get_sample_averages(fetched["norm_counts_df"], fetched["sample_rep_dict"]),
        'class_counts_df': lambda: get_class_counts(fetched["raw_counts_df"]),
        'avg_view_lims': lambda: aqplt.get_scatter_view_lims(fetched["norm_counts_avg_df"], args.dge_min, args.dge_max),
        'norm_view_lims': lambda: aqplt.get_scatter_view_lims(fetched["norm_counts_df"].select_dtypes(['number']))
    }

    # Input getters are executed here
    for plot in args.plots:
        for req in required_inputs[plot]:
            if req is not None and req not in fetched:
                fetched[req] = input_getters[req]()

    return fetched


@report_execution_time("Plotter runtime")
def main():

    args = get_args()
    validate_inputs(args)

    global aqplt, RASTER
    aqplt = plotterlib(args.style_sheet)
    RASTER = not args.vector_scatter
    inputs = setup(args)

    # Assemble work units for multiprocessing pool
    itinerary = []

    # generate plots requested
    for plot in args.plots:
        if plot == 'len_dist':
            func = len_dist_plots
            arg = (inputs["ld_matrices_dict"], args.out_prefix, args.len_dist_min, args.len_dist_max)
            kwd = {}
        elif plot == 'class_charts':
            func = class_charts
            arg = (inputs["class_counts_df"], inputs['mapped_reads_ds'], args.out_prefix, args.unassigned_class)
            kwd = {}
        elif plot == 'rule_charts':
            func = rule_charts
            arg = (inputs["rule_counts_df"], args.out_prefix)
            kwd = {}
        elif plot == 'replicate_scatter':
            func = scatter_replicates
            arg = (inputs["norm_counts_df"], inputs["sample_rep_dict"], args.out_prefix, inputs["norm_view_lims"])
            kwd = {}
        elif plot == 'sample_avg_scatter_by_dge':
            func = scatter_by_dge
            arg = (inputs["norm_counts_avg_df"], inputs["de_table_df"], args.out_prefix, inputs["avg_view_lims"])
            kwd = {"pval": args.p_value}
        elif plot == 'sample_avg_scatter_by_dge_class':
            func = scatter_by_dge_class
            arg = (inputs["norm_counts_avg_df"], inputs["de_table_df"], args.out_prefix, inputs["avg_view_lims"])
            kwd = {"pval": args.p_value, "include": args.classes_include, "exclude": args.classes_exclude}
        else:
            print('Plot type %s not recognized, please check the -p/--plot arguments' % plot)
            continue

        itinerary.append((func, arg, kwd))

    if len(itinerary) > 1 and not aqplt.is_debug_mode():
        mp.set_start_method('fork')
        with mp.Pool(len(itinerary)) as pool:
            results = []
            for task, args, kwds in itinerary:
                sentry = ExceptionManager(task)
                results.append(pool.apply_async(task, args, kwds, error_callback=sentry))
            for result in results:
                result.wait()
    else:
        # Don't use multiprocessing if only one plot type requested
        # or if in debug mode (for matplotlib compatibility)
        for task, args, kwds in itinerary:
            try:
                task(*args, **kwds)
            except Exception as e:
                ExceptionManager.add(task, e)

    ExceptionManager.print_exceptions()


class ExceptionManager:
    """Handles exception formatting for more user-friendly logging with cwltool

    In multiprocessing mode, you should create an instance for each task
    (plot type) and exceptions will be stored at the class level for ALL tasks.
    In sequential mode, you should use the add() method

    Exception tracebacks are printed to stdout as soon as they happen, and since
    the CWL CommandLineTool for tiny-plot captures stdout, this goes to the log
    file rather than terminal. Once plotting is complete, the user-friendly
    error is printed to stderr which isn't captured, so the user sees it.
    The message includes instructions for `tiny replot` followed by an
    exception summary (sans noisy traceback), organized by task."""

    excs = defaultdict(list)

    def __init__(self, task):
        self.task = task

    def __call__(self, e):
        """The multiprocessing error_callback target"""
        self.add(self.task, e, from_mp_worker=True)

    @classmethod
    def add(cls, task, e, from_mp_worker=False):
        """Prints task's traceback to stdout and stores exceptions for summary"""

        if from_mp_worker:
            print(e.__cause__.tb)
        else:
            ex_tuple = (type(e), e, e.__traceback__)
            traceback.print_exception(*ex_tuple)

        cls.excs[task].extend(traceback.format_exception_only(type(e), e))

    @classmethod
    def print_exceptions(cls):
        """Prints exception summary to stderr for all tasks"""

        if not cls.excs: return
        print('\n'.join(['', '=' * 75, '=' * 75]), file=sys.stderr)
        print("\nPlotter encountered an error. Don't worry! You don't have to start over.\n"
              "You can resume the pipeline at tiny-plot. To do so:\n\t"
              "1. cd into your Run Directory\n\t"
              '2. Run "tiny replot --config your_run_config.yml"\n\t'
              '   (that\'s the processed run config) ^^^\n', file=sys.stderr)

        ex_sum = sum(len(ex) for ex in cls.excs.values())
        header = "The following {} reported:"
        plural = "exceptions were" if ex_sum > 1 else "exception was"

        exc_list = [header.format(plural)]
        for task, task_exceptions in cls.excs.items():
            exc_list.append('\t' + f"In function {task.__name__}():")
            exc_list.append('\t\t'.join(['', *task_exceptions]))

        print('\n'.join(exc_list), file=sys.stderr)


if __name__ == '__main__':
    main()
