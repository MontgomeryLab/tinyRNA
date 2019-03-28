#! /usr/bin/env python
"""
Merge the outputs of counts per sample into one large table.

This script takes the outputs from counter.py and combines them into a larger,
merged file in order to view 1) all the counts in one table 2) all the stats in
one table. The outputs can then be used for further analysis (DEG, plots, etc).
"""
import argparse
import pandas as pd

def get_args():
    """
    Get input arguments from the user/command line.

    Requires the user provide a list of count files from aquatx-count,
    associated sample names, output file name, and mode of merging.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-files', nargs='+', metavar='FILE', required=True,
                        help='input count files to merge for processing')
    parser.add_argument('-o', '--output-file', metavar='OUTPUT', required=True,
                        help='output filename')
    parser.add_argument('-s', '--sample-names', nargs='+', metavar='NAMES', required=True,
                        help='associated sample names for input files given')
    parser.add_argument('-m', '--mode', metavar='MODE', required=True,
                        help='mode for merging: counts, stats, or unique seq files')

    args = parser.parse_args()

    return args

def merge_counts(counts_files, samples):
    """
    Takes in list of feature count files and merges them together into
    one table for further processing and analysis.

    Inputs:
        counts_files: A list of files to merge
        samples: Sample names, ordered the same as counts_files

    Outputs:
        count_df: The final, merged data frame of feature counts
    """
    # Create the first data frame to build on
    temp_counts = pd.read_csv(counts_files[0], sep='\t', header=None, index_col=0)

    # Create an empty data frame based on input dimensions to fill in
    count_df = pd.DataFrame(index=temp_counts.index, columns=samples)
    count_df.loc[:, samples[0]] = temp_counts.values

    count = 1
    for cf in counts_files[1:]:
        temp_counts = pd.read_csv(cf, sep='\t', header=None, index_col=0).reindex(count_df.index)
        count_df.loc[:, samples[count]] = temp_counts.values
        count += 1

    count_df.index.rename('feature', inplace=True)

    return count_df

def merge_stats(stats_files, samples):
    """
    Takes in a list of stats files and merges them together into one
    summary statistics table for the entire run.

    Inputs:
        stats_files: A list of files to merge
        samples: Sample names, ordered the same as stats_files

    Outputs:
        align_df: Overall alignment statistics data frame
        feature_df: Feature counting statistics data frame
    """
    # Define the stats we want
    align_stats = ['_unique_sequences_aligned', '_aligned_reads',
                   '_aligned_reads_multi_mapping', '_aligned_reads_unique_mapping']
    feature_stats = ['_no_feature', '_ambiguous_alignments_classes',
                     '_ambiguous_reads_classes', '_ambiguous_alignments_features',
                     '_ambiguous_reads_features', '_alignments_unique_features',
                     '_reads_unique_features']

    # Create the first dataframe to build on
    temp_stats = pd.read_csv(stats_files[0], sep='\t', header=None, skiprows=1, index_col=0)

    # Create an empty data frame based on input dimensions to fill in
    stat_df = pd.DataFrame(index=temp_stats.index, columns=samples)
    stat_df.loc[:, samples[0]] = temp_stats.values

    count = 1
    for sf in stats_files[1:]:
        temp_stats = pd.read_csv(sf, sep='\t', header=None, skiprows=1, index_col=0).reindex(stat_df.index)
        stat_df.loc[:, samples[count]] = temp_stats.values
        count += 1

    align_df = stat_df.loc[stat_df.index.isin(align_stats)]
    feature_df = stat_df.loc[stat_df.index.isin(feature_stats)]

    return align_df, feature_df

def main():
    """ Main routine """
    # Step 1: Get the command line arguments
    args = get_args()

    # Step 2: Determine the merge mode
    if args.mode == 'counts':
        count_df = merge_counts(args.input_files, args.sample_names)
        count_df.to_csv(args.output_file)

    elif args.mode == 'stats':
        align_stat, feat_stat = merge_stats(args.input_files, args.sample_names)
        align_stat.index.name = 'Alignment Statistics'
        align_stat.to_csv(args.output_file, header=True, sep='\t')
        with open(args.output_file, 'a') as stat_out:
            stat_out.write('\n')
            feat_stat.index.name = 'Feature Statistics'
            feat_stat.to_csv(stat_out, header=True, sep='\t')

if __name__ == '__main__':
    main()
