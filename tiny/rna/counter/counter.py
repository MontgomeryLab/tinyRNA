"""This submodule assigns feature counts for SAM alignments using a Feature Sheet ruleset.

If you find that you are sourcing all of your input files from a prior run, we recommend
that you instead run `tiny recount` within that run's directory.
"""

import multiprocessing as mp
import traceback
import argparse
import sys
import os

from collections import defaultdict
from typing import Tuple, List, Dict

from tiny.rna.counter.features import Features, FeatureCounter
from tiny.rna.counter.statistics import SummaryStats
from tiny.rna.util import report_execution_time, from_here
from tiny.rna.configuration import CSVReader

# Global variables for multiprocessing
feature_counter: FeatureCounter


def get_args():
    """Get input arguments from the user/command line."""

    arg_parser = argparse.ArgumentParser(description=__doc__)
    required_group = arg_parser.add_argument_group("required arguments")

    # Required arguments
    required_group.add_argument('-i', '--input-csv', metavar='SAMPLES', required=True,
                        help='your Samples Sheet')
    required_group.add_argument('-c', '--config', metavar='CONFIGFILE', required=True,
                        help='your Features Sheet')
    required_group.add_argument('-o', '--out-prefix', metavar='OUTPUTPREFIX', required=True,
                        help='output prefix to use for file names')

    # Optional arguments
    arg_parser.add_argument('-sf', '--source-filter', metavar='SOURCE', nargs='*', default=[],
                        help='Only produce counts for features whose '
                             'GFF column 2 matches the source(s) listed')
    arg_parser.add_argument('-tf', '--type-filter', metavar='TYPE', nargs='*', default=[],
                        help='Only produce counts for features whose '
                             'GFF column 3 matches the type(s) listed')
    arg_parser.add_argument('-nn', '--no-normalize', action='store_true',
                        help='Do not normalize counts by genomic hits '
                             'and (selected) overlapping feature counts.')
    arg_parser.add_argument('-a', '--all-features', action='store_true',
                        help='Represent all features in output counts table, '
                             'regardless of counts or identity rules.')
    arg_parser.add_argument('-p', '--is-pipeline', action='store_true',
                        help='Indicates that counter was invoked as part of a pipeline run '
                             'and that input files should be sourced as such.')
    arg_parser.add_argument('-d', '--report-diags', action='store_true',
                        help='Produce diagnostic information about uncounted/eliminated '
                             'selection elements.')

    return arg_parser.parse_args()


def load_samples(samples_csv: str, is_pipeline: bool) -> List[Dict[str, str]]:
    """Parses the Samples Sheet to determine library names and alignment files for counting

    Sample files may have a .fastq(.gz) extension (i.e. when Counter is called as part of a
    pipeline run) or a .sam extension (i.e. when Counter is called as a standalone tool).

    Args:
        samples_csv: a csv file which defines sample group, replicate, and file location
        is_pipeline: helps locate sample SAM files. If true, files are assumed to reside
            in the working directory. If false, files are assumed to reside in the same
            directory as their source FASTQ files with '_aligned_seqs.sam' appended
            (i.e. /dir/sample1.fastq -> /dir/sample1_aligned_seqs.sam).

    Returns:
        inputs: a list of dictionaries for each library, where each dictionary defines the
        library name and the location of its SAM file for counting.
    """

    def get_library_filename(csv_row_file: str, samples_csv: str) -> str:
        """The input samples.csv may contain either fastq or sam files"""

        file_ext = os.path.splitext(csv_row_file)[1].lower()

        # If the sample file has a fastq(.gz) extension, infer the name of its pipeline-produced .sam file
        if file_ext in [".fastq", ".gz"]:
            # Fix relative paths to be relative to sample_csv's path, rather than relative to cwd
            csv_row_file = os.path.basename(csv_row_file) if is_pipeline else from_here(samples_csv, csv_row_file)
            csv_row_file = os.path.splitext(csv_row_file)[0] + "_aligned_seqs.sam"
        elif file_ext == ".sam":
            if not os.path.isabs(csv_row_file):
                raise ValueError("The following file must be expressed as an absolute path:\n%s" % (csv_row_file,))
        else:
            raise ValueError("The filenames defined in your Samples Sheet must have a .fastq(.gz) or .sam extension.\n"
                             "The following filename contained neither:\n%s" % (csv_row_file,))
        return csv_row_file

    inputs = list()

    for row in CSVReader(samples_csv, "Samples Sheet").rows():
        library_name = f"{row['Group']}_rep_{row['Replicate']}"
        library_file_name = get_library_filename(row['File'], samples_csv)
        record = {"Name": library_name, "File": library_file_name}

        if record not in inputs: inputs.append(record)

    return inputs


def load_config(features_csv: str, is_pipeline: bool) -> Tuple[List[dict], Dict[str, list]]:
    """Parses the Features Sheet to provide inputs to FeatureSelector and build_reference_tables

    Args:
        features_csv: a csv file which defines feature sources and selection rules
        is_pipeline: helps locate GFF files defined in the Features Sheet. If true,
            GFF files are assumed to reside in the working directory.

    Returns:
        rules: a list of dictionaries, each representing a parsed row from input.
            Note that these are just rule definitions which FeatureSelector will
            further digest to produce its rules table.
        gff_files: a dict of GFF files and associated Name Attribute preferences
    """

    rules, gff_files = list(), defaultdict(list)

    for row in CSVReader(features_csv, "Features Sheet").rows():
        rule = {col: row[col] for col in ["Strand", "Hierarchy", "nt5end", "Length", "Strict"]}
        rule['nt5end'] = rule['nt5end'].upper().translate({ord('U'): 'T'})  # Convert RNA base to cDNA base
        rule['Identity'] = (row['Key'], row['Value'])                       # Create identity tuple
        rule['Hierarchy'] = int(rule['Hierarchy'])                          # Convert hierarchy to number
        rule['Strict'] = rule['Strict'] == 'Full'                           # Convert strict intersection to boolean

        gff = os.path.basename(row['Source']) if is_pipeline else from_here(features_csv, row['Source'])

        # Duplicate Name Attributes and rule entries are not allowed
        if row['Name'] not in ["ID", *gff_files[gff]]: gff_files[gff].append(row['Name'])
        if rule not in rules: rules.append(rule)

    rules.sort(key=lambda x: x['Hierarchy'])
    return rules, gff_files


@report_execution_time("Counting and merging")
def map_and_reduce(libraries):
    """Assigns one worker process per library and merges the statistics they report"""

    # SummaryStats handles final output files, regardless of multiprocessing status
    summary = SummaryStats(Features, FeatureCounter.out_prefix, FeatureCounter.run_diags)

    # Use a multiprocessing pool if multiple sam files were provided
    if len(libraries) > 1:
        with mp.Pool(len(libraries)) as pool:
            async_results = pool.imap_unordered(feature_counter.count_reads, libraries)

            for stats_result in async_results:
                summary.add_library(stats_result)
    else:
        # Only one library, multiprocessing not beneficial for task
        summary.add_library(feature_counter.count_reads(libraries[0]))

    return summary


@report_execution_time("Counter's overall runtime")
def main():
    # Get command line arguments.
    args = get_args()

    try:
        # Determine SAM inputs and their associated library names
        libraries = load_samples(args.input_csv, args.is_pipeline)

        # Load selection rules and feature sources from config
        selection_rules, gff_file_set = load_config(args.config, args.is_pipeline)

        # global for multiprocessing
        global feature_counter
        feature_counter = FeatureCounter(gff_file_set, selection_rules, **vars(args))

        # Assign and count features using multiprocessing and merge results
        merged_counts = map_and_reduce(libraries)

        # Write final outputs
        merged_counts.write_report_files(Features.aliases)
    except:
        traceback.print_exception(*sys.exc_info())
        print("\n\nCounter encountered an error. Don't worry! You don't have to start over.\n"
              "You can resume the pipeline at Counter. To do so:\n\t"
              "1. cd into your Run Directory\n\t"
              '2. Run "tiny recount --config your_run_config.yml"\n\t'
              '   (that\'s the processed run config) ^^^\n\n', file=sys.stderr)


if __name__ == '__main__':
    main()
