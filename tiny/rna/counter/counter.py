import multiprocessing as mp
import traceback
import argparse
import HTSeq
import queue
import sys
import csv
import os

from typing import Tuple, List
from collections import defaultdict

import tiny.rna.counter.hts_parsing as parser
from tiny.rna.counter.feature_selector import FeatureSelector
from tiny.rna.counter.statistics import LibraryStats, SummaryStats
from tiny.rna.counter.hts_parsing import SelectionRules, FeatureSources
from tiny.rna.util import report_execution_time, from_here

# Global variables for multiprocessing
feature_counter = None
is_pipeline = False
run_diags = False

# Type aliases for human readability
AssignedFeatures = set
N_Candidates = int


def get_args():
    """Get input arguments from the user/command line."""

    arg_parser = argparse.ArgumentParser()
    required_group = arg_parser.add_argument_group("required arguments")

    # Required arguments
    required_group.add_argument('-i', '--input-csv', metavar='SAMPLES', required=True,
                        help='the csv samples file/library list')
    required_group.add_argument('-c', '--config', metavar='CONFIGFILE', required=True,
                        help='the csv features configuration file')
    required_group.add_argument('-o', '--out-prefix', metavar='OUTPUTPREFIX', required=True,
                        help='output prefix to use for file names')

    # Optional arguments
    arg_parser.add_argument('-p', '--is-pipeline', action='store_true',
                        help='Indicates that counter was invoked from the tinyrna pipeline '
                             'and that input files should be sources as such.')
    arg_parser.add_argument('-d', '--diagnostics', action='store_true',
                        help='Produce diagnostic information about uncounted/eliminated '
                             'selection elements.')

    return arg_parser.parse_args()


@report_execution_time("Counting and merging")
def map_and_reduce(libraries, work_args, ret_queue):
    """Assigns one worker process per library and merges the statistics they report"""

    # Use a multiprocessing pool if multiple sam files were provided
    # Otherwise perform counts in this process
    if len(libraries) > 1:
        with mp.Pool(len(libraries)) as pool:
            pool.starmap(feature_counter.count_reads, work_args)
    else:
        feature_counter.count_reads(*work_args[0])

    # Collect counts from all pool workers and merge
    out_prefix = work_args[0][-1]
    summary = SummaryStats(out_prefix, feature_counter.run_diags)
    for _ in libraries:
        lib_stats = ret_queue.get()
        summary.add_library(lib_stats)

    return summary


def load_samples(samples_csv: str) -> list:
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
            raise ValueError(
                "The filenames defined in your Samples Sheet must have a .fastq(.gz) or .sam extension.\n"
                "The following filename contained neither:\n%s" % (csv_row_file,))
        return csv_row_file

    inputs = list()

    with open(os.path.expanduser(samples_csv), 'r', encoding='utf-8-sig') as f:
        fieldnames = ("File", "Group", "Replicate")
        csv_reader = csv.DictReader(f, fieldnames=fieldnames, delimiter=',')

        next(csv_reader)  # Skip header line
        for row in csv_reader:
            library_name = f"{row['Group']}_rep_{row['Replicate']}"
            library_file_name = get_library_filename(row['File'], samples_csv)
            record = {"Name": library_name, "File": library_file_name}

            if record not in inputs: inputs.append(record)

    return inputs


def load_config(features_csv: str) -> Tuple[SelectionRules, FeatureSources]:
    """Parses features.csv to provide inputs to FeatureSelector and build_reference_tables

    Args:
        features_csv: a csv file which defines feature sources and selection rules

    Returns:
        rules: a list of dictionaries, each representing a parsed row from input
        gff_files: a dict of GFF files and associated Name Attribute preferences
    """

    rules, gff_files = list(), defaultdict(list)
    convert_strand = {'sense': tuple('+'), 'antisense': tuple('-'), 'both': ('+', '-')}

    with open(os.path.expanduser(features_csv), 'r', encoding='utf-8-sig') as f:
        fieldnames = ("Name", "Key", "Value", "Hierarchy", "Strand", "nt5", "Length", "Strict", "Source")
        csv_reader = csv.DictReader(f, fieldnames=fieldnames, delimiter=',')

        next(csv_reader)  # Skip header line
        for row in csv_reader:
            rule = {col: row[col] for col in ["Strand", "Hierarchy", "nt5", "Length", "Strict"]}
            rule['nt5'] = rule['nt5'].upper().translate({ord('U'): 'T'})  # Convert RNA base to cDNA base
            rule['Strand'] = convert_strand[rule['Strand'].lower()]  # Convert sense/antisense to +/-
            rule['Identity'] = (row['Key'], row['Value'])  # Create identity tuple
            rule['Hierarchy'] = int(rule['Hierarchy'])  # Convert hierarchy to number
            rule['Strict'] = rule['Strict'] == 'Full'  # Convert strict intersection to boolean

            gff = os.path.basename(row['Source']) if is_pipeline else from_here(features_csv, row['Source'])

            # Duplicate Name Attributes and rule entries are not allowed
            if row['Name'] not in ["ID", *gff_files[gff]]: gff_files[gff].append(row['Name'])
            if rule not in rules: rules.append(rule)

    return rules, gff_files


class FeatureCounter:
    # Reference Tables
    features: HTSeq.GenomicArrayOfSets
    attributes: dict
    alias: dict

    selection_rules: List[dict]
    out_prefix: str
    run_diags: bool

    def __init__(self, gff_file_set, selection_rules, run_diags, out_prefix):
        reference_tables = parser.build_reference_tables(gff_file_set, selection_rules)
        FeatureCounter.features = reference_tables[0]
        FeatureCounter.attributes = reference_tables[1]
        FeatureCounter.alias = reference_tables[2]

        FeatureCounter.selection_rules = selection_rules
        FeatureCounter.out_prefix = out_prefix
        FeatureCounter.run_diags = run_diags

        self.chrom_misses = set()
        self.selector: FeatureSelector

    def assign_features(self, alignment: 'parser.Alignment') -> Tuple[AssignedFeatures, N_Candidates]:
        """Determines features associated with the interval then performs rule-based feature selection"""

        feat_matches, assignment = list(), set()
        iv = alignment.iv

        try:
            # Resolve features from alignment interval on both strands, regardless of alignment strand
            feat_matches = [match for match in
                            (self.features.chrom_vectors[iv.chrom]['.']  # GenomicArrayOfSets -> ChromVector
                             .array[iv.start:iv.end]                     # ChromVector -> StepVector
                             .get_steps(merge_steps=True))               # StepVector -> (iv_start, iv_end, {features})]
                            # If an alignment does not map to a feature, an empty set is returned at tuple position 2
                            if len(match[2]) != 0]
        except KeyError as ke:
            self.chrom_misses.add(ke.args[0])

        # If features are associated with the alignment interval, perform selection
        if len(feat_matches):
            assignment = self.selector.choose(feat_matches, alignment)

        return assignment, len(feat_matches)

    def count_reads(self, library: dict, return_queue: mp.Queue):
        """Collects statistics on features assigned to each alignment associated with each read"""

        # For complete SAM records (slower):
        # 1. Change the following line to HTSeq.BAM_Reader(sam_file)
        # 2. Change FeatureSelector.choose() to assign nt5end from chr(alignment.read.seq[0])
        read_seq = parser.read_SAM(library["File"])
        stats = LibraryStats(library, self.out_prefix, self.run_diags)
        self.selector = FeatureSelector(self.selection_rules, self.attributes, stats, self.run_diags)

        # For each sequence in the sam file...
        # Note: HTSeq only performs bundling. The alignments are our own Alignment objects
        for bundle in HTSeq.bundle_multiple_alignments(read_seq):
            bundle_stats = stats.count_bundle(bundle)

            # For each alignment of the given sequence...
            alignment: parser.Alignment
            for alignment in bundle:
                hits, n_candidates = self.assign_features(alignment)
                stats.count_bundle_alignments(bundle_stats, alignment, hits, n_candidates)

            stats.finalize_bundle(bundle_stats)

        # Place results in the multiprocessing queue to be merged by parent process
        return_queue.put(stats)

        # While stats are being merged, write intermediate file
        if run_diags:
            stats.write_intermediate_file()


@report_execution_time("Counter's overall runtime")
def main():
    # Get command line arguments.
    args = get_args()

    try:
        # Determine SAM inputs and their associated library names
        libraries = load_samples(args.input_csv)

        # Load selection rules and feature sources from config
        selection_rules, gff_file_set = load_config(args.config)

        # global for multiprocessing
        global feature_counter
        feature_counter = FeatureCounter(gff_file_set, selection_rules, args.run_diags, args.out_prefix)

        # Prepare for multiprocessing pool
        ret_queue = mp.Manager().Queue() if len(libraries) > 1 else queue.Queue()
        work_args = [(assignment, ret_queue) for assignment in libraries]

        # Assign and count features using multiprocessing and merge results
        merged_counts = map_and_reduce(libraries, work_args, ret_queue)

        # Write final outputs
        merged_counts.write_report_files(feature_counter.alias, args.out_prefix)
    except:
        traceback.print_exception(*sys.exc_info())
        print("\n\nCounter encountered an error. Don't worry! You don't have to start over.\n"
              "You can resume the pipeline at Counter. To do so:\n\t"
              "1. cd into your Run Directory\n\t"
              '2. Run "tiny recount --config your_run_config.yml"\n\t'
              '   (that\'s the processed run config) ^^^\n\n', file=sys.stderr)


if __name__ == '__main__':
    main()
