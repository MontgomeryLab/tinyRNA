import multiprocessing as mp
import argparse
import HTSeq
import queue
import time
import csv
import os

import aquatx.srna.hts_parsing as parser
from aquatx.srna.FeatureSelector import FeatureSelector
from aquatx.srna.statistics import LibraryStats, SummaryStats
from aquatx.srna.hts_parsing import SelectionRules, FeatureSources
from typing import Tuple

# Global variables for copy-on-write multiprocessing
features: 'HTSeq.GenomicArrayOfSets' = HTSeq.GenomicArrayOfSets("auto", stranded=True)
selector: 'FeatureSelector' = None

# Type aliases for human readability
AssignedFeatures = set
N_Candidates = int

# Todo: better way of handling Library titles. Filename or Report Title from run config?
#  input-file + out prefix combo? "-i Lib303.sam N2_rep_1, Lib311.sam mut-16(pk710)"
#  just pass in samples.csv? "-i samples.csv" however filenames will have changed by pipeline...
def get_args():
    """
    Get input arguments from the user/command line.

    Requires the user provide a SAM file, a GFF file, and an output
    prefix to save output tables and plots using.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-csv', metavar='SAMPLES', required=True,
                        help='the csv samples file/library list')
    parser.add_argument('-c', '--config', metavar='CONFIGFILE', required=True,
                        help='the csv features configuration file')
    parser.add_argument('-o', '--out-prefix', metavar='OUTPUTPREFIX', required=True,
                        help='output prefix to use for file names')
    parser.add_argument('-t', '--intermed-file', action='store_true',
                        help='Save the intermediate file containing all alignments and'
                             'associated features.')

    return parser.parse_args()


def load_samples(samples_csv: str):
    def get_input_filename(csv_row_file):
        """The input samples.csv may contain either fastq or sam files"""

        sample_file_ext = os.path.splitext(csv_row_file)[1].lower()

        if sample_file_ext == ".fastq":
            # Note: these fastq files MUST be included as an input (but not argument) for the CWL step
            return os.path.splitext(csv_row_file)[0] + "_aligned_seqs.sam"
        elif sample_file_ext == ".sam":
            if not os.path.isabs(csv_row_file):
                raise ValueError("The following file must be expressed as an absolute path:\n%s" % (csv_row_file,))
            return csv_row_file
        else:
            raise ValueError("The filenames defined in your samples CSV file must have a .fastq or .sam extension.\n"
                             "The following filename contained neither:\n%s" % (csv_row_file,))

    inputs = list()

    with open(samples_csv, 'r', encoding='utf-8-sig') as f:
        fieldnames = ("File", "Group", "Replicate")
        csv_reader = csv.DictReader(f, fieldnames=fieldnames, delimiter=',')

        next(csv_reader)  # Skip header line
        for row in csv_reader:
            library_name = f"{row['Group']}_replicate_{row['Replicate']}"
            library_file_name = get_input_filename(row['File'])
            record = {"Name": library_name, "File": library_file_name}

            if record not in inputs: inputs.append(record)

    return inputs


def load_config(features_csv: str) -> Tuple[SelectionRules, FeatureSources]:
    """Parses features.csv to provide inputs to FeatureSelector and build_reference_tables

    Args:
        features_csv: a csv file which defines feature sources and selection rules

    Returns:
        rules: a list of dictionaries, each representing a parsed row from input
        gff_files: a set of unique GFF files and associated ID attribute preferences
    """

    gff_files, rules = set(), list()
    convert_strand = {'sense': tuple('+'), 'antisense': tuple('-'), 'both': ('+', '-')}

    with open(features_csv, 'r', encoding='utf-8-sig') as f:
        fieldnames = ("ID", "Key", "Value", "Hierarchy", "Strand", "nt5", "Length", "Strict", "Source")
        csv_reader = csv.DictReader(f, fieldnames=fieldnames, delimiter=',')

        next(csv_reader)  # Skip header line
        for row in csv_reader:
            rule = {col: row[col] for col in ["Strand", "Hierarchy", "nt5", "Length", "Strict"]}
            rule['nt5'] = rule['nt5'].upper().translate({ord('U'): 'T'})  # Convert RNA base to cDNA base
            rule['Strand'] = convert_strand[rule['Strand'].lower()]       # Convert sense/antisense to +/-
            rule['Identity'] = (row['Key'], row['Value'])                 # Create identity tuple
            rule['Hierarchy'] = int(rule['Hierarchy'])                    # Convert hierarchy to number
            rule['Strict'] = rule['Strict'] == 'Full'                     # Convert strict intersection to boolean

            # Duplicate rule entries are not allowed
            if rule not in rules: rules.append(rule)
            gff_files.add((row['Source'], row['ID']))

    return rules, gff_files


def assign_features(alignment: 'parser.Alignment') -> Tuple[AssignedFeatures, N_Candidates]:

    feat_matches, assignment = list(), set()
    iv = alignment.iv

    feat_matches = [match for match in
                    (features.chrom_vectors[iv.chrom][iv.strand]  # GenomicArrayOfSets -> ChromVector
                             .array[iv.start:iv.end]              # ChromVector -> StepVector
                             .get_steps(merge_steps=True))        # StepVector -> (iv_start, iv_end, {features})]
                    # If an alignment does not map to a feature, an empty set is returned at tuple position 2
                    if len(match[2]) != 0]

    # If features are associated with the alignment interval, perform selection
    if len(feat_matches):
        assignment = selector.choose(feat_matches, alignment)

    return assignment, len(feat_matches)


def count_reads(library: dict, return_queue: mp.Queue, intermediate_file: bool = False, out_prefix: str = None):

    # For complete SAM records (slower):
    # 1. Change the following line to HTSeq.BAM_Reader(sam_file)
    # 2. Change FeatureSelector.choose() to assign nt5end from chr(alignment.read.seq[0])
    read_seq = parser.read_SAM(library["File"])
    stats = LibraryStats(library, out_prefix, intermediate_file)

    # For each sequence in the sam file...
    # Note: HTSeq only performs bundling. The alignments are our own Alignment objects
    for bundle in HTSeq.bundle_multiple_alignments(read_seq):
        bundle_stats = stats.count_bundle(bundle)

        # For each alignment of the given sequence...
        alignment: parser.Alignment
        for alignment in bundle:
            hits, n_candidates = assign_features(alignment)
            stats.count_bundle_alignments(bundle_stats, alignment, hits)

        stats.finalize_bundle(bundle_stats)

    # Place results in the multiprocessing queue to be merged by parent process
    return_queue.put(stats)

    if intermediate_file:
        stats.write_intermediate_file()


def map_and_reduce(libraries, work_args, ret_queue):
    mapred_start = time.time()

    # Use a multiprocessing pool if multiple sam files were provided
    # Otherwise perform counts in this process
    if len(libraries) > 1:
        with mp.Pool(len(libraries)) as pool:
            pool.starmap(count_reads, work_args)
    else:
        count_reads(*work_args[0])

    print("Counting took %.2f seconds" % (time.time() - mapred_start))
    merge_start = time.time()

    # Collect counts from all pool workers and merge
    summary = SummaryStats(libraries, work_args[0][-1])
    for _ in libraries:
        lib_stats = ret_queue.get()
        summary.add_library(lib_stats)

    print("Merging took %.2f seconds" % (time.time() - merge_start))
    return summary


def main():
    start = time.time()
    # Get command line arguments.
    args = get_args()

    # Determine SAM inputs and their associated library names
    libraries = load_samples(args.input_csv)

    # Load selection rules and feature sources from config
    selection_rules, gff_file_set = load_config(args.config)

    # Build features table and selector, global for multiprocessing
    global features, selector
    features, attributes, alias = parser.build_reference_tables(gff_file_set, selection_rules)
    selector = FeatureSelector(selection_rules, attributes)

    # Prepare for multiprocessing pool
    ret_queue = mp.Manager().Queue() if len(libraries) > 1 else queue.Queue()
    work_args = [(assignment, ret_queue, args.intermed_file, args.out_prefix) for assignment in libraries]

    # Assign and count features using multiprocessing and merge results
    merged_counts = map_and_reduce(libraries, work_args, ret_queue)

    # Write final outputs
    merged_counts.write_report_files(alias, args.out_prefix)
    print("Overall runtime took %.2f seconds" % (time.time() - start))


if __name__ == '__main__':
    main()
