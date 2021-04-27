import argparse
import queue
import time
import os

import multiprocessing as mp
import numpy as np
import csv

from aquatx.srna.FeatureSelector import FeatureSelector, StatsCollector
from aquatx.srna.hts_parsing import *
from typing import Tuple, Set, List, Dict

# Global variables for copy-on-write multiprocessing
features: 'HTSeq.GenomicArrayOfSets' = HTSeq.GenomicArrayOfSets("auto", stranded=True)
selector: 'FeatureSelector' = None

# Type aliases for human readability
Features = HTSeq.GenomicArrayOfSets  # interval -> set of associated features
Attributes = Dict[str, list]  # feature -> feature attributes
FeatureSources = Set[Tuple[str, str]]
SelectionRules = List[dict]
AssignedFeatures = set
N_Candidates = int
Alias = dict

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
    parser.add_argument('-i', '--input-files', metavar='SAMFILES', required=True,
                        help='comma separated list of input sam files to count features for')
    parser.add_argument('-c', '--config', metavar='CONFIGFILE', required=True,
                        help='the csv features configuration file')
    parser.add_argument('-o', '--out-prefix', metavar='OUTPUTPREFIX', required=True,
                        help='output prefix to use for file names')
    parser.add_argument('-t', '--intermed-file', action='store_true',
                        help='Save the intermediate file containing all alignments and'
                             'associated features.')

    return parser.parse_args()


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
            rule['nt5'] = rule['nt5'].upper().translate({'U': 'T'})   # Convert RNA base to cDNA base
            rule['Strand'] = convert_strand[rule['Strand'].lower()]   # Convert sense/antisense to +/-
            rule['Identity'] = (row['Key'], row['Value'])             # Create identity tuple
            rule['Hierarchy'] = int(rule['Hierarchy'])                # Convert hierarchy to number
            rule['Strict'] = rule['Strict'] == 'Full'                 # Convert strict intersection to boolean

            # Duplicate rule entries are not allowed
            if rule not in rules: rules.append(rule)
            gff_files.add((row['Source'], row['ID']))

    return rules, gff_files


def build_reference_tables(gff_files: FeatureSources, rules: SelectionRules) -> Tuple[Features, Attributes, Alias]:
    """A simplified and slightly modified version of HTSeq.create_genomicarrayofsets

    This modification changes the information stored in an interval's step vector
    within the features table. This stores (feature ID, feature type) tuples rather
    than the feature ID alone. It also allows for cataloguing features by any attribute,
    not just by ID, on a per GFF file basis.

    Note: at this time if the same feature is defined in multiple GFF files using
    different ID attributes, the feature will
    """

    start_time = time.time()

    # Patch the GFF attribute parser to support comma separated attribute value lists
    setattr(HTSeq.features, 'parse_GFF_attribute_string', parse_GFF_attribute_string)

    # Obtain an ordered list of unique attributes of interest from selection rules
    attrs_of_interest = list(np.unique([attr['Identity'][0] for attr in rules]))

    feats = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    attrs, alias = {}, {}

    for file, preferred_id in gff_files:
        gff = HTSeq.GFF_Reader(file)
        for row in gff:
            if row.iv.strand == ".":
                raise ValueError(f"Feature {row.name} at {row.iv} in {row.file} has no strand information.")

            try:
                # Add feature_id -> feature_interval record
                feature_id = row.attr["ID"][0]
                feats[row.iv] += feature_id
                row_attrs = [(interest, row.attr[interest]) for interest in attrs_of_interest]

                if preferred_id != "ID":
                    # If an alias already exists for this feature, append to feature's aliases
                    alias[feature_id] = alias.get(feature_id, ()) + row.attr[preferred_id]
            except KeyError as ke:
                raise ValueError(f"Feature {row.name} does not contain a {ke} attribute in {file}")

            if feature_id in attrs and row_attrs != attrs[feature_id]:
                # If an attribute record already exists for this feature, and this row provides new attributes,
                #  append the new attribute values to the existing values
                cur_attrs = attrs[feature_id]
                row_attrs = [(cur[0], cur[1] + new[1]) for cur, new in zip(cur_attrs, row_attrs)]

            # Add feature_id -> feature_alias_tuple record
            attrs[feature_id] = row_attrs

    print("GFF parsing took %.2f seconds" % (time.time() - start_time))
    return feats, attrs, alias


def assign_features(alignment: 'HTSeq.SAM_Alignment') -> Tuple[AssignedFeatures, N_Candidates]:

    feat_matches, assignment = list(), set()
    iv = alignment.iv

    for match_tuple in (features.chrom_vectors[iv.chrom][iv.strand]  # GenomicArrayOfSets -> ChromVector
                                .array[iv.start:iv.end]              # ChromVector -> StepVector
                                .get_steps(merge_steps=True)):       # StepVector -> (iv_start, iv_end, {features})
        if len(match_tuple[2]):
            feat_matches.append(match_tuple)

    if len(feat_matches):
        assignment, uncounted = selector.choose(feat_matches, alignment)

    return assignment, len(feat_matches)


def get_library_name(sam_file):
    filename = os.path.splitext(os.path.basename(sam_file))[0]
    return filename.replace("_aligned_seqs", "")


def count_reads(sam_file: str, return_queue: mp.Queue, intermediate_file: bool = False, out_prefix: str = None):

    # For complete SAM records (slower):
    # 1. Change the following line to HTSeq.BAM_Reader(sam_file)
    # 2. Change FeatureSelector.choose() to assign nt5end from chr(alignment.read.seq[0])
    read_seq = read_SAM(sam_file)
    lib_name = get_library_name(sam_file)
    stats = StatsCollector(lib_name, out_prefix, intermediate_file)

    # For each sequence in the sam file...
    for bundle in HTSeq.bundle_multiple_alignments(read_seq):
        bundle_stat = stats.count_bundle(bundle)

        # For each alignment of the given sequence...
        for alignment in bundle:
            hits, n_candidates = assign_features(alignment)
            if n_candidates > 0 and len(hits):
                stats.count_bundle_alignments(bundle_stat, alignment, hits)

        stats.finalize_bundle(bundle_stat)
    # Place results in the multiprocessing queue to be merged by parent process
    # We only want to return pertinent/minimal stats since this will be pickled
    return_queue.put(stats)

    if intermediate_file:
        stats.write_intermediate_file()


def map_and_reduce(sam_files, work_args, ret_queue):
    mapred_start = time.time()

    # Use a multiprocessing pool if multiple sam files were provided
    # Otherwise perform counts in this process
    if len(sam_files) > 1:
        with mp.Pool(len(sam_files)) as pool:
            pool.starmap(count_reads, work_args)
    else:
        count_reads(*work_args.pop())

    print("Counting took %.2f seconds" % (time.time() - mapred_start))
    merge_start = time.time()

    # Collect counts from all pool workers and merge
    merged_counts = StatsCollector("Summary")
    for _ in sam_files:
        sam_stats = ret_queue.get()
        merged_counts.merge(sam_stats)

    print("Merging took %.2f seconds" % (time.time() - merge_start))
    return merged_counts


def main():
    start = time.time()
    # Get command line arguments.
    args = get_args()

    # Load selection rules and feature sources from config
    selection_rules, gff_file_set = load_config(args.config)

    # Build features table and selector, global for multiprocessing
    global features, selector
    features, attributes, alias = build_reference_tables(gff_file_set, selection_rules)
    selector = FeatureSelector(selection_rules, attributes)

    # Prepare for multiprocessing pool
    sam_files = [sam.strip() for sam in args.input_files.split(',')]
    ret_queue = mp.Manager().Queue() if len(sam_files) > 1 else queue.Queue()
    work_args = [(sam, ret_queue, args.intermed_file, args.out_prefix) for sam in sam_files]

    # Assign and count features using multiprocessing and merge results
    merged_counts = map_and_reduce(sam_files, work_args, ret_queue)

    # Write final outputs
    merged_counts.write_report_files(alias, args.out_prefix)
    print("Overall runtime took %.2f seconds" % (time.time() - start))


if __name__ == '__main__':
    main()
