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

def get_args():
    """
    Get input arguments from the user/command liobjectne.

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
    gff_files = set()
    rules = []

    with open(features_csv, 'r', encoding='utf-8-sig') as f:
        fieldnames = ("ID", "Key", "Value", "Strand", "Source", "Hierarchy", "5pnt", "Length")
        convert_strand = {'sense': tuple('+'), 'antisense': tuple('-'), 'both': ('+', '-')}
        csv_reader = csv.DictReader(f, fieldnames=fieldnames, delimiter=',')

        next(csv_reader)  # Skip header line
        for row in csv_reader:
            rule = {col: row[col] for col in ["Strand", "Hierarchy", "5pnt", "Length"]}
            rule['Hierarchy'] = int(rule['Hierarchy'])
            rule['Identity'] = (row['Key'], row['Value'])
            rule['Strand'] = convert_strand[rule['Strand'].lower()]
            rules.append(rule)

            gff_files.add((row['Source'], row['ID']))

    return rules, gff_files


def build_reference_tables(gff_files: FeatureSources, rules: SelectionRules) -> Tuple[Features, Attributes]:
    """A simplified and slightly modified version of HTSeq.create_genomicarrayofsets

    This modification changes the information stored in an interval's step vector
    within the features table. This stores (feature ID, feature type) tuples rather
    than the feature ID alone. It also allows for cataloguing features by any attribute,
    not just by ID, on a per GFF file basis.

    Note: at this time if the same feature is defined in multiple GFF files using
    different ID attributes, the feature will
    """

    b4_4 = time.time()  # Cloud Atlas

    # Patch the GFF attribute parser to support comma separated attribute value lists
    setattr(HTSeq.features, 'parse_GFF_attribute_string', parse_GFF_attribute_string)

    # Obtain an ordered list of unique attributes of interest from selection rules
    attrs_of_interest = list(np.unique([attr['Identity'][0] for attr in rules]))

    feats = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    attrs = {}

    for file, preferred_id in gff_files:
        gff = HTSeq.GFF_Reader(file)
        for f in gff:
            try:
                if f.iv.strand == ".":
                    raise ValueError("Feature %s at %s in %s does not have strand information." % (f.name, f.iv, file))

                feature_id = f.attr[preferred_id]
                feats[f.iv] += (feature_id, f.type)
                attrs[feature_id] = [(attr, f.attr[attr]) for attr in attrs_of_interest]
            except KeyError as ke:
                raise ValueError("Feature %s does not contain a '%s' attribute in %s" % (f.name, ke, file))

    print("GFF parsing took %.2f seconds" % (time.time() - b4_4))
    return feats, attrs


def assign_features(alignment: 'HTSeq.SAM_Alignment') -> Tuple[AssignedFeatures, N_Candidates]:
    # CIGAR match characters
    com = ('M', '=', 'X')

    feature_set, assignment = set(), set()

    iv = alignment.iv
    for fs2 in features[iv].array[iv.start:iv.end].get_steps(values_only=True):  # Not as pretty, but much faster!
        feature_set |= fs2

    if len(feature_set):
        strand = alignment.iv.strand
        nt5end = chr(alignment.read.seq[0])
        length = len(alignment.read)
        choices, uncounted = selector.choose(feature_set, strand, nt5end, length)
        assignment |= choices
    else:
        assignment.add(selector.empty)

    return assignment, len(feature_set)


def count_reads(sam_file: str, return_queue: mp.Queue, intermediate_file: bool = False, out_prefix: str = None):

    # Change the following line to HTSeq.BAM_Reader(sam_file) for complete SAM records (slower)
    read_seq = read_SAM(sam_file)
    stats = StatsCollector(out_prefix, intermediate_file)

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
        im_prefix = os.path.splitext(os.path.basename(sam_file))[0]
        stats.write_intermediate_file(im_prefix)


def map_and_reduce(sam_files, work_args, ret_queue):
    b4_4 = time.time()

    # Use a multiprocessing pool if multiple sam files were provided
    # Otherwise perform counts in this process
    if len(sam_files) > 1:
        with mp.Pool(len(sam_files)) as pool:
            pool.starmap(count_reads, work_args)
    else:
        count_reads(*work_args[0])

    # Collect counts from all pool workers and merge
    merged_counts = StatsCollector()
    for _ in sam_files:
        sam_stats = ret_queue.get()
        merged_counts.merge(sam_stats)

    print("Counting and merging took %.2f seconds" % (time.time() - b4_4))
    return merged_counts


def main():
    start = time.time()
    # Get command line arguments.
    args = get_args()

    # Load selection rules and feature sources from config
    selection_rules, gff_file_set = load_config(args.config)

    # Build features table and selector, global for multiprocessing
    global features, selector
    features, attributes = build_reference_tables(gff_file_set, selection_rules)
    selector = FeatureSelector(selection_rules, attributes)

    # Prepare for multiprocessing pool
    sam_files = [sam.strip() for sam in args.input_files.split(',')]
    ret_queue = mp.Manager().Queue() if len(sam_files) > 1 else queue.Queue()
    work_args = [(sam, ret_queue, args.intermed_file, args.out_prefix) for sam in sam_files]

    # Assign and count features using multiprocessing and merge results
    merged_counts = map_and_reduce(sam_files, work_args, ret_queue)

    # Write final outputs
    merged_counts.write_report_files(args.out_prefix)
    print("Overall runtime took %.2f seconds" % (time.time() - start))


if __name__ == '__main__':
    main()
