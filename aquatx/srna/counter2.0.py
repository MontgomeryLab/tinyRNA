import itertools
import argparse
import os
import queue
import time

import multiprocessing as mp
import numpy as np
import HTSeq
import sys
import csv
import re

from aquatx.srna.FeatureSelector import FeatureSelector, StatsCollector
from typing import Tuple, Set

# For parse_GFF_attribute_string()
_re_attr_main = re.compile(r"\s*([^\s=]+)[\s=]+(.*)")
_re_attr_empty = re.compile(r"^\s*$")

# Global variables for copy-on-write multiprocessing
features: 'HTSeq.GenomicArrayOfSets' = HTSeq.GenomicArrayOfSets("auto", stranded=True)
selector: 'FeatureSelector' = None
attributes: dict = {}


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

    global args
    args = parser.parse_args()

    return args


def load_config(file: str) -> Set[Tuple[str, str]]:
    gff_files = set()
    rules = []

    with open(file, 'r', encoding='utf-8-sig') as f:
        fieldnames = ("ID", "Key", "Value", "Strand", "Source", "Hierarchy", "5pnt", "Length")
        csv_reader = csv.DictReader(f, fieldnames=fieldnames, delimiter=',')
        next(csv_reader)  # Skip header line
        for row in csv_reader:
            rule = {col: row[col] for col in row if col not in ["ID", "Key", "Value", "Source"]}
            rule['Identity'] = (row['Key'], row['Value'])
            rules.append(rule)
            gff_files.add((row['Source'], row['ID']))

    convert_strand = {'sense': tuple('+'), 'antisense': tuple('-'), 'both': ('+', '-')}
    for rule in rules:
        rule['Hierarchy'] = int(rule['Hierarchy'])
        rule['Strand'] = convert_strand[rule['Strand'].lower()]

    global selector
    selector = FeatureSelector(rules)
    return gff_files


def parse_GFF_attribute_string(attrStr, extra_return_first_value=False):
    """Parses a GFF attribute string and returns it as a dictionary.

    This is a slight modification of the same method found in HTSeq.features.
    It has been adapted to allow features to have multiple classes, which are
    split and stored as a tuple rather than a comma separated string. This
    should save some CPU time down the road.

    If 'extra_return_first_value' is set, a pair is returned: the dictionary
    and the value of the first attribute. This might be useful if this is the
    ID.
    """
    attribute_dict = {}
    first_val = "_unnamed_"
    for i, attr in enumerate(HTSeq._HTSeq.quotesafe_split(attrStr.rstrip().encode())):
        attr = attr.decode()
        if _re_attr_empty.match(attr):
            continue
        if attr.count('"') not in (0, 2):
            raise ValueError(
                "The attribute string seems to contain mismatched  quotes.")
        mo = _re_attr_main.match(attr)
        if not mo:
            raise ValueError("Failure parsing GFF attribute line")
        idt = mo.group(1)
        val = mo.group(2)
        if val.startswith('"') and val.endswith('"'):
            val = val[1:-1]
        # Modification: allow for comma separated multiple-class assignment
        attribute_dict[sys.intern(idt)] = (sys.intern(val),) \
            if ',' not in val \
            else tuple(c.strip() for c in val.split(','))
        if extra_return_first_value and i == 0:
            first_val = val
    if extra_return_first_value:
        return attribute_dict, first_val
    else:
        return attribute_dict


def build_reference_tables(gff_files: Set[Tuple[str, str]]) -> None:
    """A simplified and slightly modified version of HTSeq.create_genomicarrayofsets

    This modification changes the information stored in an interval's step vector
    within the features table. This stores (feature ID, feature type) tuples rather
    than the feature ID alone.
    """

    b4_4 = time.time()  # Cloud Atlas

    # Patch the GFF attribute parser to support comma separated attribute value lists
    setattr(HTSeq.features, 'parse_GFF_attribute_string', parse_GFF_attribute_string)

    # Obtain an ordered list of unique attributes of interest from selector's Identity rules
    attrs_of_interest = list(np.unique([attr['Identity'][0] for attr in selector.rules_table]))

    feats = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    attributes = {}

    for file, preferred_id in gff_files:
        gff = HTSeq.GFF_Reader(file)
        for f in gff:
            try:
                feature_id = f.attr[preferred_id]
                if f.iv.strand == ".":
                    raise ValueError("Feature %s at %s in %s does not have strand information." % (f.name, f.iv, file))

                feats[f.iv] += (feature_id, f.type)
                attributes[feature_id] = [(attr, f.attr.get(attr, '')) for attr in attrs_of_interest]
            except KeyError as ke:
                raise ValueError("Feature %s does not contain a '%s' attribute in %s" % (f.name, ke, file))

    global features
    features = feats
    selector.set_attributes_table(attributes, attrs_of_interest)

    print("GFF parsing took %.2f seconds" % (time.time() - b4_4))


def assign_features(alignment) -> Tuple[set, int]:
    # CIGAR match characters
    com = ('M', '=', 'X')
    # Generator for all match intervals between the read and genome
    # This ought to be a single short interval for our purposes, but just in case...
    iv_seq = (co.ref_iv for co in alignment.cigar if co.type in com and co.size > 0)

    feature_set, assignment = set(), set()

    # Build set of matching features based on match intervals
    for iv in iv_seq:
        # Check that features[] even contains this alignment's chromosome
        if iv.chrom not in features.chrom_vectors:
            assignment.add(selector.empty)
            continue
        for iv2, fs2 in features[iv].steps():
            feature_set |= fs2

    if len(feature_set):
        strand = alignment.iv.strand
        nt5end = chr(alignment.read.seq[0])
        length = len(alignment.read)
        assignment |= selector.choose(feature_set, strand, nt5end, length)
    else:
        assignment.add(selector.empty)

    return assignment, len(feature_set)


def count_reads(sam_file: str, return_queue: mp.Queue, intermediate_file: bool = False, out_prefix: str = None):

    read_seq = HTSeq.BAM_Reader(sam_file)
    stats = StatsCollector(out_prefix, intermediate_file)

    # For each sequence in the sam file...
    for bundle in HTSeq.bundle_multiple_alignments(read_seq):
        stats.count_bundle(bundle)

        # For each alignment of the given sequence...
        for alignment in bundle:
            hits, n_candidates = assign_features(alignment)
            if n_candidates > 0 and len(hits):
                stats.count_alignment(alignment, hits, n_candidates)

    stats.finalize_bundles()
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

    # Load config csv, selector stored globally
    gff_file_set = load_config(args.config)

    # Build features table, global for multiprocessing
    build_reference_tables(gff_file_set)

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
