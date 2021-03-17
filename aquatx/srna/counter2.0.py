import itertools
import argparse
import time

import operator
import numpy as np
import HTSeq
import sys
import csv
import re

from aquatx.srna.FeatureSelector import FeatureSelector
from multiprocessing import Pool, Queue
from collections import Counter, OrderedDict
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
    parser.add_argument('-o', '--out-prefix', metavar='OUTPUTPREFIX',
                        help='output prefix to use for file names')
    parser.add_argument('-t', '--intermed-file', action='store_true',
                        help='Save the intermediate file containing all alignments and'
                             'associated features.')

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
            assignment.add('Empty')
            continue
        for iv2, fs2 in features[iv].steps():
            feature_set |= fs2

    if len(feature_set):
        strand = alignment.iv.strand
        nt5end = chr(alignment.read.seq[0])
        length = len(alignment.read)
        assignment |= selector.choose(feature_set, strand, nt5end, length)

    return assignment, len(feature_set)


def count_reads(sam_file: str):
    read_seq = HTSeq.BAM_Reader(sam_file)

    counts = {feat: 0 for feat in sorted(selector.attributes.keys(), key=lambda x: x[0])}
    nt_len_mat = {nt: Counter() for nt in ['A', 'T', 'G', 'C']}
    stats_counts = {stat: 0 for stat in
                    ['_aligned_reads', '_aligned_reads_unique_mapping', '_aligned_reads_multi_mapping',
                     '_unique_sequences_aligned', '_reads_unique_features', '_alignments_unique_features',
                     '_ambiguous_alignments_classes', '_ambiguous_reads_classes', '_ambiguous_alignments_features',
                     '_ambiguous_reads_features', '_no_feature']}

    for bundle in HTSeq.bundle_multiple_alignments(read_seq):
        # Calculate counts for multimapping
        dup_counts = int(bundle[0].read.name.split('=')[1])
        cor_counts = dup_counts / len(bundle)
        stats_counts['_unique_sequences_aligned'] += 1
        stats_counts['_aligned_reads'] += dup_counts
        if len(bundle) > 1:
            stats_counts['_aligned_reads_multi_mapping'] += dup_counts
        else:
            stats_counts['_aligned_reads_unique_mapping'] += dup_counts

        # fill in 5p nt/length matrix
        nt_len_mat[str(bundle[0].read)[0]][len(bundle[0].read)] += dup_counts

        # bundle counts
        bundle_feats = Counter()
        bundle_class = Counter()

        alignment: HTSeq.SAM_Alignment
        for alignment in bundle:
            # Perform selective assignment
            selected_features, total_feat = assign_features(alignment)

            if len({x[selector.feat] for x in selected_features}) > 1:
                bundle_class["ambiguous"] += cor_counts
            if selected_features is None or len(selected_features) == 0:  # No feature
                stats_counts['_no_feature'] += cor_counts
            else:
                for feat in selected_features:
                    bundle_feats[feat] += cor_counts / total_feat

            # Finally, count
            for fsi in selected_features:
                counts[fsi[selector.feat]] += 1

    # Todo: return via multiprocessing queue
    #  Then assemble/write intermediate outputs here while the main process merges counts from queue
    return {
        'counts': counts,
        'sam_file': sam_file
    }

def disperse_and_aggregate():
    pass


def main():
    start = time.time()
    # Get command line arguments.
    args = get_args()

    # Load config csv, selector stored globally
    gff_file_set = load_config(args.config)

    # Build features table, global for multiprocessing
    b4 = time.time()
    build_reference_tables(gff_file_set)
    print("GFF parsing took %.2f seconds" % (time.time() - b4))

    # Prepare for multiprocessing pool
    sam_files = [sam.strip() for sam in args.input_files.split(',')]
    pool_args = [[sam] for sam in sam_files]

    # Perform counts
    b4 = time.time()
    if len(sam_files) > 1:
        with Pool(len(sam_files)) as pool:
            results = pool.starmap(count_reads, pool_args)
        results.sort(key=operator.itemgetter('sam_file'))
    else:
        results = list(itertools.starmap(count_reads, pool_args))
    print("Counting took %.2f seconds" % (time.time() - b4))

    out = {}
    for file in results:
        print(file['sam_file'])
        for k,v in file['counts'].items():
            if v:
                out[k] = v
        print(out)
    print("Overall runtime took %.2f seconds" % (time.time() - start))


if __name__ == '__main__':
    main()
