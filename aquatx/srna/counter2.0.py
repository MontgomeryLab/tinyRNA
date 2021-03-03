import itertools
import multiprocessing
import argparse
import operator
import re
import sys
import csv
from collections import Counter
from typing import Tuple, List, Set

import HTSeq

_re_attr_main = re.compile("\s*([^\s=]+)[\s=]+(.*)")
_re_attr_empty = re.compile("^\s*$")
features: 'HTSeq.GenomicArrayOfSets'
attributes = {}

class FeatureSelector:
    def __init__(self, selector):
        pass



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
    parser.add_argument('-o', '--out-prefix', metavar='OUTPUTPREFIX',
                        help='output prefix to use for file names')
    parser.add_argument('-t', '--intermed-file', action='store_true',
                        help='Save the intermediate file containing all alignments and'
                             'associated features.')
    parser.add_argument('-s', '--shared-memory', action='store_true',
                        help='')

    args = parser.parse_args()

    return args


def load_config(file: str) -> Tuple[Set[str], List[dict]]:
    gff_files = set()
    selectors = []

    short = {
        'Strand (sense/antisense/both)': 'Strand',
        'Feature Source': 'Source',
        "5' End Nucleotide": '5pnt'
    }

    with open(file, 'r', encoding='utf-8-sig') as f:
        csv_reader = csv.DictReader(f, delimiter=',')
        for row in csv_reader:
            selector = {short.get(col, col): row[col] for col in row}
            gff_files.add(selector['Source'])
            selectors.append(selector)

    return gff_files, selectors


def parse_GFF_attribute_string(attrStr, extra_return_first_value = False):
    """Parses a GFF attribute string and returns it as a dictionary.

    This is a slight modification of the same method found in HTSeq.features.
    It has been adapted to allow features to have multiple classes, which are
    stored as a list rather than a comma separated string. This should save
    some CPU time down the road.

    If 'extra_return_first_value' is set, a pair is returned: the dictionary
    and the value of the first attribute. This might be useful if this is the
    ID.
    """
    attribute_dict = {}
    first_val = "_unnamed_"
    for i, attr in enumerate(HTSeq._HTSeq.quotesafe_split(attrStr.rstrip('\n').encode())):
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
        attribute_dict[sys.intern(idt)] = sys.intern(val) \
            if idt != 'Class' \
            else [c.strip() for c in val.split(',')]
        if extra_return_first_value and i == 0:
            first_val = val
    if extra_return_first_value:
        return attribute_dict, first_val
    else:
        return attribute_dict


def parse_unified_gff(gff_files: set) -> Tuple[dict, dict]:
    # Some features have multiple classes. Need to patch the GFF attribute parser to parse them into lists
    setattr(HTSeq.features, 'parse_GFF_attribute_string', parse_GFF_attribute_string)

    # HTSeq's GenomicArrays cannot be merged, but we want a unified set
    # So instead...trick GFF_Reader into thinking there's just one file
    open_files = [open(gff, 'r') for gff in gff_files]
    file_chain = itertools.chain.from_iterable(open_files)

    # Open the "unified" GFF file and process features & attributes
    gff = HTSeq.GFF_Reader(file_chain)
    feature_scan = HTSeq.make_feature_genomicarrayofsets(
        gff, 'ID', additional_attributes=['Class', 'biotype'], stranded=True,
    )

    # Ensure entire file was read, then close it
    for handle in open_files:
        assert handle.tell() == handle.seek(0, 2)  # index 0 relative to end (2)
        handle.close()

    # Global to ensure that per-sam-file processes have a copy at launch
    global features, attributes
    features = feature_scan['features']
    attributes = feature_scan['attributes']

    return features, attributes


def assign_features(iv_seq, selectors) -> Tuple[list, list, int]:
    seq_feat = set()
    empty = 0

    # Build set of matching features based on match intervals
    for iv in iv_seq:
        # Check that our features table even contains the interval in question
        if iv.chrom not in features.chrom_vectors:
            empty += 1
        for iv2, fs2 in features[iv].steps():
            seq_feat = seq_feat.union(fs2)

    if len(seq_feat) == 0:
        pass
    else:
        [ident for ident in features if ident in selectors['Identity']]

        if selectors['Identity'] in seq_feat:
            if seq_feat[selectors['Identity']] == selectors['Feature']:



    return [],[],empty


def count_reads(sam_file, selectors):

    # CIGAR match characters (including alignment match, sequence match, and
    # sequence mismatch
    com = ('M', '=', 'X')
    counts = {feat: 0 for feat in sorted(attributes.keys())}
    ambiguous = notaligned = nonunique = 0
    read_seq = HTSeq.BAM_Reader(sam_file)

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
            # Obtain all matching intervals between the read and genome
            # This ought to be a single short interval for our purposes, but just in case...
            iv_seq = (co.ref_iv for co in alignment.cigar if co.type in com and co.size > 0)

            # Then perform selective assignment
            selected_features, selected_classes, empty = assign_features(iv_seq, selectors)
            stats_counts['_no_feature'] += empty * cor_counts

            if len(selected_classes) > 1:
                bundle_class["ambiguous"] += cor_counts
            if selected_features is None or len(selected_features) == 0:  # No feature
                stats_counts['_no_feature'] += cor_counts
            else:
                for feat in selected_features:
                    bundle_feats[feat] += cor_counts / len(selected_features)

            # Finally, count
            for fsi in list(selected_features):
                counts[fsi] += 1
                print(fsi)

        if len(bundle_class) != 0:
            class_counts["ambiguous"] += sum(bundle_class.values())
            stats_counts['_ambiguous_alignments_classes'] += 1
            stats_counts['_ambiguous_reads_classes'] += dup_counts

    # print(counts)
    print(f"empty {empty} ambiguous {ambiguous} notaligned {notaligned} nonunique {nonunique}")


def main():
    # Get command line arguments.
    args = get_args()

    # Load selection configuration
    gff_file_set, selectors = load_config(args.config)

    # Build features table, global for multiprocessing
    parse_unified_gff(gff_file_set)

    # Prepare for multiprocessing pool
    sam_files = [sam.strip() for sam in args.input_files.split(',')]
    pool_args = [[sam, selectors] for sam in sam_files]

    # Perform counts
    if len(sam_files) > 1:
        with multiprocessing.Pool(len(sam_files)) as pool:
            results = pool.starmap(count_reads, pool_args)
        results.sort(key=operator.itemgetter('sam_file'))
    else:
        results = list(itertools.starmap(count_reads, pool_args))

    print(results)


if __name__ == '__main__':
    main()