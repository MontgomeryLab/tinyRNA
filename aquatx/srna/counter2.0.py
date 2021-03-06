import itertools
import multiprocessing
import argparse
import operator
import re
import sys
import csv
from collections import Counter, defaultdict
from operator import itemgetter
from typing import Tuple, List, Set

import HTSeq

# For parse_GFF_attribute_string()
_re_attr_main = re.compile("\s*([^\s=]+)[\s=]+(.*)")
_re_attr_empty = re.compile("^\s*$")

# Global variables for copy-on-write multiprocessing
features: 'HTSeq.GenomicArrayOfSets' = None
selector: 'Selector' = None
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
    parser.add_argument('-s', '--shared-memory', action='store_true',
                        help='')

    args = parser.parse_args()

    return args


def load_config(file: str) -> Set[str]:
    gff_files = set()
    rules = []

    with open(file, 'r', encoding='utf-8-sig') as f:
        fieldnames = ("Identifier", "Feature", "Strand", "Source", "Hierarchy", "5pnt", "Length")
        csv_reader = csv.DictReader(f, fieldnames=fieldnames, delimiter=',')
        next(csv_reader)  # Skip header line
        for row in csv_reader:
            rule = {col: row[col] for col in row if col not in ["Identifier", "Feature", "Source"]}
            rule['Identity'] = (row['Identifier'], row['Feature'])
            rules.append(rule)
            gff_files.add(row['Source'])

    for rule in rules:
        rule['Hierarchy'] = int(rule['Hierarchy'])

    global selector
    selector = Selector(rules)
    return gff_files


def parse_GFF_attribute_string(attrStr, extra_return_first_value=False):
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
        gff, 'ID', additional_attributes=[k[0] for k in Selector.ident_idx], stranded=True
    )

    # Ensure entire file was read, then close it
    for handle in open_files:
        assert handle.tell() == handle.seek(0, 2)  # index 0 relative to end (2)
        handle.close()

    # Global to ensure that multi-processing workers have a reference
    global features, attributes
    features = feature_scan['features']
    attributes = feature_scan['attributes']

    return features, attributes


class Selector:
    # Keys determine which attributes are extracted when parsing GFF files (case must match GFF)
    # Values indicate corresponding index within each feature in the resulting attributes[] table
    ident_idx = [('Class', 0), ('biotype', 1)]

    # filter types: tuple membership, range, list, wildcard
    # Identifier/Class: need to lookup in attributes table
    # Strand, 5pnt, Length: provided by assign_features()

    def __init__(self, rules: List[dict]):
        self.interest = ['Identity', 'Strand', '5pnt', 'Length']
        # The ruleset provided, sorted by hierarchy
        self.rule_def = sorted(rules, key=lambda x: int(x['Hierarchy']))
        # For each interest, create dict of preference: [associated rule indexes]
        self.int_sets = []
        for step in self.interest:
            inverted_rules = defaultdict(list)
            for index, rule in enumerate(self.rule_def):
                inverted_rules[rule[step]].append(index)
            self.int_sets.append(dict(inverted_rules))

    def identity_choice(self, feat_set, interests):
        # For each feature, match across all identity interests
        identity_hits = [(feat, interests[(id, attr)])
                         for feat in feat_set
                         for id, i in Selector.ident_idx
                         for attr in attributes[feat][i]
                         if (id, attr) in interests]
        # -> [(feature, [matched rule indexes,...]), ...]

        choices, finalists = set(), set()

        if len(identity_hits) < len(feat_set):
            choices.add('Unknown')
        if len(identity_hits) == 1 and len(identity_hits[0][1]) == 1:
            # Only one feature matched only one rule
            finalists.add(identity_hits[0][0])
        else:
            # Perform any possible hierarchy-based eliminations
            ranks = [self.rule_def[rule]['Hierarchy'] for hit in identity_hits for rule in hit[1]]
            rank_set = set(ranks)
            # Flatten identity_hits and add hierarchy to each tuple
            get_rank = lambda x: x[0]
            hit_rank = [z for z in zip(ranks,
                        (f[0] for f in identity_hits for _ in f),
                        (r for i in identity_hits for r in i[1]))]
            # -> [(hierarchy, feature, rule), ...]

            if len(hit_rank) == len(rank_set):
                finalists.add(min(hit_rank, key=get_rank))
            else:
                # Two or more hits share the same hierarchy.
                min_rank = min(rank_set)
                finalists.add(hit for hit in hit_rank if get_rank(hit) == min_rank)

        if len(finalists) < len(identity_hits): choices.add('Unknown')
        return choices, finalists

    def is_complete(self, choices: set, finalists: set) -> bool:
        if len(finalists) == 1:
            choices.add(finalists.pop())
            return True
        else:
            return len(finalists) == 0

    def choose(self, feat_set, strand, endnt, length) -> set:
        step = iter(self.int_sets)

        # Identity
        identity_interests = next(step)
        choices, finalists = self.identity_choice(feat_set, identity_interests)
        if not finalists: return choices

        # Strand
        strand_interests = next(step)
        for feat in finalists:
            # Todo: convert features.csv input value to HTSeq's native representation
            if strand != self.rule_def[feat[2]]['Strand']:
                finalists.remove(feat)
                choices.add('Unknown')
        if not finalists: return choices


def assign_features(alignment) -> Tuple[list, list, int]:
    # CIGAR match characters
    com = ('M', '=', 'X')
    # Generator for all match intervals between the read and genome
    # This ought to be a single short interval for our purposes, but just in case...
    iv_seq = (co.ref_iv for co in alignment.cigar if co.type in com and co.size > 0)

    feature_set = set()
    empty = 0

    # Build set of matching features based on match intervals
    for iv in iv_seq:
        # Check that our features table even contains the interval in question
        if iv.chrom not in features.chrom_vectors:
            empty += 1
        for iv2, fs2 in features[iv].steps():
            feature_set = feature_set.union(fs2)

    if len(feature_set) == 0:
        pass
    else:
        strand = alignment.iv.strand
        nt5end = chr(alignment.read.seq[0])
        length = len(alignment.read)
        assignment = selector.choose(feature_set, strand, nt5end, length)

    return [], [], empty


def count_reads(sam_file):
    read_seq = HTSeq.BAM_Reader(sam_file)

    counts = {feat: 0 for feat in sorted(attributes.keys())}
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
            selected_features, selected_classes, empty = assign_features(alignment)
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

    print(counts)


def main():
    # Get command line arguments.
    args = get_args()

    # Load config csv, selector stored globally
    gff_file_set = load_config(args.config)

    # Build features table, global for multiprocessing
    parse_unified_gff(gff_file_set)

    # Prepare for multiprocessing pool
    sam_files = [sam.strip() for sam in args.input_files.split(',')]
    pool_args = [[sam] for sam in sam_files]

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
