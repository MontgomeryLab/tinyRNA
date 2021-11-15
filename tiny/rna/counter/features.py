import itertools
import HTSeq

from collections import defaultdict
from typing import List, Tuple, Set, Dict, Iterator

import tiny.rna.counter.hts_parsing as parser
from .matching import Wildcard, StrandMatch, NumericalMatch, NtMatch
from .statistics import LibraryStats

# Type aliases for human readability
Hit = Tuple[int, int, str]

BOTH_STRANDS = ('+', '-')

# Global indexes for Hits produced by choose_identity()
RANK, RULE, FEAT = 0, 1, 2


class Features:
    chrom_vectors: HTSeq.ChromVector
    attributes: dict
    intervals: dict
    aliases: dict

    _instance = None  # Singleton

    def __init__(self, features: HTSeq.GenomicArrayOfSets, attributes: dict, aliases: dict, intervals: dict):
        if Features._instance is None:
            Features.chrom_vectors = features.chrom_vectors  # For interval -> feature ID lookups
            Features.attributes = attributes                 # For feature ID -> GFF column 9 attribute lookups
            Features.aliases = aliases                       # For feature ID -> preferred feature name lookups
            Features.intervals = intervals                   # For feature ID -> interval lookups
            Features._instance = self


class FeatureCounter:
    out_prefix: str
    run_diags: bool

    def __init__(self, gff_file_set, selection_rules, **prefs):
        reference_tables = parser.ReferenceTables(gff_file_set, selection_rules, **prefs)
        Features(*reference_tables.get())

        FeatureCounter.out_prefix = prefs['out_prefix']
        FeatureCounter.run_diags = prefs['report_diags']

        self.stats = LibraryStats(**prefs)
        self.selector = FeatureSelector(selection_rules, self.stats, **prefs)

    def assign_features(self, alignment: 'parser.Alignment') -> Tuple[set, int]:
        """Determines features associated with the interval then performs rule-based feature selection"""

        feat_matches, assignment = list(), set()
        iv = alignment.iv

        try:
            # Resolve features from alignment interval on both strands, regardless of alignment strand
            feat_matches = [match for strand in BOTH_STRANDS for match in
                            (Features.chrom_vectors[iv.chrom][strand]  # GenomicArrayOfSets -> ChromVector
                                     .array[iv.start:iv.end]           # ChromVector -> StepVector
                                     .get_steps(merge_steps=True))     # StepVector -> (iv_start, iv_end, {features})
                            # If an alignment does not map to a feature, an empty set is returned at tuple pos 2 ^^^
                            if len(match[2]) != 0]
        except KeyError as ke:
            self.stats.chrom_misses[ke.args[0]] += 1

        # If features are associated with the alignment interval, perform selection
        if len(feat_matches):
            assignment = self.selector.choose(feat_matches, alignment)

        return assignment, len(feat_matches)

    def count_reads(self, library: dict):
        """Collects statistics on features assigned to each alignment associated with each read"""

        # For complete SAM records (slower):
        # 1. Change the following line to HTSeq.BAM_Reader(library["File"])
        # 2. Change FeatureSelector.choose() to assign nt5end from chr(alignment.read.seq[0])
        read_seq = parser.read_SAM(library["File"])
        self.stats.assign_library(library)

        # For each sequence in the sam file...
        # Note: HTSeq only performs bundling. The alignments are our own Alignment objects
        for bundle in HTSeq.bundle_multiple_alignments(read_seq):
            bundle_stats = self.stats.count_bundle(bundle)

            # For each alignment of the given sequence...
            alignment: parser.Alignment
            for alignment in bundle:
                hits, n_candidates = self.assign_features(alignment)
                self.stats.count_bundle_alignments(bundle_stats, alignment, hits, n_candidates)

            self.stats.finalize_bundle(bundle_stats)

        # While stats are being merged, write intermediate file
        if FeatureCounter.run_diags:
            self.stats.diags.write_intermediate_file(library["Name"])

        return self.stats


class FeatureSelector:
    """Performs hierarchical selection given a set of candidate features for a locus

    Two sources of data serve as targets for selection: feature attributes (sourced from
    input GFF files), and sequence attributes (sourced from input SAM files).

    The first round of selection is performed against each candidate feature's attributes.
    The target for this stage is feature attribute key-value pairs, referred to here as Identities.
    A candidate may match multiple identities. Each match is referred to as a Hit. If more
    than one Hit is produced, elimination is performed using each Hit's hierarchy/rank value.
    Hits are tuples for performance reasons, and are of the format:
        (hierarchy, rule, feature_id)

    If more than one Hit remains following first round selection, a second round of selection
    is performed against sequence attributes: strand, 5' end nucleotide, and length. Rules for
    5' end nucleotides support lists (e.g. C,G,U) and wildcards (e.g. "all"). Rules for length
    support lists, wildcards, and ranges (i.e. 20-27) which may be intermixed in the same rule.
    Lengths may be specified as "strict", meaning that the feature must be completely contained
    by the alignment interval.
    """

    rules_table = dict()
    inv_ident = dict()

    def __init__(self, rules: List[dict], libstats: 'LibraryStats', report_diags=False, **kwargs):
        FeatureSelector.rules_table = self.build_selectors(rules)
        FeatureSelector.inv_ident = self.build_inverted_identities(FeatureSelector.rules_table)

        self.report_eliminations = report_diags
        if report_diags: self.elim_stats = libstats.diags.selection_diags

    def choose(self, feat_list: List[Tuple[int, int, Set[str]]], alignment: 'parser.Alignment') -> Set[str]:
        # Perform hierarchy-based first round of selection for identities
        finalists = self.choose_identities(feat_list, alignment.iv)
        if not finalists: return set()

        read_aln_attrs = {
            'Strand': (alignment.iv.strand,),
            'nt5end': alignment.read.nt5,
            'Length': len(alignment.read)
        }

        eliminated = set()
        for selector, read in read_aln_attrs.items():
            for hit in finalists:
                if selector == "Strand":
                    feat_strand = Features.intervals[hit[FEAT]][0].strand
                    read = (read[0], feat_strand)
                if read not in FeatureSelector.rules_table[hit[RULE]][selector]:
                    eliminated.add(hit)
                    if self.report_eliminations:
                        feat_class = Features.attributes[hit[FEAT]][0][1][0]
                        self.elim_stats[feat_class][f"{selector}={read}"] += 1

            finalists -= eliminated
            eliminated.clear()

            if not finalists: return set()

        # Remaining finalists have passed all filters
        return {choice[FEAT] for choice in finalists}

    @staticmethod
    def is_perfect_iv_match(feat_start, feat_end, aln_iv) -> bool:
        # Returns true if the alignment interval is fully enclosed by the feature interval
        return feat_start <= aln_iv.start and feat_end >= aln_iv.end

    def choose_identities(self, feats_list: List[Tuple[int, int, Set[str]]], aln_iv: 'HTSeq.GenomicInterval') -> Set[Hit]:
        """Performs the initial selection on the basis of identity rules: attribute (key, value)

        Feature candidates are supplied to this function via feats_list. This is a list
        of tuples, each representing features associated with an interval which
        overlapped the alignment interval. The interval in this tuple may be a partial
        (incomplete) overlap with the alignment.

        The feats_list takes the following form:
            [(iv_A_start, iv_A_end, {features, associated, with, iv_A, ... }),
             (iv_B_start, iv_B_end, {features, associated, with, iv_B, ... }), ... ]
            Where iv_A and iv_B overlap aln_iv by at least 1 base

        Args:
            feats_list: a list of tuples, each representing features associated with
                an interval which overlapped the alignment interval. See above.
            aln_iv: the GenomicInterval of the alignment to which we are trying to
                assign features.

        Returns:
            identity_hits: a set of Hits for features which matched an identity rule,
            after performing elimination by hierarchy.
        """

        finalists, identity_hits = set(), list()

        for iv_feats in feats_list:
            # Check for perfect interval match only once per IntervalFeatures
            perfect_iv_match = self.is_perfect_iv_match(iv_feats[0], iv_feats[1], aln_iv)
            for feat in iv_feats[2]:
                for ident in Features.attributes[feat]:
                    for rule in self.get_identity_matches(ident):
                        if not perfect_iv_match and FeatureSelector.rules_table[rule]['Strict']:
                            continue
                        identity_hits.append((FeatureSelector.rules_table[rule]['Hierarchy'], rule, feat))
        # -> identity_hits: [(hierarchy, rule, feature), ...]

        # Perform hierarchy-based elimination
        if len(identity_hits) == 1:
            finalists.add(identity_hits[0])
        elif len(identity_hits) > 1:
            uniq_ranks = {hit[RANK] for hit in identity_hits}

            if len(identity_hits) == len(uniq_ranks):
                finalists.add(min(identity_hits, key=lambda x: x[RANK]))
            else:
                # Two or more hits share the same hierarchy.
                min_rank = min(uniq_ranks)
                finalists.update(hit for hit in identity_hits if hit[RANK] == min_rank)

        return finalists

    @classmethod
    def get_identity_matches(cls, ident) -> Iterator[int]:
        """Returns indexes of rules associated with a feature attribute record

        An attribute record is of the form ('key', ('value1', 'value2', ...)) where
        one or more values may be associated with the key. If there are multiple values,
        key-value products are formed, i.e. ('key', 'value1'), ('key', 'value2'), ...,
        for lookup in the inverted identities table. In this way multiple identities may
        be produced from a single attribute record. If any identity product matches a
        rule in the inverted identities table, the indexes of the associated selection
        rules are returned
        """

        try:
            # Check if rules are defined for this feature identity
            for identity_rule in cls.inv_ident[ident]:
                yield identity_rule
        except KeyError:
            pass

    @classmethod
    def get_all_identity_matches(cls) -> Set[str]:
        """Returns all features which match an identity rule in the rules table"""

        matches = set()
        for feat, feat_attrs in Features.attributes.items():
            for attr in feat_attrs:
                ident_matches = cls.get_identity_matches(attr)
                try:
                    next(ident_matches)
                    matches.add(feat)
                except StopIteration:
                    # No rules defined for this identity
                    pass

        return matches

    @staticmethod
    def build_selectors(rules_table) -> List[dict]:
        """Builds single/list/range/wildcard membership-matching selectors.

        Applies to: strand, 5' end nucleotide, and length

        This function replaces text-based selector definitions in the rules_table with
        their corresponding selector classes. Selector evaluation is then performed via
        the membership operator (keyword `in`) which is handled by the selector class'
        __contains__() method.
        """

        rules_table = sorted(rules_table, key=lambda x: x['Hierarchy'])

        selector_builders = {"Strand": StrandMatch, "nt5end": NtMatch, "Length": NumericalMatch}
        for row in rules_table:
            for selector, build_fn in selector_builders.items():
                defn = row[selector]
                if type(defn) is str and any([wc in defn.lower() for wc in ['all', 'both']]):
                    row[selector] = Wildcard()
                else:
                    row[selector] = build_fn(defn)

        return rules_table

    @staticmethod
    def build_inverted_identities(rules_table) -> Dict[Tuple[str, str], List[int]]:
        """Builds inverted identity rules for fast matching in phase 1 selection

        The resulting dictionary has (Attrib Key, Attrib Val) as key, [associated rule indexes] as val
        """

        inverted_identities = defaultdict(list)
        for i, rule in enumerate(rules_table):
            inverted_identities[rule['Identity']].append(i)

        return dict(inverted_identities)

    @staticmethod
    def get_hit_indexes() -> Tuple[int, int, int]:
        """Hits are stored as tuples for performance. This returns a human friendly index map for the tuple."""
        return RANK, RULE, FEAT