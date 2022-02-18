import HTSeq
import sys

from collections import defaultdict
from typing import List, Tuple, Set, Dict

from tiny.rna.counter.hts_parsing import ReferenceTables, SAM_reader
from .statistics import LibraryStats
from .matching import *

# Type aliases for human readability
match_tuple = Tuple[int, int, bool]
feature_record_tuple = Tuple[str, int, int, str, Tuple[match_tuple]]


class Features:
    chrom_vectors: HTSeq.ChromVector
    classes: dict
    aliases: dict

    _instance = None  # Singleton

    def __init__(self, features: HTSeq.GenomicArrayOfSets, aliases: dict, classes: dict):
        if Features._instance is None:
            Features.chrom_vectors = features.chrom_vectors  # For interval -> feature ID lookups
            Features.aliases = aliases                       # For feature ID -> preferred feature name lookups
            Features.classes = classes                       # For feature ID -> class lookups
            Features._instance = self


class FeatureCounter:
    out_prefix: str
    run_diags: bool

    def __init__(self, gff_file_set, selection_rules, **prefs):
        self.stats = LibraryStats(**prefs)
        self.sam_reader = SAM_reader(**prefs)
        self.selector = FeatureSelector(selection_rules, self.stats, **prefs)

        reference_tables = ReferenceTables(gff_file_set, self.selector, **prefs)
        Features(*reference_tables.get())

        FeatureCounter.out_prefix = prefs['out_prefix']
        FeatureCounter.run_diags = prefs['report_diags']

    def assign_features(self, al: dict) -> Tuple[set, int]:
        """Determines features associated with the interval then performs rule-based feature selection"""

        feat_matches, assignment = set(), set()

        try:
            feat_matches = feat_matches.union(*(match for match in
                    (Features.chrom_vectors[al['chrom']]['.']       # GenomicArrayOfSets -> ChromVector
                             .array[al['start']:al['end']]          # ChromVector -> StepVector
                             .get_steps(values_only=True))          # StepVector -> {features}
                    # If an alignment does not map to a feature, an empty set is returned
                    if len(match) != 0))
        except KeyError as ke:
            self.stats.chrom_misses[ke.args[0]] += 1

        # If features are associated with the alignment interval, perform selection
        if len(feat_matches):
            assignment = self.selector.choose(feat_matches, al)

        return assignment, len(feat_matches)

    def count_reads(self, library: dict):
        """Collects statistics on features assigned to each alignment associated with each read"""

        read_seq = self.sam_reader.bundle_multi_alignments(library["File"])
        self.stats.assign_library(library)

        # For each sequence in the sam file...
        for bundle in read_seq:
            bstat = self.stats.count_bundle(bundle)

            # For each alignment of the given sequence...
            for alignment in bundle:
                hits, n_candidates = self.assign_features(alignment)
                self.stats.count_bundle_assignments(bstat, alignment, hits, n_candidates)

            self.stats.finalize_bundle(bstat)

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
        (hierarchy, rule, feature_id, strand)

    If more than one Hit remains following first round selection, a second round of selection
    is performed against sequence attributes: strand, 5' end nucleotide, and length. Rules for
    5' end nucleotides support lists (e.g. C,G,U) and wildcards (e.g. "all"). Rules for length
    support lists, wildcards, and ranges (i.e. 20-27) which may be intermixed in the same rule.
    Lengths may be specified as "strict", meaning that the feature must be completely contained
    by the alignment interval.
    """

    rules_table: List[dict]
    inv_ident: Dict[tuple, List[int]]

    def __init__(self, rules: List[dict], libstats: 'LibraryStats', report_diags=False, **kwargs):
        FeatureSelector.rules_table = self.build_selectors(rules)
        FeatureSelector.inv_ident = self.build_inverted_identities(FeatureSelector.rules_table)

        self.report_eliminations = report_diags
        if report_diags: self.elim_stats = libstats.diags.selection_diags

    @classmethod
    def choose(cls, candidates: Set[feature_record_tuple], alignment: dict) -> Set[str]:
        """Selects features according to the selection rules provided at construction

        Feature candidates are supplied to this function via feats_list. This is a list
        of tuples, each representing features associated with an interval which
        overlapped the alignment interval. The interval in this tuple may be a partial
        (incomplete) overlap with the alignment.

        The feats_list takes the following form:
            [(iv_A_start, iv_A_end, {features associated with iv_A, ... }),
             (iv_B_start, iv_B_end, {features associated with iv_A, ... }), ... ]

            Each feature is represented as:
                (featureID, start, stop, strand, (match-tuple, ... ))

            Each match-tuple represents a rule which matched the feature on identity.
                (rule, rank, strict)

        Args:
            feats_list: a list of tuples, each representing features associated with
                an interval which overlapped the alignment interval. See above.
            alignment: the alignment to which features are being selected for assignment.

        Returns:
            selections: a list of features which passed selection
        """

        identity_hits, min_rank = [], sys.maxsize

        for feat, strand, matches in candidates:
            for rule, rank, iv_match in matches:
                if rank > min_rank: continue
                if alignment not in iv_match: continue
                if rank < min_rank: min_rank = rank
                identity_hits.append((rank, rule, feat, strand))
        # -> identity_hits: [(hierarchy, rule, feature, strand), ...]

        if not identity_hits: return set()

        selections = set()
        for hit in identity_hits:
            if hit[0] != min_rank: continue
            _, rule, feat, strand = hit

            strand = (alignment['strand'], strand)
            nt5end = alignment['nt5']
            length = alignment['len']

            rule = FeatureSelector.rules_table[rule]
            if strand not in rule["Strand"]: continue
            if nt5end not in rule["nt5end"]: continue
            if length not in rule["Length"]: continue
            selections.add(feat)

        return selections

    @staticmethod
    def build_selectors(rules_table) -> List[dict]:
        """Builds single/list/range/wildcard membership-matching selectors.

        Applies to: strand, 5' end nucleotide, and length

        This function replaces text-based selector definitions in the rules_table with
        their corresponding selector classes. Selector evaluation is then performed via
        the membership operator (keyword `in`) which is handled by the selector class'
        __contains__() method.
        """

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
    def build_interval_selectors(iv, match_tuples):
        """Builds partial/full/exact/3' anchored/5' anchored interval selectors

        Unlike build_selectors() and build_inverted_identities(), this function
        is not called at construction time. Instead, it is called when finalizing
        match tuples in ReferenceTables. This is because interval selectors are
        created for each feature (requiring start/stop/strand to be known) for
        each of the feature's identity matches (each match tuple).

        Args:
            iv: The interval of the feature from which each selector is built
            match_tuples: A list of tuples representing the feature's identity
                matches. Each tuple index 2 defines and is replaced by the selector.
        """

        selector_factory = {
            'full': lambda: IntervalFullMatch(iv),
            'exact': lambda: IntervalExactMatch(iv),
            'partial': lambda: IntervalPartialMatch(iv),
            "5' anchored": lambda: Interval5pMatch(iv),
            "3' anchored": lambda: Interval3pMatch(iv),
        }

        for i in range(len(match_tuples)):
            try:
                match = match_tuples[i]
                selector = selector_factory[match[2]]()
                match_tuples[i] = (match[0], match[1], selector)
            except KeyError:
                raise ValueError(f"Unrecognized match type: '{match_tuples[i][2]}'")

        return match_tuples

    @staticmethod
    def build_inverted_identities(rules_table) -> Dict[Tuple[str, str], List[int]]:
        """Builds inverted identity rules for fast matching in phase 1 selection

        The resulting dictionary has (Attrib Key, Attrib Val) as key, [associated rule indexes] as val
        """

        inverted_identities = defaultdict(list)
        for i, rule in enumerate(rules_table):
            inverted_identities[rule['Identity']].append(i)

        return dict(inverted_identities)