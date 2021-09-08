import itertools
import HTSeq

from collections import defaultdict
from typing import List, Tuple, Set

import tiny.rna.counter.hts_parsing as parser
from .selectors import Wildcard, StrandMatch, NumericalMatch, NtMatch
from .statistics import LibraryStats

# Type aliases for human readability
IntervalFeatures = Tuple[int, int, Set[str]]  # A set of features associated with an interval
AssignedFeatures = set
N_Candidates = int

BOTH_STRANDS = ('+', '-')

# Global indexes for Hits produced by choose_identity()
RANK, RULE, FEAT = 0, 1, 2


class Features:
    chrom_vectors: HTSeq.ChromVector
    attributes: dict
    intervals: dict
    _instance = None  # Singleton

    def __init__(self, features: HTSeq.GenomicArrayOfSets, attributes: dict, intervals: dict):
        if Features._instance is None:
            Features.chrom_vectors = features.chrom_vectors  # For interval -> feature ID lookups
            Features.attributes = attributes                 # For feature ID -> GFF column 9 attribute lookups
            Features.intervals = intervals                   # For feature ID -> interval lookups
            Features._instance = self


class FeatureCounter:
    out_prefix: str
    run_diags: bool
    alias: dict

    def __init__(self, gff_file_set, selection_rules, run_diags, out_prefix):
        reference_tables = parser.build_reference_tables(gff_file_set, selection_rules)
        Features(reference_tables[0], reference_tables[1], reference_tables[3])

        FeatureCounter.alias = reference_tables[2]
        FeatureCounter.out_prefix = out_prefix
        FeatureCounter.run_diags = run_diags

        self.stats = LibraryStats(out_prefix, report_diags=run_diags)
        self.selector = FeatureSelector(selection_rules, self.stats, diags=run_diags)

    def assign_features(self, alignment: 'parser.Alignment') -> Tuple[AssignedFeatures, N_Candidates]:
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

    If more than one hit remains following first round selection, a second round of selection
    is performed against sequence attributes: strand, 5' end nucleotide, and length. Rules for
    5' end nucleotides support lists (e.g. C,G,U) and wildcards (e.g. "all"). Rules for length
    support lists, wildcards, and ranges (i.e. 20-27) which may be intermixed in the same rule.
    Lengths may be specified as "strict", meaning that the feature must be completely contained
    by the alignment interval.
    """

    def __init__(self, rules: List[dict], libstats: 'LibraryStats', diags=False):
        self.rules_table = self.build_selectors(rules)
        self.inv_ident = self.build_inverted_identities(self.rules_table)

        self.phase1_candidates = libstats.identity_roster
        self.report_eliminations = diags
        if diags: self.elim_stats = libstats.diags.selection_diags

    def choose(self, feat_list: List[IntervalFeatures], alignment: 'parser.Alignment') -> set:
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
                    feat_strand = Features.intervals[hit[FEAT]].strand
                    read = (read[0], feat_strand)
                if read not in self.rules_table[hit[RULE]][selector]:
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
    def is_perfect_iv_match(feat_start, feat_end, aln_iv):
        # Only accept perfect interval matches for rules requiring such
        return feat_start <= aln_iv.start and feat_end >= aln_iv.end

    def choose_identities(self, feats_list: List[IntervalFeatures], aln_iv: 'HTSeq.GenomicInterval'):
        """Performs the initial selection on the basis of identity rules: attribute (key, value)

        Feature candidates are supplied to this function via feats_list. This is a list
        of tuples, each representing features associated with an in interval which
        overlapped the alignment interval. The interval in this tuple may be a partial
        (incomplete) overlap with the alignment.

        The list of IntervalFeatures takes the following form:
            [(iv_A_start, iv_A_end, {features, associated, with, iv_A, ... }),
             (iv_B_start, iv_B_end, {features, associated, with, iv_B, ... }), ... ]
            Where iv_A and iv_B overlap aln_iv by at least 1 base

        Args:
            feats_list: a list of tuples, each representing features associated with
                an interval which overlapped the alignment interval. See above.
            aln_iv: the GenomicInterval of the alignment to which we are trying to
                assign features.

        Returns:

        """

        finalists, identity_hits = set(), list()
        start, end, features = 0, 1, 2  # IntervalFeatures tuple indexes

        for iv_feats in feats_list:
            # Check for perfect interval match only once per IntervalFeatures
            perfect_iv_match = self.is_perfect_iv_match(iv_feats[start], iv_feats[end], aln_iv)
            for feat in iv_feats[features]:
                for attrib in Features.attributes[feat]:
                    # If multiple values are associated with the attribute key, create their key/value products
                    for feat_ident in itertools.product([attrib[0]], attrib[1]):
                        try:
                            # Check if rules are defined for this feature identity
                            for rule in self.inv_ident[feat_ident]:
                                if not perfect_iv_match and self.rules_table[rule]['Strict']:
                                    continue
                                identity_hits.append((self.rules_table[rule]['Hierarchy'], rule, feat))
                        except KeyError:
                            pass
        # -> identity_hits: [(hierarchy, rule, feature), ...]

        self.phase1_candidates.update(hit[FEAT] for hit in identity_hits)

        # Only one feature matched only one rule
        if len(identity_hits) == 1:
            finalists.add(identity_hits[0])
        # Perform any possible hierarchy-based eliminations
        elif len(identity_hits) > 1:
            uniq_ranks = {hit[RANK] for hit in identity_hits}

            if len(identity_hits) == len(uniq_ranks):
                finalists.add(min(identity_hits, key=lambda x: x[RANK]))
            else:
                # Two or more hits share the same hierarchy.
                min_rank = min(uniq_ranks)
                finalists.update(hit for hit in identity_hits if hit[RANK] == min_rank)

        return finalists

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
    def build_inverted_identities(rules_table) -> dict:
        """Builds inverted identity rules for fast matching in phase 1 selection

        The resulting dictionary has (Attrib Key, Attrib Val) as key, [associated rule indexes] as val
        """

        inverted_identities = defaultdict(list)
        for i, rule in enumerate(rules_table):
            inverted_identities[rule['Identity']].append(i)

        return dict(inverted_identities)

    @classmethod
    def get_hit_indexes(cls):
        """Hits are stored as tuples for performance. This returns a human friendly index map for the tuple."""
        return RANK, RULE, FEAT
