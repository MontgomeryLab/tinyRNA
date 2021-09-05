import itertools
import HTSeq
import re

from collections import defaultdict, Counter
from typing import List, Tuple, FrozenSet, Dict, Set

from .hts_parsing import Alignment
from .statistics import LibraryStats

# Type aliases for human readability
IntervalFeatures = Tuple[int, int, Set[str]]  # A set of features associated with an interval

# Global indexes for Hits produced by choose_identity()
RANK, RULE, FEAT = 0, 1, 2


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

    attributes = {}
    intervals = {}

    def __init__(self, rules: List[dict], reference_table: Dict, intervals: Dict, libstats: 'LibraryStats', diags=False):
        FeatureSelector.attributes = reference_table
        FeatureSelector.intervals = intervals
        self.rules_table = sorted(rules, key=lambda x: x['Hierarchy'])
        self.build_selectors()

        # Inverted ident rules: (Attrib Key, Attrib Val) as key, [associated rules] as val
        inverted_identities = defaultdict(list)
        for i, rule in enumerate(self.rules_table):
            inverted_identities[rule['Identity']].append(i)
        self.inv_ident = dict(inverted_identities)

        self.set_stats_collectors(libstats, diags)

    def choose(self, feat_list: List[IntervalFeatures], alignment: 'Alignment') -> set:
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
                    feat_strand = FeatureSelector.intervals[hit[FEAT]].strand
                    read = (read[0], feat_strand)
                if read not in self.rules_table[hit[RULE]][selector]:
                    eliminated.add(hit)
                    if self.report_eliminations:
                        feat_class = FeatureSelector.attributes[hit[FEAT]][0][1][0]
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
                for attrib in FeatureSelector.attributes[feat]:
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

    def build_selectors(self):
        """Builds single/list/range/wildcard membership-matching selectors."""

        selector_builders = {"Strand": StrandMatch, "nt5end": NtMatch, "Length": NumericalMatch}
        for row in self.rules_table:
            for selector, build_fn in selector_builders.items():
                defn = row[selector]
                if type(defn) is str and any([wc in defn.lower() for wc in ['all', 'both']]):
                    row[selector] = Wildcard()
                else:
                    row[selector] = build_fn(defn)

    def set_stats_collectors(self, libstats: 'LibraryStats', diags=False):
        """Enables diagnostic reporting for eliminations made in selection phase 2"""

        self.phase1_candidates = libstats.identity_roster
        if diags:
            self.elim_stats = libstats.selection_diags
            self.report_eliminations = diags

    @classmethod
    def get_hit_indexes(cls):
        """Hits are stored as tuples for performance. This returns a human friendly index map for the tuple."""
        return RANK, RULE, FEAT


class Wildcard:
    @staticmethod
    def __contains__(x): return True
    def __repr__(self): return "all"


class StrandMatch:
    """Evaluates BOTH the alignment's strand and the feature's strand for a match

    If sense: alignment's strand == feature strand for a match
    If antisense: alignment's strand != feature strand for a match
    """
    def __init__(self, strand):
        self.strand = strand.lower()
        self.select = (self.strand == 'sense')

    def __contains__(self, x):
        return self.select == (x[0] == x[1])

    def __repr__(self): return str(self.strand)


class NtMatch(tuple):
    """For evaluating a single nucleotide against a list of desired bases"""

    def __new__(cls, nts):
        nts = nts.split(',')
        return super().__new__(cls, map(lambda x: x.strip().upper(), nts))


class NumericalMatch(frozenset):
    """For evaluating sequence length against a list and/or range of desired values"""

    def __new__(cls, lengths):
        # Supports intermixed lists and ranges
        rule, lengths = lengths.split(','), []
        for piece in rule:
            if '-' in piece:
                for lo, hi in re.findall(r"(\d+)-(\d+)", piece):
                    lengths.extend([*range(int(lo), int(hi) + 1)])
            else:
                lengths.append(int(piece))

        return super().__new__(cls, lengths)