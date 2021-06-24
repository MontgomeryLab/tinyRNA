import itertools
import HTSeq
import re

from collections import defaultdict
from typing import List, Tuple, FrozenSet, Dict, Set

# Type aliases for human readability
IntervalFeatures = Tuple[int, int, Set[str]]  # A set of features associated with an interval


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

    rank, rule, feat = 0, 1, 2
    attributes = {}

    def __init__(self, rules: List[dict], reference_table: Dict):
        FeatureSelector.attributes = reference_table
        self.interest = ('Identity', 'Strand', 'nt5', 'Length')
        self.rules_table = sorted(rules, key=lambda x: x['Hierarchy'])
        self.build_filters()

        # Inverted ident rules: (Attrib Key, Attrib Val) as key, [associated rules] as val
        inverted_identities = defaultdict(list)
        for i, rule in enumerate(self.rules_table):
            inverted_identities[rule['Identity']].append(i)
        self.inv_ident = dict(inverted_identities)

    def choose(self, feat_set, alignment) -> set:
        # Perform hierarchy-based first round of selection for identities
        finalists = self.choose_identities(feat_set, alignment.iv)
        if not finalists: return set()

        strand = alignment.iv.strand
        nt5end = alignment.read.nt5
        length = len(alignment.read)

        eliminated = set()
        for step, read in zip(self.interest[1:], (strand, nt5end, length)):
            for hit in finalists:
                if read not in self.rules_table[hit[self.rule]][step]:
                    eliminated.add(hit)

            finalists -= eliminated
            eliminated.clear()

            if not finalists: return set()

        # Remaining finalists have passed all filters
        return {choice[self.feat] for choice in finalists}

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

        # Only one feature matched only one rule
        if len(identity_hits) == 1:
            finalists.add(identity_hits[0])
        # Perform any possible hierarchy-based eliminations
        elif len(identity_hits) > 1:
            uniq_ranks = {hit[self.rank] for hit in identity_hits}

            if len(identity_hits) == len(uniq_ranks):
                finalists.add(min(identity_hits, key=lambda x: x[self.rank]))
            else:
                # Two or more hits share the same hierarchy.
                min_rank = min(uniq_ranks)
                finalists.update(hit for hit in identity_hits if hit[self.rank] == min_rank)

        return finalists

    def build_filters(self):
        """Builds single/list/range/wildcard membership-matching filters

        Any combination of the above filter types may be present in each rule.
        These filter types are only supported for 5' End Nucleotide and Length.
        """

        class Wildcard:
            @staticmethod
            def __contains__(x): return True
            def __repr__(self): return "all"

        def wildcard(step) -> bool:
            if "all" in row[step].lower():
                row[step] = Wildcard()
                return True

        def nt_filter() -> Tuple:
            rule = row["nt5"].split(',')
            return tuple(map(lambda x: x.strip().upper(), rule))

        def numerical_filter() -> FrozenSet[int]:
            # Supports intermixed lists and ranges
            rule, lengths = row["Length"].split(','), []
            for piece in rule:
                if '-' in piece:
                    lo, hi = re.findall(r"(\d+)-(\d+)", piece)[0]
                    lengths.extend([*range(int(lo), int(hi) + 1)])
                else:
                    lengths.append(int(piece))

            return frozenset(lengths)

        filters = [("nt5", nt_filter), ("Length", numerical_filter)]
        for row in self.rules_table:
            for step, filt in filters:
                if not wildcard(step):
                    row[step] = filt()

    @classmethod
    def get_hit_indexes(cls):
        """Hits are stored as tuples for performance. This returns a human friendly index map for the tuple."""
        return cls.rank, cls.rule, cls.feat