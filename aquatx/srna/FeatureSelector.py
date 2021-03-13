from collections import defaultdict
from typing import List
import re


class Wildcard:
    @staticmethod
    def __contains__(val): return True
    def __repr__(self): return "all"


class FeatureSelector:

    # filter types: membership, range, list, wildcard
    # Identifier/Class: need to lookup in attributes table
    # Strand, 5pnt, Length: provided by assign_features()

    def __init__(self, rules: List[dict]):
        self.interest = ('Identity', 'Strand', '5pnt', 'Length')
        self.rules_table = sorted(rules, key=lambda x: int(x['Hierarchy']))
        self.build_filters()
        self.attributes = []

        # Inverted ident rules: (Identifier,Feature) as key, [rule index(es)] as val
        # This allows us to do O(1) identity matching for each candidate feature
        # Subsequent matching is then performed only against identity-matched rows
        self.inv_ident = defaultdict(list)
        for i, rule in enumerate(self.rules_table):
            self.inv_ident[rule['Identity']].append(i)

    def choose(self, feat_set, strand, endnt, length) -> set:
        # Identity
        choices, finalists = self.choose_identities(feat_set)
        if not finalists: return choices

        # Strand, 5pnt, and Length filtering uses simpler logic
        for step, read in zip(self.interest[1:], (strand, endnt, length)):
            eliminations = set()
            for feat in finalists:
                matched_rule_idx = feat[2]
                if read not in self.rules_table[matched_rule_idx][step]:
                    eliminations.add(feat)
                    choices.add('Filtered')
            finalists -= eliminations
            if not finalists:
                return choices

        choices.update({f[1] for f in finalists})
        return choices

    def choose_identities(self, feat_set):
        """Performs the initial selection using identity rules (Identifier, Feature)"""

        # For each feature, match across all identity interests
        identity_hits = [(feat, self.inv_ident[(key, attr)])
                         for feat, type in feat_set
                         for i, key in self.attr_order
                         for attr in self.attributes[feat][i]
                         if (key, attr) in self.inv_ident]
        # -> [(feature, [matched rule indexes, ...]), ...]

        choices, finalists = set(), set()

        if len(identity_hits) < len(feat_set):
            choices.add('Unknown')
        if len(identity_hits) == 1 and len(identity_hits[0][1]) == 1:
            # Only one feature matched only one rule
            feat = identity_hits[0][0]
            rule = identity_hits[0][1][0]
            rank = self.rules_table[rule]['Hierarchy']
            finalists.add((rank, feat, rule))
        elif len(identity_hits) > 1:
            # Perform any possible hierarchy-based eliminations
            ranks = [self.rules_table[rule]['Hierarchy'] for hit in identity_hits for rule in hit[1]]
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
                finalists.update(hit for hit in hit_rank if get_rank(hit) == min_rank)

        if len(finalists) < len(identity_hits): choices.add('Filtered')
        return choices, finalists

    def build_filters(self):
        """Builds single/list/range/wildcard membership-matching filters

        Any combination of the above filter types may be present in each rule.
        These filter types are only supported for 5' End Nucleotide and Length.
        """

        range_match = re.compile(r"(\d+)-(\d+)")
        for row in self.rules_table:
            for step in ["5pnt", "Length"]:
                rule = row[step]
                if "all" in rule.lower():
                    # Infinite membership ~
                    row[step] = Wildcard()
                    continue

                rule = rule.split(',')
                if step == "5pnt":
                    rule = map(lambda x: x.strip().upper(), rule)
                else:
                    # Length supports intermixed lists and ranges
                    lengths = []
                    for piece in rule:
                        if '-' not in piece:
                            lengths.append(int(piece))
                            continue
                        lo, hi = range_match.findall(piece)[0]
                        lengths.extend([*range(int(lo), int(hi) + 1)])
                    rule = map(int, lengths)

                row[step] = tuple(rule)

    def set_attributes_table(self, attributes, attrs_of_interest):
        """Attributes are not available at construction time; add later"""
        self.attributes = attributes
        self.attr_order = tuple(enumerate(attrs_of_interest))
        self.attributes['Filtered'] = {}
        self.attributes['Unknown'] = {}
