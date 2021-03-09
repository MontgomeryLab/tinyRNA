import itertools
import re
from collections import defaultdict
from typing import List, Iterator


class FiltExec:
    """Executable filter unit to allow for wildcards and pipelining"""
    def __init__(self, rule: Iterator, step, index):
        self.rule = tuple(rule) if rule is not None else None
        self.exec = lambda x: x in self.rule if rule is not None else lambda x: True
        self.step = step
        self.idx  = index

    def __call__(self, query, *args, **kwargs):
        return self.exec(query)

    def __iter__(self):
        return iter(self.rule) if self.rule is not None else iter({True})

    def __hash__(self):
        return hash(self.rule) + hash(self.idx)

    def __repr__(self):
        return "%s: rule [%d]" % (str(self.rule), self.idx)


class FeatureSelector:
    # Keys determine which attributes are extracted when parsing GFF files (case must match GFF)
    # Values indicate corresponding index within each feature in the resulting attributes[] table
    ident_idx = [('Class', 0), ('biotype', 1)]

    # filter types: membership, range, list, wildcard
    # Identifier/Class: need to lookup in attributes table
    # Strand, 5pnt, Length: provided by assign_features()

    def __init__(self, rules: List[dict]):
        self.interest = ('Identity', 'Strand', '5pnt', 'Length')
        self.attributes, self.int_sets = [], []
        self.rules_table = sorted(rules, key=lambda x: int(x['Hierarchy']))
        self.compile_filters()
        # For each interest, create dict of preference: [associated rule indexes]
        for step in self.interest:
            inverted_rules = defaultdict(list)
            for row, rule in enumerate(self.rules_table):
                inverted_rules[rule[step]].append(row)
            self.int_sets.append(dict(inverted_rules))

    def choose(self, feat_set, strand, endnt, length) -> set:
        step = iter(self.int_sets)

        # Identity
        identity_interests = next(step)
        choices, finalists = self.identity_choice(feat_set, identity_interests)
        if not finalists: return choices

        # Strand, 5pnt, and Length filtering uses simpler logic
        for step, read in zip(self.interest[1:], [strand, endnt, length]):
            eliminations = set()
            for feat in finalists:
                rule_at_feat = feat[2]
                if read not in self.rules_table[rule_at_feat][step]:
                    eliminations.add(feat)
                    choices.add('Filtered')
            finalists -= eliminations
            if not finalists: return choices

    def compile_filters(self):
        range_match = re.compile(r"(\d+)-(\d+)")
        for i,row in enumerate(self.rules_table):
            for step in ["5pnt", "Length"]:
                rule = row[step]
                if "all" in rule.lower():
                    # Wildcard filter
                    row[step] = lambda x: True
                    continue

                rule = rule.split(',')
                if step is "5pnt":
                    rule = map(lambda x: x.strip().upper(), rule)
                else:
                    # Length supports intermixed lists and ranges
                    lengths = []
                    for piece in rule:
                        if '-' not in piece:
                            lengths.append(int(piece))
                            continue
                        lo, hi = range_match.findall(piece)[0]
                        lengths.extend([*range(int(lo), int(hi)+1)])
                    rule = map(int, lengths)

                row[step] = FiltExec(rule, step, i)

    def identity_choice(self, feat_set, interests):
        # For each feature, match across all identity interests
        identity_hits = [(feat, interests[(id, attr)])
                         for feat in feat_set
                         for id, i in self.ident_idx
                         for attr in self.attributes[feat][i]
                         if (id, attr) in interests]
        # -> [(feature, [matched rule indexes, ...]), ...]

        choices, finalists = set(), set()

        if len(identity_hits) < len(feat_set):
            choices.add('Unknown')
        if len(identity_hits) == 1 and len(identity_hits[0][1]) == 1:
            # Only one feature matched only one rule
            finalists.add(identity_hits[0][0])
        else:
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

        if len(finalists) < len(identity_hits): choices.add('Unknown')
        return choices, finalists

    def set_attributes_table(self, attributes):
        self.attributes = attributes
