from collections import defaultdict
from typing import List, Tuple, FrozenSet
import re


class Wildcard:
    @staticmethod
    def __contains__(val): return True
    def __repr__(self): return "all"


class FeatureSelector:

    # filter types: membership, range, list, wildcard
    # Identifier/Class: provided by attributes table
    # Strand, 5pnt, Length: provided by assign_features()
    range_match = re.compile(r"(\d+)-(\d+)")
    rank,rule,feat,type = 0,1,2,3
    filtered = (None, None, ('Filtered',), None)
    unknown  = (None, None, ('Unknown',), None)

    def __init__(self, rules: List[dict]):
        self.interest = ('Identity', 'Strand', '5pnt', 'Length')
        self.rules_table = sorted(rules, key=lambda x: int(x['Hierarchy']))
        self.build_filters()
        self.attributes = []

        # Inverted ident rules: (Attrib Key, Attrib Val) as key, [rule matches] as val
        # This allows us to do O(1) identity matching for each candidate feature
        # Subsequent matching is then performed only against identity-matched rows
        inverted_identities = defaultdict(list)
        for i, rule in enumerate(self.rules_table):
            inverted_identities[rule['Identity']].append(i)
        self.inv_ident = dict(inverted_identities)

    def choose(self, feat_set, strand, endnt, length) -> set:
        # Perform an efficient first-round
        choices, finalists = self.choose_identities(feat_set)
        if not finalists: return choices

        # Strand, 5pnt, and Length filtering uses simpler logic
        eliminated = set()
        for step, read in zip(self.interest[1:], (strand, endnt, length)):
            for hit in finalists:
                if read not in self.rules_table[hit[self.rule]][step]:
                    eliminated.add(hit)
                    choices.add(self.filtered)
            finalists -= eliminated
            eliminated.clear()
            if not finalists:
                return choices

        choices.update(finalists)
        return choices

    def choose_identities(self, feat_set):
        """Performs the initial selection using identity rules (Identifier, Feature)"""

        choices, finalists = set(), set()

        # For each feature, match across all identity interests
        identity_hits = [(self.rules_table[rule]['Hierarchy'], rule, feat, type)
                         for feat, type in feat_set
                         for akey, av in self.attributes[feat]
                         for attr_val in av
                         for rule in self.inv_ident.get((akey, attr_val), ())
                         if (akey, attr_val) in self.inv_ident]
        # -> [(hierarchy, rule, feature, type), ...]

        if len(identity_hits) == 1:
            # Only one feature matched only one rule
            finalists.add(identity_hits[0])
        elif len(identity_hits) > 1:
            # Perform any possible hierarchy-based eliminations
            rank_set = {hit[self.rank] for hit in identity_hits}

            if len(identity_hits) == len(rank_set):
                finalists.add(min(identity_hits, key=lambda x: x[self.rank]))
            else:
                # Two or more hits share the same hierarchy.
                min_rank = min(rank_set)
                finalists.update(hit for hit in identity_hits if hit[self.rank] == min_rank)

        ih_unique = {hit[self.feat] for hit in identity_hits}
        if len(finalists) < len(ih_unique): choices.add(self.filtered)
        if len(ih_unique) < len(feat_set): choices.add(self.unknown)
        return choices, finalists

    def build_filters(self):
        """Builds single/list/range/wildcard membership-matching filters

        Any combination of the above filter types may be present in each rule.
        These filter types are only supported for 5' End Nucleotide and Length.
        """

        def wildcard(step) -> bool:
            if "all" in row[step].lower():
                row[step] = Wildcard()
                return True

        def nt_filter() -> Tuple:
            rule = row["5pnt"].split(',')
            return tuple(map(lambda x: x.strip().upper(), rule))

        def numerical_filter() -> FrozenSet[int]:
            # Supports intermixed lists and ranges
            rule, lengths = row["Length"].split(','), []
            for piece in rule:
                if '-' in piece:
                    lo, hi = self.range_match.findall(piece)[0]
                    lengths.extend([*range(int(lo), int(hi) + 1)])
                else:
                    lengths.append(int(piece))

            return frozenset(lengths)

        filters = [("5pnt", nt_filter), ("Length", numerical_filter)]
        for row in self.rules_table:
            for step, filt in filters:
                if not wildcard(step):
                    row[step] = filt()

    def set_attributes_table(self, attributes, attrs_of_interest):
        """Attributes are not available at construction time; add later"""
        self.attributes = attributes
        self.attr_order = tuple(enumerate(attrs_of_interest))
        self.attributes[('Filtered',)] = {}
        self.attributes[('Unknown',)] = {}

class StatsCollector:
    def __init__(self):
        pass

    """
    stats_counts:
	• Before feature assignment (upon acquisition of each bundle)
		○ _unique_sequences_aligned: increment once per unique read (bundle)
		○ _aligned_reads: increment by dup counts (fasta header) for each bundle, regardless
		○ _aligned_reads_multi_mapping: increment by dup counts if bundle length is > 1
		○ _aligned_reads_unique_mapping: increment by dup counts if bundle is singular
		○ _no_feature: increment by cor counts if no features were assigned
	• After feature assignment (after processing each bundle completely)
		○ _reads_unique_features: increment by 1 if bundle_class and bundle_feats lengths are each == 1
		○ _alignments_unique_features: increment by 1 if bundle_class and bundle_feats lengths are each == 1
		○ _ambiguous_alignments_classes: increment by 1 if bundle_class length is > 1
		○ _ambiguous_reads_classes: increment by dup counts if bundle_class length is > 1
		○ _ambiguous_alignments_features: increment by 1 if bundle_feats length is > 1
        ○_ambiguous_reads_features: increment by dup counts if bundle_feats length is > 1
        
    bundle_feats (per-bundle)
	    • feature: at aln_feats.item(), increment by cor counts (this is for the case that assigned features is 1)

    bundle_class (per-bundle)
        • ambiguous: for all (unique) classes, if length > 1 then increment by 1
        • class: at aln_classes.item(), increment by corcounts if assigned classes length > 1
    
    class_counts (after processing each bundle)
        • ambiguous: increment by sum of bundle_class values if bundle_class length is greater than 1
        
    STATS FILES:
	• If intermediate files (args.out_prefix + '_out_aln_table.txt'):
		○ For each alignment in the bundle, write tab delimited:
		    Aln.read	Cor_counts  Aln.iv.strand	Aln.iv.start	Aln.iv.end  C1;C2;… F1;F2
	• At conclusion of counting (args.out_prefix + '_stats.txt')
		○ Header: Summary Statistics
		○ Tab delimited stat -> count, one per line
		○ Tab delimited _no_feature -> count, one per line
	• Final CSV outputs from pandas DataFrames
		○ Class_counts (args.out_prefix + '_out_class_counts.csv')
		○ Feature_counts (args.out_prefix + '_out_feature_counts.txt')
        ○ Nt_len_mat (args.out_prefix + '_out_nt_len_dist.csv')
        
    """
