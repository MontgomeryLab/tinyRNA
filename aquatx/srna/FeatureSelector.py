from collections import defaultdict, Counter
from typing import List, Tuple, FrozenSet
from dataclasses import dataclass
from HTSeq import SAM_Alignment
import re


class FeatureSelector:

    # filter types: membership, range, list, wildcard
    # Identifier/Class: provided by attributes table
    # Strand, 5pnt, Length: provided by assign_features()

    range_match = re.compile(r"(\d+)-(\d+)")
    rank, rule, feat, type = 0, 1, 2, 3
    filtered = (None, None, ('Filtered',), None)
    unknown  = (None, None, ('Unknown',), None)

    def __init__(self, rules: List[dict]):
        # Todo: now that filters are type-specific rather than column specific,
        #  we can specify filter type in the input rules list such that feature selection
        #  is no longer bound to this specific interest set.
        #  This will improve the project's maintainability.
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

    def choose_identities(self, feats_set):
        """Performs the initial selection using identity rules (attribute key, value)"""

        choices, finalists = set(), set()

        # For each feature, match across all identity interests
        identity_hits = [(self.rules_table[rule]['Hierarchy'], rule, feat, type)
                         for feat, type in feats_set
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
        if len(ih_unique) < len(feats_set): choices.add(self.unknown)
        return choices, finalists

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

        attributes[('Filtered',)] = {}
        attributes[('Unknown',)] = {}
        self.attributes = attributes


class StatsCollector:
    @dataclass
    class Bundle:
        feat_counts = Counter()
        type_counts = Counter()

        bundle_len: int
        dup_counts: int
        cor_counts: float

    def __init__(self, selector):
        self.select = selector
        self.nt_len_mat = {nt: Counter() for nt in ['A', 'T', 'G', 'C']}
        self.counts = {feat: 0 for feat in sorted(selector.attributes.keys(), key=lambda x: x[0])}
        self.bundle_counts = {}
        self.class_counts = Counter()
        self.stats_counts = {stat: 0 for stat in
                        ['_aligned_reads', '_aligned_reads_unique_mapping', '_aligned_reads_multi_mapping',
                         '_unique_sequences_aligned', '_reads_unique_features', '_alignments_unique_features',
                         '_ambiguous_alignments_classes', '_ambiguous_reads_classes', '_ambiguous_alignments_features',
                         '_ambiguous_reads_features', '_no_feature']}

    def count_bundle(self, aln_bundle: iter) -> int:
        bundle_size = len(aln_bundle)
        bundle_hash = hash(aln_bundle)
        sequence = aln_bundle[0].read

        # Calculate counts for multi-mapping
        dup_counts = int(sequence.name.split('=')[1])
        cor_counts = dup_counts / bundle_size

        # fill in 5p nt/length matrix
        self.nt_len_mat[chr(sequence[0])][len(sequence)] += dup_counts

        self.stats_counts['_aligned_reads'] += dup_counts
        self.stats_counts['_unique_sequences_aligned'] += 1
        self.stats_counts['_aligned_reads_multi_mapping'] += dup_counts * (bundle_size > 1)
        self.stats_counts['_aligned_reads_unique_mapping'] += dup_counts * (bundle_size == 1)

        # Save bundle stats to later update with per-alignment stats
        self.bundle_counts[bundle_hash] = self.Bundle(bundle_size, dup_counts, cor_counts)

        return bundle_hash

    # Needs work
    def count_alignment(self, bundle_hash, selected_feats, n_candidates):
        bundle = self.bundle_counts[bundle_hash]
        selection_size = len(selected_feats)
        unique_feats = {f[self.select.feat] for f in selected_feats}
        unique_types = {t[self.select.type] for t in selected_feats}

        if len(unique_types) > 1:
            bundle.type_counts["ambiguous"] += bundle.cor_counts
        elif len(unique_feats) > 1:
            for feat in selected_feats:
                feat = feat[self.select.feat]
                bundle.feat_counts[feat] += bundle.cor_counts / selection_size
        else:
            single = selected_feats.pop()
            if selection_size == 1 and single[self.select.feat] == 'Unknown':
                self.stats_counts['_no_feature'] += bundle.cor_counts
            else:
                bundle.type_counts[single[self.select.type]] += bundle.cor_counts
                bundle.feat_counts[single[self.select.feat]] += bundle.cor_counts

    def finalize_bundle(self, bundle_hash):
        bundle = self.bundle_counts[bundle_hash]

        if len(bundle.type_counts) > 1:
            self.class_counts["ambiguous"] += sum(bundle.type_counts.values())
            self.stats_counts['_ambiguous_alignments_classes'] += 1
            self.stats_counts['_ambiguous_reads_classes'] += bundle.dup_counts
        elif len(bundle.feat_counts) > 1:
            # Needs further conversion
            key = next(iter(bundle.type_counts))
            self.stats_counts['_ambiguous_alignments_features'] += 1
            self.stats_counts['_ambiguous_reads_features'] += bundle.dup_counts
            self.class_counts[key] += bundle.type_counts[key]
            for key, value in bundle.feat_counts.items():
                bundle.feat_counts += value

    def write_reports(self, prefix):
        pass

    def merge(self, other):
        pass

    """
	stats_counts:
	• Before feature assignment (upon acquisition of each bundle)
		○ _unique_sequences_aligned: increment once per unique read (bundle)
		○ _aligned_reads: increment by dup counts (fasta header) for each bundle, regardless
		○ _aligned_reads_multi_mapping: increment by dup counts if bundle length is > 1
		○ _aligned_reads_unique_mapping: increment by dup counts if bundle is singular
		○ _no_feature: increment by cor counts if no features were assigned
    • After feature assignment (per-bundle)
		○ See bundle_feats and bundle_class below
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
		• class: at aln_classes.item(), increment by cor counts if assigned classes length > 1

	class_counts (after processing each bundle)
		• ambiguous: increment by sum of bundle_class values if bundle_class length is greater than 1

	STATS FILES:
	• If intermediate files (args.out_prefix + '_out_aln_table.txt'):
		○ For each alignment in the bundle, write tab delimited:
			Aln.read	Cor_counts	Aln.iv.strand	Aln.iv.start	Aln.iv.end	C1;C2;...	F1;F2;...
	• At conclusion of counting (args.out_prefix + '_stats.txt')
		○ Header: Summary Statistics
		○ Tab delimited stat -> count, one per line
		○ Tab delimited _no_feature -> count, one per line
	• Final CSV outputs from pandas DataFrames
		○ Class_counts (args.out_prefix + '_out_class_counts.csv')
		○ Feature_counts (args.out_prefix + '_out_feature_counts.txt')
		○ Nt_len_mat (args.out_prefix + '_out_nt_len_dist.csv')

    """
