from collections import defaultdict, Counter
from typing import List, Tuple, FrozenSet, Dict

import HTSeq
import pandas as pd
import re


class FeatureSelector:
    """Performs hierarchical selection given a set of candidate features for a locus

    Two sources of data serve as targets for selection: feature attributes (sourced from
    input GFF files), and sequence attributes (sourced from in input SAM files).

    The first round of selection is performed against each candidate feature's attributes.
    The target for this stage is attribute key-value pairs, referred to here as Identities.
    A candidate may match multiple identities. Each match is referred to as a Hit. If more
    than one Hit is produced, elimination is performed using each Hit rule's hierarchy value.

    If more than one hit remains following first round selection, a second round of selection
    is performed against sequence attributes: strand, 5' end nucleotide, and length. Rules for
    5' end nucleotides support lists (e.g. C,G,U) and wildcards (e.g. "all"). Rules for length
    support lists, wildcards, and ranges (i.e. 20-27) which may be intermixed in the same rule.
    """

    rank, rule, feat = 0, 1, 2

    def __init__(self, rules: List[dict], reference_table: Dict):
        # Todo: now that filters are type-specific rather than column specific,
        #  we can specify filter type in the input rules list such that feature selection
        #  is no longer bound to this specific interest set.
        #  This will improve the project's maintainability.
        self.interest = ('Identity', 'Strand', '5pnt', 'Length')
        self.rules_table = sorted(rules, key=lambda x: x['Hierarchy'])
        self.attributes = reference_table
        self.build_filters()

        # Inverted ident rules: (Attrib Key, Attrib Val) as key, [rule matches] as val
        # This allows us to do O(1) identity matching for each candidate feature
        # Subsequent matching (also O(1)) is performed only against identity-matched rows
        inverted_identities = defaultdict(list)
        for i, rule in enumerate(self.rules_table):
            inverted_identities[rule['Identity']].append(i)
        self.inv_ident = dict(inverted_identities)

    def choose(self, feat_set, strand, endnt, length) -> Tuple[set, set]:
        # Perform an efficient first-round
        finalists, uncounted = self.choose_identities(feat_set)
        if not finalists: return finalists, uncounted

        # Strand, 5pnt, and Length filtering uses simpler logic
        choices, eliminated = set(), set()
        for step, read in zip(self.interest[1:], (strand, endnt, length)):
            for hit in finalists:
                if read not in self.rules_table[hit[self.rule]][step]:
                    eliminated.add(hit)
                    uncounted.add("Filtered")

            finalists -= eliminated
            eliminated.clear()

            if not finalists: return choices, uncounted

        # Remaining finalists have passed all filters
        choices.update(finalists)
        return choices, uncounted

    def choose_identities(self, feats_set):
        """Performs the initial selection using identity rules (attribute key, value)"""

        finalists, uncounted = set(), set()

        # For each feature, match across all identity interests
        identity_hits = [(self.rules_table[rule]['Hierarchy'], rule, feat)
                         for feat in feats_set
                         for attr_key, av in self.attributes[feat]
                         for attr_val in av
                         for rule in self.inv_ident.get((attr_key, attr_val), ())
                         if (attr_key, attr_val) in self.inv_ident]
        # -> [(hierarchy, rule, feature), ...]

        # Only one feature matched only one rule
        if len(identity_hits) == 1:
            finalists.add(identity_hits[0])
        # Perform any possible hierarchy-based eliminations
        elif len(identity_hits) > 1:
            rank_set = {hit[self.rank] for hit in identity_hits}
            feat_set = {hit[self.feat] for hit in identity_hits}

            # Scenarios:
            # 1. len(identity_hits) == len(rank_set) == len(feat_set): unique features and ranks -- covered
            # 2. len(identity_hits) == len(rank_set) != len(feat_set): unique ranks, 1 or more shared feat -- covered
            # 3. len(identity_hits) == len(feat_set) != len(rank_set): unique feats, 1 or more shared rank -- covered
            # 4. len(identity_hits) != len(rank_set) and != len(feat_set): shared feats, shared ranks, different number
            # 5. len(identity_hits) != len(rank_set) and != len(feat_set): shared feats, shared ranks, same number
                # Scenarios 4 & 5: needs further investigation/testing

            if len(identity_hits) == len(rank_set) == len(feat_set):
                finalists.add(min(identity_hits, key=lambda x: x[self.rank]))
            elif len(identity_hits) == len(rank_set) != len(feat_set):
                min_feat = min(identity_hits, key=lambda x: x[self.rank])
                finalists.add(min_feat)
                finalists.update(hit for hit in identity_hits if hit[self.feat] == min_feat[self.feat])
            else:
                # Two or more hits share the same hierarchy.
                min_rank = min(rank_set)
                finalists.update(hit for hit in identity_hits if hit[self.rank] == min_rank)

        ih_unique = {hit[self.feat] for hit in identity_hits}
        if len(finalists) < len(ih_unique): uncounted.add("Filtered")  # Todo: this logic is incorrect. Is relevant?
        if len(ih_unique) < len(feats_set): uncounted.add("Unknown")
        return finalists, uncounted

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
                    lo, hi = re.findall(r"(\d+)-(\d+)", piece)[0]
                    lengths.extend([*range(int(lo), int(hi) + 1)])
                else:
                    lengths.append(int(piece))

            return frozenset(lengths)

        filters = [("5pnt", nt_filter), ("Length", numerical_filter)]
        for row in self.rules_table:
            for step, filt in filters:
                if not wildcard(step):
                    row[step] = filt()

    @classmethod
    def get_hit_indexes(cls):
        return cls.rank, cls.rule, cls.feat


class StatsCollector:
    class Bundle:
        def __init__(self, size, dup_counts, cor_counts):
            self.feat_count = 0
            self.bundle_len = size
            self.dup_counts = dup_counts
            self.cor_counts = cor_counts

    rank, rule, feat = FeatureSelector.get_hit_indexes()

    def __init__(self, out_prefix: str = None, save_intermediate_file: bool = False):
        self.feat_counts = Counter()
        self.nt_len_mat = {nt: Counter() for nt in ['A', 'T', 'G', 'C']}
        self.stats_counts = {stat: 0 for stat in
                        ['_aligned_reads', '_aligned_reads_unique_mapping', '_aligned_reads_multi_mapping',
                         '_unique_sequences_aligned', '_reads_unique_features', '_alignments_unique_features',
                         '_ambiguous_alignments_features', '_ambiguous_reads_features', '_no_feature']}
        self.save_intermediate_file = save_intermediate_file
        self.out_prefix = out_prefix
        self.alignments = []

    def count_bundle(self, aln_bundle: iter) -> 'Bundle':
        bundle_read = aln_bundle[0].read
        bundle_size = len(aln_bundle)

        # Calculate counts for multi-mapping
        dup_counts = int(bundle_read.name.split('=')[1])
        cor_counts = dup_counts / bundle_size
        bundle_seq = bundle_read.seq

        # fill in 5p nt/length matrix
        self.nt_len_mat[chr(bundle_seq[0])][len(bundle_seq)] += dup_counts

        self.stats_counts['_aligned_reads'] += dup_counts
        self.stats_counts['_unique_sequences_aligned'] += 1
        self.stats_counts['_aligned_reads_multi_mapping'] += dup_counts * (bundle_size > 1)
        self.stats_counts['_aligned_reads_unique_mapping'] += dup_counts * (bundle_size == 1)  # Made ya look

        return self.Bundle(bundle_size, dup_counts, cor_counts)

    def count_bundle_alignments(self, bundle: 'Bundle', aln: HTSeq.SAM_Alignment, hits: set) -> None:
        for hit in hits:
            feature_corrected_count = bundle.cor_counts / len(hits)
            if hit[self.feat] not in self.feat_counts: bundle.feat_count += 1
            self.feat_counts[hit[self.feat]] += feature_corrected_count

        if self.save_intermediate_file:
            # sequence, cor_counts, strand, start, end, feat1;feat2a/feat2b/feat2c;feat3...
            self.alignments.append((aln.read, bundle.cor_counts, aln.iv.strand, aln.iv.start, aln.iv.end,
                                    ';'.join(map('/'.join, map(lambda x: x[self.feat], hits)))))

    def finalize_bundle(self, bundle) -> None:
        if bundle.feat_count > 1:
            self.stats_counts['_ambiguous_alignments_features'] += 1
            self.stats_counts['_ambiguous_reads_features'] += bundle.dup_counts
        else:
            self.stats_counts['_alignments_unique_features'] += 1
            self.stats_counts['_reads_unique_features'] += 1

    def write_intermediate_file(self, prefix):
        with open(f"{self.out_prefix}_{prefix}_out_aln_table.txt", 'w') as imf:
            imf.writelines(
                # sequence, cor_counts, strand, start, end, feat1a/feat1b;feat2;...
                map(lambda rec: "%s\t%f\t%c\t%d\t%d\t%s\n" % rec, self.alignments)
            )

    def write_report_files(self, alias, prefix=None):
        if prefix is None: prefix = self.out_prefix
        self.write_summary_statistics(prefix)
        self.write_feat_counts(alias, prefix)
        self.write_nt_len_mat(prefix)

    def write_summary_statistics(self, prefix):
        with open(prefix + '_stats.txt', 'w') as out:
            out.write('Summary Statistics\n')
            out.writelines("%s\t%d\n" % (key, val) for key, val in self.stats_counts.items())

    def write_feat_counts(self, alias, prefix):
        """Writes selected features and their associated counts to `prefix`_out_feature_counts.txt

        If a features.csv rule defined an ID Attribute other than "ID", then the associated features will
        be aliased to their corresponding ID Attribute value and the feature's true "ID" will be indicated
        next to it in parentheses. If an aliased or unaliased feature identifier is of list type, then the
        list of names will be joined on forward slash.

        For example, if the rule contained an ID Attribute which aliases gene1 to abc123,def456,123.456,
        then the feature column of the output file for this feature will read:
            abc123/def456/123.456 (gene1)
        """

        def list_and_alias(feat):
            if feat in alias:
                return '/'.join(alias[feat]) + f" ({feat[0]})"
            else:
                return '/'.join(feat)

        feat_counts_df = pd.DataFrame.from_dict(self.feat_counts, orient='index').reset_index()
        feat_counts_df['index'] = feat_counts_df['index'].apply(list_and_alias)
        feat_counts_df.to_csv(prefix + '_out_feature_counts.txt', sep='\t', index=False, header=False)

    def write_nt_len_mat(self, prefix):
        pd.DataFrame(self.nt_len_mat).to_csv(prefix + '_out_nt_len_dist.csv')

    def get_mp_return(self) -> Tuple[Counter, dict, dict]:
        """Produces a tuple of only relevant stats, for efficiency of serialization upon MP return"""
        return self.feat_counts, self.stats_counts, self.nt_len_mat

    def merge(self, other: 'StatsCollector'):
        self.feat_counts.update(other.feat_counts)

        for stat, count in other.stats_counts.items():
            self.stats_counts[stat] += count

        for nt, counter in other.nt_len_mat.items():
            self.nt_len_mat[nt].update(counter)

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
