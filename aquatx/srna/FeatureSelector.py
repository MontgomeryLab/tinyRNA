import itertools
from collections import defaultdict, Counter
from typing import List, Tuple, FrozenSet, Dict, Set

import HTSeq
import pandas as pd
import json
import os
import re

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
    """

    rank, rule, feat = 0, 1, 2
    attributes = {}

    def __init__(self, rules: List[dict], reference_table: Dict):
        # Todo: now that filters are type-specific rather than column specific,
        #  we can specify filter type in the input rules list such that feature selection
        #  is no longer bound to this specific interest set.
        #  This will improve the project's maintainability.
        FeatureSelector.attributes = reference_table
        self.interest = ('Identity', 'Strand', 'nt5', 'Length')
        self.rules_table = sorted(rules, key=lambda x: x['Hierarchy'])
        self.build_filters()

        # Inverted ident rules: (Attrib Key, Attrib Val) as key, [rule matches] as val
        # This allows us to do O(1) identity matching for each candidate feature
        # Subsequent matching (also O(1)) is performed only against identity-matched rows
        inverted_identities = defaultdict(list)
        for i, rule in enumerate(self.rules_table):
            inverted_identities[rule['Identity']].append(i)
        self.inv_ident = dict(inverted_identities)

    def choose(self, feat_set, alignment) -> Tuple[set, set]:
        # Perform an efficient first-round
        finalists, uncounted = self.choose_identities(feat_set, alignment.iv)
        if not finalists: return finalists, uncounted

        strand = alignment.iv.strand
        nt5end = alignment.read.nt5
        length = len(alignment.read)

        # Strand, nt5, and Length filtering uses simpler logic
        choices, eliminated = set(), set()
        for step, read in zip(self.interest[1:], (strand, nt5end, length)):
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

        finalists, uncounted, identity_hits = set(), set(), list()
        start, end, features = 0, 1, 2
        key, value = 0, 1

        for iv_feats in feats_list:
            # Check for perfect interval match only once per feature_set/alignment tuple
            perfect_iv_match = self.is_perfect_iv_match(iv_feats[start], iv_feats[end], aln_iv)
            for feat in iv_feats[features]:
                for attrib in FeatureSelector.attributes[feat]:
                    # If multiple values are associated with the attribute key, create their tuple products
                    for feat_ident in itertools.product([attrib[key]], attrib[value]):
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
            uniq_feats = {hit[self.feat] for hit in identity_hits}

            # Scenarios:
            # 1. len(identity_hits) == len(uniq_ranks) == len(uniq_feats): unique features and ranks -- covered
            # 2. len(identity_hits) == len(uniq_ranks) != len(uniq_feats): unique ranks, 1 or more shared feat -- covered
            # 3. len(identity_hits) == len(uniq_feats) != len(uniq_ranks): unique feats, 1 or more shared rank -- covered
            # 4. len(i_h) != len(u_r) and len(i_h) != len(u_f) and len(u_r) != len(u_f): shared feats, shared ranks, different number
            # 5. len(i_h) != len(u_r) and len(i_h) != len(u_f) and len(u_r) == len(u_f): shared feats, shared ranks, same number
                # Scenarios 4 & 5: needs further investigation/testing

            if len(identity_hits) == len(uniq_ranks) == len(uniq_feats):
                finalists.add(min(identity_hits, key=lambda x: x[self.rank]))
            elif len(identity_hits) == len(uniq_ranks) != len(uniq_feats):
                min_feat = min(identity_hits, key=lambda x: x[self.rank])
                finalists.add(min_feat)
                finalists.update(hit for hit in identity_hits if hit[self.feat] == min_feat[self.feat])
            else:
                # Two or more hits share the same hierarchy.
                min_rank = min(uniq_ranks)
                finalists.update(hit for hit in identity_hits if hit[self.rank] == min_rank)

        ih_unique = {hit[self.feat] for hit in identity_hits}
        if len(finalists) < len(ih_unique): uncounted.add("Filtered")  # Todo: this logic is incorrect. Is relevant?
        if len(ih_unique) < len(feats_list): uncounted.add("Unknown")
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
        return cls.rank, cls.rule, cls.feat


class LibraryStats:
    class Bundle:
        def __init__(self, loci_count, read_count, corr_count):
            self.loci_count = loci_count
            self.read_count = read_count
            self.corr_count = corr_count
            self.assignments = set()
            self.feat_count = 0

    rank, rule, feat = FeatureSelector.get_hit_indexes()

    summary_categories = ['Total Assigned Reads', 'Total Assigned Sequences',
                          'Assigned Single-Mapping Reads', 'Assigned Multi-Mapping Reads',
                          'Reads Assigned to Single Feature', 'Sequences Assigned to Single Feature',
                          'Reads Assigned to Multiple Features', 'Sequences Assigned to Multiple Features',
                          'Total Unassigned Reads', 'Total Unassigned Sequences']

    def __init__(self, library: dict, out_prefix: str = None, save_intermediate_file: bool = False):
        self.library = library
        self.out_prefix = out_prefix
        self.save_intermediate_file = save_intermediate_file

        self.feat_counts = Counter()
        self.nt_len_mat = {nt: Counter() for nt in ['A', 'T', 'G', 'C']}
        self.library_stats = {stat: 0 for stat in LibraryStats.summary_categories}
        self.alignments = []

    def count_bundle(self, aln_bundle: iter) -> 'Bundle':
        bundle_read = aln_bundle[0].read
        loci_counts = len(aln_bundle)

        # Calculate counts for multi-mapping
        read_counts = int(bundle_read.name.split('=')[1])
        corr_counts = read_counts / loci_counts

        # Fill in 5p nt/length matrix
        sequence = bundle_read.seq
        self.nt_len_mat[bundle_read.nt5][len(sequence)] += read_counts

        return self.Bundle(loci_counts, read_counts, corr_counts)

    def count_bundle_alignments(self, bundle: 'Bundle', aln: HTSeq.SAM_Alignment, assignments: set) -> None:
        assigned_count = len(assignments)

        if assigned_count == 0:
            self.library_stats['Total Unassigned Reads'] += bundle.corr_count
        if assigned_count == 1:
            self.library_stats['Reads Assigned to Single Feature'] += bundle.corr_count
        if assigned_count > 1:
            self.library_stats['Reads Assigned to Multiple Features'] += bundle.corr_count
        if assigned_count > 0:
            self.library_stats['Total Assigned Reads'] += bundle.corr_count

            feature_corrected_count = bundle.corr_count / assigned_count
            bundle.feat_count += len(assignments)
            bundle.assignments |= assignments

            for hit in assignments:
                self.feat_counts[hit[self.feat]] += feature_corrected_count

        # Todo: this will record antisense sequences as reverse complement
        #  Do we need intermediate files? Is it worth converting?
        if self.save_intermediate_file:
            # sequence, cor_counts, strand, start, end, feat1;feat2;feat3
            self.alignments.append((aln.read, bundle.corr_count, aln.iv.strand, aln.iv.start, aln.iv.end,
                                    ';'.join(map(lambda x: x[self.feat], assignments))))

    def finalize_bundle(self, bundle) -> None:
        assignment_count = len(bundle.assignments)

        if assignment_count == 0:
            self.library_stats['Total Unassigned Sequences'] += 1
        else:
            self.library_stats['Total Assigned Sequences'] += 1

            self.library_stats['Sequences Assigned to Single Feature'] += 1 * (assignment_count == 1)
            self.library_stats['Sequences Assigned to Multiple Features'] += 1 * (bundle.feat_count > 1)
            self.library_stats['Assigned Single-Mapping Reads'] += bundle.read_count * (bundle.loci_count == 1)
            self.library_stats['Assigned Multi-Mapping Reads'] += bundle.read_count * (bundle.loci_count > 1)

    def write_intermediate_file(self):
        with open(f"{self.out_prefix}_{self.library['Name']}_aln_table.txt", 'w') as imf:
            imf.writelines(
                # sequence, cor_counts, strand, start, end, feat1a/feat1b;feat2;...
                map(lambda rec: "%s\t%f\t%c\t%d\t%d\t%s\n" % rec, self.alignments)
            )


class SummaryStats:

    summary_categories = ["Total Reads", "Retained Reads", "Unique Sequences", "Mapped Sequences", "Aligned Reads"]

    def __init__(self, libraries, out_prefix):
        self.out_prefix = out_prefix
        self.libraries = libraries

        self.pipeline_stats_df = pd.DataFrame(index=SummaryStats.summary_categories)
        self.feat_counts_df = pd.DataFrame(index=FeatureSelector.attributes.keys())
        self.lib_stats_df = pd.DataFrame(index=LibraryStats.summary_categories)
        self.nt_len_mat = {nt: Counter() for nt in ['A', 'T', 'G', 'C']}
        self.report_pipeline_stats = self.is_pipeline_invocation()

    def is_pipeline_invocation(self):
        """Check for pipeline outputs from previous steps in working directory"""

        for library in self.libraries:
            sam_basename = os.path.splitext(os.path.basename(library['File']))[0]
            lib_basename = sam_basename.replace("_aligned_seqs", "")
            fastp_logfile = lib_basename + "_qc.json"
            collapsed_fa = lib_basename + "_collapsed.fa"

            if not os.path.isfile(fastp_logfile) or not os.path.isfile(collapsed_fa):
                return False
            else:
                library['fastp_log'] = fastp_logfile
                library['collapsed'] = collapsed_fa

        return True

    def write_report_files(self, alias, prefix=None):
        if prefix is None: prefix = self.out_prefix
        self.write_alignment_statistics(prefix)
        self.write_pipeline_statistics(prefix)
        self.write_feat_counts(alias, prefix)
        self.write_nt_len_mat(prefix)

    def write_alignment_statistics(self, prefix):
        # Sort columns by title and round all counts to 2 decimal places
        self.lib_stats_df = self.sort_cols_and_round(self.lib_stats_df)
        self.lib_stats_df.to_csv(prefix + '_alignment_stats.csv')

    def write_pipeline_statistics(self, prefix):
        if self.report_pipeline_stats:
            # Sort columns by title and round all counts to 2 decimal places
            self.pipeline_stats_df = self.sort_cols_and_round(self.pipeline_stats_df)
            self.pipeline_stats_df.to_csv(prefix + '_summary_stats.csv')

    def write_feat_counts(self, alias, prefix):
        """Writes selected features and their associated counts to {prefix}_out_feature_counts.csv

        The resulting table will have the following columns:
            - Feature ID: the "ID" attribute of the

        If a features.csv rule defined an ID Attribute other than "ID", then the associated features will
        be aliased to their corresponding ID Attribute value and displayed in the Feature Name column, and
        the feature's true "ID" will be indicated in the Feature ID column. If multiple aliases exist for
        a feature then they will be joined by ", " in the Feature Name column.

        For example, if the rule contained an ID Attribute which aliases gene1 to abc123,def456,123.456,
        then the feature column of the output file for this feature will read:
            abc123, def456, 123.456 	gene1
        """

        # Sort columns by title and round all counts to 2 decimal places
        summary = self.sort_cols_and_round(self.feat_counts_df)
        # Add Feature Name column, which is the feature alias (default is Feature ID if no alias exists)
        summary.insert(0, "Feature Name", summary.index.map(lambda feat: ', '.join(alias.get(feat, feat))))
        # Add Classes column for classes associated with the given feature
        feat_class_map = lambda feat: ', '.join(FeatureSelector.attributes[feat][0][1])
        summary.insert(1, "Feature Class", summary.index.map(feat_class_map))
        # Sort by index, make index its own column, and rename it to Feature ID
        summary = summary.sort_index().reset_index().rename(columns={"index": "Feature ID"})

        summary.to_csv(prefix + '_feature_counts.csv', index=False)

    def write_nt_len_mat(self, prefix):
        pd.DataFrame(self.nt_len_mat).to_csv(prefix + '_nt_len_dist.csv')

    def add_library(self, other: 'LibraryStats'):
        # Add incoming feature counts as a new column of the data frame
        # Since other.feat_counts is a Counter object, unrecorded features default to 0 on lookup
        self.feat_counts_df[other.library["Name"]] = self.feat_counts_df.index.map(other.feat_counts)
        self.lib_stats_df[other.library["Name"]] = self.lib_stats_df.index.map(other.library_stats)

        for nt, counter in other.nt_len_mat.items():
            self.nt_len_mat[nt].update(counter)

        if self.report_pipeline_stats:
            total_reads, retained_reads = self.get_fastp_stats(other.library['fastp_log'])
            unique_sequences = self.get_collapser_stats(other.library['collapsed'])
            mapped_seqs = other.library_stats["Total Assigned Sequences"] + other.library_stats["Total Unassigned Sequences"]
            aligned_reads = other.library_stats["Total Assigned Reads"]

            other_summary = {
                "Total Reads": total_reads,
                "Retained Reads": retained_reads,
                "Unique Sequences": unique_sequences,
                "Mapped Sequences": mapped_seqs,
                "Aligned Reads": aligned_reads
            }

            self.pipeline_stats_df[other.library["Name"]] = self.pipeline_stats_df.index.map(other_summary)

    @staticmethod
    def get_fastp_stats(logfile):
        with open(logfile, 'r') as f:
            fastp_summary = json.load(f)['summary']

        total_reads = fastp_summary['before_filtering']['total_reads']
        retained_reads = fastp_summary['after_filtering']['total_reads']

        return total_reads, retained_reads

    @staticmethod
    def get_collapser_stats(collapsed):
        with open(collapsed, 'r') as f:
            # Get file size and seek to 75 bytes from end
            size = f.seek(0, 2)
            f.seek(size - 75)

            # Read the last 75 bytes of the collapsed fasta
            tail = f.read(75)

        # Parse unique sequence count from the final fasta header
        from_pos = tail.rfind(">") + 1
        to_pos = tail.rfind("_count=")
        return tail[from_pos:to_pos]

    @staticmethod
    def sort_cols_and_round(df):
        sorted_columns = sorted(df.columns)
        return df.round(decimals=2).reindex(sorted_columns, axis="columns")