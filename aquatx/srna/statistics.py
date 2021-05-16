import pandas as pd
import json
import sys
import os

from typing import Tuple, Union
from collections import Counter

from aquatx.srna.FeatureSelector import FeatureSelector
from aquatx.srna.hts_parsing import Alignment


class LibraryStats:
    class Bundle:
        def __init__(self, loci_count, read_count, corr_count):
            self.loci_count = loci_count
            self.read_count = read_count
            self.corr_count = corr_count
            self.assignments = set()
            self.feat_count = 0

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

    def count_bundle(self, aln_bundle: iter) -> Bundle:
        bundle_read = aln_bundle[0].read
        loci_counts = len(aln_bundle)

        # Calculate counts for multi-mapping
        read_counts = int(bundle_read.name.split('=')[1])
        corr_counts = read_counts / loci_counts

        # Fill in 5p nt/length matrix
        sequence = bundle_read.seq
        self.nt_len_mat[bundle_read.nt5][len(sequence)] += read_counts

        return self.Bundle(loci_counts, read_counts, corr_counts)

    def count_bundle_alignments(self, bundle: Bundle, aln: Alignment, assignments: set) -> None:
        assigned_count = len(assignments)

        if assigned_count == 0:
            self.library_stats['Total Unassigned Reads'] += bundle.corr_count
        if assigned_count == 1:
            self.library_stats['Reads Assigned to Single Feature'] += bundle.corr_count
            # print(f"Single: {assignments}")
        if assigned_count > 1:
            self.library_stats['Reads Assigned to Multiple Features'] += bundle.corr_count
            # print(f"Multi: {assignments}")
        if assigned_count > 0:
            self.library_stats['Total Assigned Reads'] += bundle.corr_count

            feature_corrected_count = bundle.corr_count / assigned_count
            bundle.feat_count += len(assignments)
            bundle.assignments |= assignments

            for feat in assignments:
                self.feat_counts[feat] += feature_corrected_count

        # Todo: this will record antisense sequences as reverse complement
        #  Do we need intermediate files? Is it worth converting?
        if self.save_intermediate_file:
            # sequence, cor_counts, strand, start, end, feat1;feat2;feat3
            self.alignments.append((aln.read, bundle.corr_count, aln.iv.strand, aln.iv.start, aln.iv.end,
                                    ';'.join(assignments)))

    def finalize_bundle(self, bundle: Bundle) -> None:
        assignment_count = len(bundle.assignments)

        if assignment_count == 0:
            self.library_stats['Total Unassigned Sequences'] += 1
        else:
            self.library_stats['Total Assigned Sequences'] += 1

            self.library_stats['Sequences Assigned to Single Feature'] += 1 * (assignment_count == 1)
            self.library_stats['Sequences Assigned to Multiple Features'] += 1 * (bundle.feat_count > 1)
            self.library_stats['Assigned Single-Mapping Reads'] += bundle.read_count * (bundle.loci_count == 1)
            self.library_stats['Assigned Multi-Mapping Reads'] += bundle.read_count * (bundle.loci_count > 1)

    def write_intermediate_file(self) -> None:
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
        self.report_pipeline_stats = self.is_pipeline_invocation()
        self.nt_len_mat = {}

    def is_pipeline_invocation(self) -> bool:
        """Check working directory for pipeline outputs from previous steps"""

        for library in self.libraries:
            sam_basename = os.path.splitext(os.path.basename(library['File']))[0]
            lib_basename = sam_basename.replace("_aligned_seqs", "")
            fastp_logfile = lib_basename + "_qc.json"
            collapsed_fa = lib_basename + "_collapsed.fa"

            if not os.path.isfile(fastp_logfile) or not os.path.isfile(collapsed_fa):
                print(f"Pipeline output for {lib_basename} not found. Skipping Summary Statistics.")
                return False
            else:
                library['fastp_log'] = fastp_logfile
                library['collapsed'] = collapsed_fa

        return True

    def write_report_files(self, alias: dict, prefix=None) -> None:
        if prefix is None: prefix = self.out_prefix
        self.write_alignment_statistics(prefix)
        self.write_pipeline_statistics(prefix)
        self.write_feat_counts(alias, prefix)
        self.write_nt_len_mat(prefix)

    def write_alignment_statistics(self, prefix: str) -> None:
        # Sort columns by title and round all counts to 2 decimal places
        self.lib_stats_df = self.sort_cols_and_round(self.lib_stats_df)
        self.lib_stats_df.index.name = "Alignment Statistics"
        self.lib_stats_df.to_csv(prefix + '_alignment_stats.csv')

    def write_pipeline_statistics(self, prefix: str) -> None:
        if self.report_pipeline_stats:
            # Sort columns by title and round all counts to 2 decimal places
            self.pipeline_stats_df = self.sort_cols_and_round(self.pipeline_stats_df)
            self.pipeline_stats_df.index.name = "Summary Statistics"
            self.pipeline_stats_df.to_csv(prefix + '_summary_stats.csv')

    def write_feat_counts(self, alias: dict, prefix: str) -> None:
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

    def write_nt_len_mat(self, prefix: str) -> None:
        for lib_name, matrix in self.nt_len_mat.items():
            sanitized_lib_name = lib_name.replace('/', '_')
            pd.DataFrame(matrix).fillna(0).to_csv(f'{prefix}_{sanitized_lib_name}_nt_len_dist.csv')

    def add_library(self, other: LibraryStats) -> None:
        # Add incoming feature counts as a new column of the data frame
        # Since other.feat_counts is a Counter object, unrecorded features default to 0 on lookup
        self.feat_counts_df[other.library["Name"]] = self.feat_counts_df.index.map(other.feat_counts)
        self.lib_stats_df[other.library["Name"]] = self.lib_stats_df.index.map(other.library_stats)
        self.nt_len_mat[other.library["Name"]] = other.nt_len_mat

        if self.report_pipeline_stats:
            mapped_seqs = other.library_stats["Total Assigned Sequences"] + other.library_stats["Total Unassigned Sequences"]
            aligned_reads = other.library_stats["Total Assigned Reads"]
            total_reads, retained_reads = self.get_fastp_stats(other)
            unique_seqs = self.get_collapser_stats(other)

            other_summary = {
                "Total Reads": total_reads,
                "Retained Reads": retained_reads,
                "Unique Sequences": unique_seqs,
                "Mapped Sequences": mapped_seqs,
                "Aligned Reads": aligned_reads
            }

            self.pipeline_stats_df[other.library["Name"]] = self.pipeline_stats_df.index.map(other_summary)

    @staticmethod
    def get_fastp_stats(other: LibraryStats) -> Tuple[Union[int,str], Union[int,str]]:
        try:
            with open(other.library['fastp_log'], 'r') as f:
                fastp_summary = json.load(f)['summary']

            total_reads = fastp_summary['before_filtering']['total_reads']
            retained_reads = fastp_summary['after_filtering']['total_reads']
        except (KeyError, json.JSONDecodeError):
            print("Unable to parse fastp json logs for Summary Statistics.", file=sys.stderr)
            print("Associated file: " + other.library['File'], file=sys.stderr)
            return "err", "err"

        return total_reads, retained_reads

    @staticmethod
    def get_collapser_stats(other: LibraryStats) -> Union[int, str]:
        with open(other.library['collapsed'], 'r') as f:
            # Get file size and seek to 75 bytes from end
            size = f.seek(0, 2)
            f.seek(size - 75)

            # Read the last 75 bytes of the collapsed fasta
            tail = f.read(75)

        # Parse unique sequence count from the final fasta header
        from_pos = tail.rfind(">") + 1
        to_pos = tail.rfind("_count=")

        try:
            return int(tail[from_pos:to_pos])
        except ValueError:
            print("Unable to parse collapsed fasta associated with for Summary Statistics.", file=sys.stderr)
            print("Associated file: " + other.library['File'], file=sys.stderr)
            return "err"

    @staticmethod
    def sort_cols_and_round(df: pd.DataFrame) -> pd.DataFrame:
        sorted_columns = sorted(df.columns)
        return df.round(decimals=2).reindex(sorted_columns, axis="columns")
