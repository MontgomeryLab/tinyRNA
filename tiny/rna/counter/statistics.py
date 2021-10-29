import pandas as pd
import json
import sys
import os

from typing import Tuple, Optional, Iterable
from collections import Counter, defaultdict

from .hts_parsing import Alignment
from ..util import make_filename


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

    def __init__(self, out_prefix: str = None, report_diags: bool = False, no_normalize = False, **kwargs):
        self.library = {'Name': 'Unassigned', 'File': 'Unassigned'}
        self.out_prefix = out_prefix
        self.diags = Diagnostics(out_prefix) if report_diags else None

        self.feat_counts = Counter()
        self.chrom_misses = Counter()
        self.nt_len_mat = {nt: Counter() for nt in ['A', 'T', 'G', 'C']}
        self.library_stats = {stat: 0 for stat in LibraryStats.summary_categories}

        if no_normalize:
            setattr(self, "norm_by_genomic_hits", lambda x, _: x)
            setattr(self, "norm_by_feature_hits", lambda x, _: x.corr_count)

    def assign_library(self, library: dict):
        self.library = library

    def norm_by_genomic_hits(self, read_count, loci_count):
        return read_count / loci_count

    def count_bundle(self, aln_bundle: iter) -> Bundle:
        """Called for each multiple-alignment bundle before it is processed"""

        bundle_read = aln_bundle[0].read
        loci_counts = len(aln_bundle)

        # Calculate counts for multi-mapping
        read_counts = int(bundle_read.name.split('=')[1])
        corr_counts = self.norm_by_genomic_hits(read_counts, loci_counts)

        # Fill in 5p nt/length matrix
        sequence = bundle_read.seq
        self.nt_len_mat[bundle_read.nt5][len(sequence)] += read_counts

        return self.Bundle(loci_counts, read_counts, corr_counts)

    def norm_by_feature_hits(self, bundle, assigned_count):
        return bundle.corr_count / assigned_count

    def count_bundle_alignments(self, bundle: Bundle, aln: Alignment, assignments: set, n_candidates: int) -> None:
        """Called for each alignment for each read"""

        assigned_count = len(assignments)

        if assigned_count == 0:
            self.library_stats['Total Unassigned Reads'] += bundle.corr_count
        if assigned_count == 1:
            self.library_stats['Reads Assigned to Single Feature'] += bundle.corr_count
        if assigned_count > 1:
            self.library_stats['Reads Assigned to Multiple Features'] += bundle.corr_count
        if assigned_count > 0:
            self.library_stats['Total Assigned Reads'] += bundle.corr_count

            feature_corrected_count = self.norm_by_feature_hits(bundle, assigned_count)
            bundle.feat_count += len(assignments)
            bundle.assignments |= assignments

            for feat in assignments:
                self.feat_counts[feat] += feature_corrected_count

        if self.diags is not None:
            self.diags.record_diagnostics(assignments, n_candidates, aln, bundle)
            self.diags.record_alignment_details(aln, bundle, assignments)

    def finalize_bundle(self, bundle: Bundle) -> None:
        """Called at the conclusion of processing each multiple-alignment bundle"""

        assignment_count = len(bundle.assignments)

        if assignment_count == 0:
            self.library_stats['Total Unassigned Sequences'] += 1
        else:
            self.library_stats['Total Assigned Sequences'] += 1

            self.library_stats['Sequences Assigned to Single Feature'] += 1 * (assignment_count == 1)
            self.library_stats['Sequences Assigned to Multiple Features'] += 1 * (bundle.feat_count > 1)
            self.library_stats['Assigned Single-Mapping Reads'] += bundle.read_count * (bundle.loci_count == 1)
            self.library_stats['Assigned Multi-Mapping Reads'] += bundle.read_count * (bundle.loci_count > 1)


class SummaryStats:

    summary_categories = ["Total Reads", "Retained Reads", "Unique Sequences",
                          "Mapped Sequences", "Mapped Reads", "Aligned Reads"]

    def __init__(self, feature_attributes: dict, out_prefix: str, report_diags: bool):
        self.feature_attributes = feature_attributes
        self.out_prefix = out_prefix
        self.report_diags = report_diags

        # Will become False if an added library lacks its corresponding Collapser and Bowtie outputs
        self.report_summary_statistics = True

        self.pipeline_stats_df = pd.DataFrame(index=SummaryStats.summary_categories)
        self.feat_counts_df = pd.DataFrame(index=feature_attributes.keys())
        self.lib_stats_df = pd.DataFrame(index=LibraryStats.summary_categories)
        self.chrom_misses = Counter()
        self.nt_len_mat = {}
        self.warnings = []

        if report_diags:
            self.aln_diags = pd.DataFrame(columns=Diagnostics.aln_diag_categories)
            self.selection_diags = {}

    def write_report_files(self, counts_idx: set, alias: dict) -> None:
        self.print_warnings()
        if self.report_diags:
            Diagnostics.write_summary(self.out_prefix, self.aln_diags, self.selection_diags)

        self.write_feat_counts(counts_idx, alias, self.out_prefix)
        self.write_alignment_statistics(self.out_prefix)
        self.write_pipeline_statistics(self.out_prefix)
        self.write_nt_len_mat(self.out_prefix)

    def write_alignment_statistics(self, prefix: str) -> None:
        # Sort columns by title and round all counts to 2 decimal places
        self.df_to_csv(self.lib_stats_df, "Alignment Statistics", prefix, 'alignment_stats')

    def write_pipeline_statistics(self, prefix: str) -> None:
        if self.report_summary_statistics:
            # Sort columns by title and round all counts to 2 decimal places
            self.df_to_csv(self.pipeline_stats_df, "Summary Statistics", prefix, "summary_stats")

    def write_feat_counts(self, display_index: Iterable, alias: dict, prefix: str) -> None:
        """Writes selected features and their associated counts to {prefix}_out_feature_counts.csv

        Only the features in display_index will be listed in the output table, regardless of count, even
        if the features had no associated alignments.

        If a features.csv rule defined a Name Attribute other than "ID", then the associated features will
        be aliased to their corresponding Name Attribute value and displayed in the Feature Name column, and
        the feature's true "ID" will be indicated in the Feature ID column. If multiple aliases exist for
        a feature then they will be joined by ", " in the Feature Name column. A Feature Class column also
        follows.

        For example, if the rule contained a Name Attribute which aliases gene1 to abc123,def456,123.456
        and is both ALG and CSR class, then the Feature ID, Feature Name, and Feature Class column of the
        output file for this feature will read:
            gene1, "abc123,def456,123.456", "ALG,CSR"

        Subsequent columns represent the counts from each library processed by Counter. Column titles are
        formatted for each library as {group}_rep_{replicate}
        """

        # Subset the counts table to only show features that were matched on identity (regardless of count)
        summary = self.feat_counts_df.loc[display_index]
        # Sort columns by title and round all counts to 2 decimal places
        summary = self.sort_cols_and_round(summary)
        # Add Feature Name column, which is the feature alias (default is Feature ID if no alias exists)
        summary.insert(0, "Feature Name", summary.index.map(lambda feat: ', '.join(alias.get(feat, [feat]))))
        # Add Classes column for classes associated with the given feature
        feat_class_map = lambda feat: ', '.join(self.feature_attributes[feat][0][1])
        summary.insert(1, "Feature Class", summary.index.map(feat_class_map))
        # Sort by index, make index its own column, and rename it to Feature ID
        summary = summary.sort_index().reset_index().rename(columns={"index": "Feature ID"})

        summary.to_csv(prefix + '_feature_counts.csv', index=False)

    def write_nt_len_mat(self, prefix: str) -> None:
        """Writes each library's 5' end nucleotide / length matrix to its own file."""

        for lib_name, matrix in self.nt_len_mat.items():
            sanitized_lib_name = lib_name.replace('/', '_')
            len_dist_df = pd.DataFrame(matrix).sort_index().fillna(0)
            len_dist_df.to_csv(f'{prefix}_{sanitized_lib_name}_nt_len_dist.csv')

    def add_library(self, other: LibraryStats) -> None:
        # Add incoming feature counts as a new column of the data frame
        # Since other.feat_counts is a Counter object, unrecorded features default to 0 on lookup
        self.feat_counts_df[other.library["Name"]] = self.feat_counts_df.index.map(other.feat_counts)
        self.lib_stats_df[other.library["Name"]] = self.lib_stats_df.index.map(other.library_stats)
        self.nt_len_mat[other.library["Name"]] = other.nt_len_mat
        self.chrom_misses.update(other.chrom_misses)

        if self.report_diags: self.add_diags(other)

        # Process pipeline step outputs for this library, if they exist, to provide Summary Statistics
        if self.report_summary_statistics and self.library_has_pipeline_outputs(other):
            self.add_summary_stats(other)

    def add_summary_stats(self, other: LibraryStats):
        """Add incoming summary stats as new column in the master table"""

        mapped_seqs = sum([other.library_stats[stat] for stat in
                           ["Total Assigned Sequences",
                            "Total Unassigned Sequences"]])
        mapped_reads = sum([other.library_stats[stat] for stat in
                            ['Assigned Multi-Mapping Reads',
                             'Assigned Single-Mapping Reads',
                             'Total Unassigned Reads']])
        aligned_reads = other.library_stats["Total Assigned Reads"]
        total_reads, retained_reads = self.get_fastp_stats(other)
        unique_seqs = self.get_collapser_stats(other)

        other_summary = {
            "Total Reads": total_reads,
            "Retained Reads": retained_reads,
            "Unique Sequences": unique_seqs,
            "Mapped Sequences": mapped_seqs,
            "Mapped Reads": mapped_reads,
            "Aligned Reads": aligned_reads
        }

        self.pipeline_stats_df[other.library["Name"]] = self.pipeline_stats_df.index.map(other_summary)

    def add_diags(self, other: LibraryStats) -> None:
        """Append alignment and selection diagnostics to the master tables"""

        other_lib = other.library['Name']
        self.selection_diags[other_lib] = other.diags.selection_diags
        self.aln_diags.loc[other_lib] = other.diags.alignment_diags

    def print_warnings(self):
        if len(self.chrom_misses) != 0:
            self.warnings.append("\nSome reads aligned to chromosomes not defined in your GFF files.")
            self.warnings.append("This may be due to a mismatch in formatting (e.g. roman vs. arabic numerals)")
            self.warnings.append("between your bowtie indexes and GFF files, or your bowtie indexes may")
            self.warnings.append("contain chromosomes not defined in your GFF files.")
            self.warnings.append("The chromosome names and their alignment counts are:")
            for chr in sorted(self.chrom_misses.keys()):
                self.warnings.append("\t" + chr + ": " + str(self.chrom_misses[chr]))
            self.warnings.append("\n")

        for warning in self.warnings:
            print(warning, file=sys.stderr)

    def library_has_pipeline_outputs(self, other: LibraryStats) -> bool:
        """Check working directory for pipeline outputs from previous steps"""

        sam_basename = os.path.splitext(os.path.basename(other.library['File']))[0]
        lib_basename = sam_basename.replace("_aligned_seqs", "")
        fastp_logfile = lib_basename + "_qc.json"
        collapsed_fa = lib_basename + "_collapsed.fa"

        if not os.path.isfile(fastp_logfile) or not os.path.isfile(collapsed_fa):
            self.warnings.append(f"Pipeline output for {lib_basename} not found. Summary Statistics were skipped.")
            self.report_summary_statistics = False
            return False
        else:
            other.library['fastp_log'] = fastp_logfile
            other.library['collapsed'] = collapsed_fa
            return True

    def get_fastp_stats(self, other: LibraryStats) -> Optional[Tuple[int, int]]:
        """Determine the total number of reads for this library, and the total number retained by fastp"""

        try:
            with open(other.library['fastp_log'], 'r') as f:
                fastp_summary = json.load(f)['summary']

            total_reads = fastp_summary['before_filtering']['total_reads']
            retained_reads = fastp_summary['after_filtering']['total_reads']
        except (KeyError, json.JSONDecodeError):
            self.warnings.append("Unable to parse fastp json logs for Summary Statistics.")
            self.warnings.append("Associated file: " + other.library['File'])
            return None

        return total_reads, retained_reads

    def get_collapser_stats(self, other: LibraryStats) -> Optional[int]:
        """Determine the total number of unique sequences (after quality filtering) in this library"""

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
            self.warnings.append("Unable to parse collapsed fasta associated with for Summary Statistics.")
            self.warnings.append("Associated file: " + other.library['File'])
            return None

    @staticmethod
    def df_to_csv(df: pd.DataFrame, idx_name:str, prefix: str, postfix: str = None, sort_axis="columns"):
        if postfix is None:
            postfix = '_'.join(map(str.lower, idx_name.split(' ')))

        out_df = SummaryStats.sort_cols_and_round(df, axis=sort_axis)
        out_df.index.name = idx_name

        file_name = make_filename([prefix, postfix])
        out_df.to_csv(file_name)

    @staticmethod
    def sort_cols_and_round(df: pd.DataFrame, axis="columns") -> pd.DataFrame:
        """Convenience function to sort columns by title and round all values to 2 decimal places"""
        return df.round(decimals=2).sort_index(axis=axis)


class Diagnostics:

    aln_diag_categories = ['Eliminated counts', 'No feature counts',
                           'Uncounted alignments (+)', 'Uncounted alignments (-)']

    complement = bytes.maketrans(b'ACGTacgt', b'TGCAtgca')

    def __init__(self, out_prefix: str):
        self.out_prefix = out_prefix
        self.alignment_diags = {stat: 0 for stat in Diagnostics.aln_diag_categories}
        self.selection_diags = defaultdict(Counter)
        self.alignments = []

    def record_alignment_details(self, aln, bundle, assignments):
        """Record detailed alignment info if user elects to save diagnostics info with the run

        This is called once per locus per read (every alignment) when the user elects to save
        diagnostics. The recorded information is later written to {library['Name']}_aln_table.txt
        after the entire SAM file has been processed."""

        # Perform reverse complement for anti-sense reads
        read = aln.read.seq \
            if aln.iv.strand == '+' \
            else aln.read.seq[::-1].translate(self.complement)

        # sequence, cor_counts, strand, start, end, feat1;feat2;feat3
        self.alignments.append((read, bundle.corr_count, aln.iv.strand, aln.iv.start, aln.iv.end,
                                ';'.join(assignments)))

    def write_intermediate_file(self, library_name: str) -> None:
        """Write all recorded alignment info for this library to its corresponding alignment table."""

        outfile = make_filename([self.out_prefix, library_name, 'aln_table'], ext='.txt')
        with open(outfile, 'w') as imf:
            imf.writelines(
                # sequence, cor_counts, strand, start, end, feat1a/feat1b;feat2;...
                map(lambda rec: "%s\t%f\t%c\t%d\t%d\t%s\n" % rec, self.alignments)
            )

    def record_diagnostics(self, assignments, n_candidates, aln, bundle):
        """Records basic diagnostic info"""

        if len(assignments) == 0:
            if aln.iv.strand == '+':
                self.alignment_diags['Uncounted alignments (+)'] += 1
            else:
                self.alignment_diags['Uncounted alignments (-)'] += 1
            if n_candidates == 0:
                self.alignment_diags['No feature counts'] += bundle.corr_count
            else:
                self.alignment_diags['Eliminated counts'] += bundle.corr_count

    @classmethod
    def write_summary(cls, out_prefix: str, alignment_summary: pd.DataFrame, selection_summary: dict) -> None:
        """Write summary diagnostics, gathered by SummaryStatistics, for all libraries

        Args:
            out_prefix: prefix for the output diagnostic files
            alignment_summary: alignment diagnostics for ALL libraries
            selection_summary: selection diagnostics for ALL libraries
        """

        SummaryStats.df_to_csv(alignment_summary, "Sample", out_prefix, "alignment_diags", sort_axis="index")

        out = []
        for lib in sorted(selection_summary.keys()):
            out.append(lib)
            for feat_class in sorted(selection_summary[lib].keys()):
                out.append('\t' + feat_class)
                for stat in sorted(selection_summary[lib][feat_class].keys()):
                    out.append("\t\t%s: %d" % (stat, selection_summary[lib][feat_class][stat]))

        selection_summary_filename = make_filename([out_prefix, 'selection_diags'], ext='.txt')
        with open(selection_summary_filename, 'w') as f:
            f.write('\n'.join(out))