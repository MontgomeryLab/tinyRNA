import pandas as pd
import mmap
import json
import csv
import sys
import os

from typing import Tuple, Optional, Iterable
from collections import Counter, defaultdict

from ..util import make_filename


class LibraryStats:

    summary_categories = ['Total Assigned Reads', 'Total Unassigned Reads',
                          'Total Assigned Sequences', 'Total Unassigned Sequences',
                          'Assigned Single-Mapping Reads', 'Assigned Multi-Mapping Reads',
                          'Reads Assigned to Single Feature', 'Sequences Assigned to Single Feature',
                          'Reads Assigned to Multiple Features', 'Sequences Assigned to Multiple Features']

    def __init__(self, out_prefix: str = None, report_diags: bool = False, normalize=True, **kwargs):
        self.library = {'Name': 'Unassigned', 'File': 'Unassigned'}
        self.out_prefix = out_prefix
        self.diags = Diagnostics(out_prefix) if report_diags else None
        self.norm = normalize

        self.feat_counts = Counter()
        self.ident_counts = Counter()
        self.chrom_misses = Counter()
        self.mapped_nt_len = {nt: Counter() for nt in ['A', 'T', 'G', 'C']}
        self.assigned_nt_len = {nt: Counter() for nt in ['A', 'T', 'G', 'C']}
        self.library_stats = {stat: 0 for stat in LibraryStats.summary_categories}

    def assign_library(self, library: dict):
        self.library = library

    def count_bundle(self, aln_bundle: iter) -> dict:
        """Called for each multiple-alignment bundle before it is processed"""

        bundle_read = aln_bundle[0]
        loci_counts = len(aln_bundle)
        nt5, seqlen = bundle_read['nt5'], len(bundle_read['seq'])

        # Calculate counts for multi-mapping
        read_counts = int(bundle_read['name'].split('=')[1])
        corr_counts = read_counts / loci_counts

        # Fill in 5p nt/length matrix
        self.mapped_nt_len[nt5][seqlen] += read_counts

        return {
            'loci_count': loci_counts,
            'read_count': read_counts,
            'corr_count': corr_counts,
            'assigned_feats': set(),
            'assigned_reads': 0,
            'nt5_counts': self.assigned_nt_len[nt5],
            'seq_length': seqlen,
            'mapping_stat':
                "Assigned Single-Mapping Reads"
                if loci_counts == 1 else
                "Assigned Multi-Mapping Reads"
        }

    def count_bundle_assignments(self, bundle: dict, aln: dict, assignments: dict, n_candidates: int) -> None:
        """Called for each alignment for each read"""

        feat_count = len(assignments)
        corr_count = bundle['corr_count']

        if feat_count == 0:
            self.library_stats['Total Unassigned Reads'] += corr_count
        else:
            fcorr_count = corr_count / feat_count if self.norm else corr_count
            bundle['assigned_reads'] += fcorr_count * feat_count
            bundle['assigned_feats'] |= assignments.keys()

            for feat, matched_rules in assignments.items():
                self.feat_counts[feat] += fcorr_count
                rcorr_count = fcorr_count / len(matched_rules)
                for rule in matched_rules:
                    self.ident_counts[rule] += rcorr_count

        if self.diags is not None:
            self.diags.record_diagnostics(assignments.keys(), n_candidates, aln, bundle)
            self.diags.record_alignment_details(aln, bundle, assignments.keys())

    def finalize_bundle(self, bundle: dict) -> None:
        """Called at the conclusion of processing each multiple-alignment bundle"""

        assigned_feat_count = len(bundle['assigned_feats'])

        if assigned_feat_count == 0:
            self.library_stats['Total Unassigned Sequences'] += 1
        else:
            self.library_stats['Total Assigned Sequences'] += 1
            assigned_reads = bundle['assigned_reads']

            self.library_stats['Total Assigned Reads'] += assigned_reads
            self.library_stats[bundle['mapping_stat']] += assigned_reads
            bundle['nt5_counts'][bundle['seq_length']] += assigned_reads

            if assigned_feat_count == 1:
                self.library_stats['Reads Assigned to Single Feature'] += assigned_reads
                self.library_stats['Sequences Assigned to Single Feature'] += 1
            else:
                self.library_stats['Reads Assigned to Multiple Features'] += assigned_reads
                self.library_stats['Sequences Assigned to Multiple Features'] += 1


# Todo: this class is messy. It needs to be refactored into purposed classes stored in a list
#  and iterated over with minimum interface:
#  - add_library()
#  - write_output_logfile()

class SummaryStats:

    summary_categories = ["Total Reads", "Retained Reads", "Unique Sequences",
                          "Mapped Sequences", "Mapped Reads", "Assigned Reads"]

    def __init__(self, classes, prefs):
        self.feature_classes = classes
        self.out_prefix = prefs.get('out_prefix', '')
        self.report_diags = prefs.get('report_diags', False)
        self.features_csv = prefs['config']

        # Will become False if an added library lacks its corresponding Collapser and Bowtie outputs
        self.report_summary_statistics = True

        self.pipeline_stats_df = pd.DataFrame(index=SummaryStats.summary_categories)
        self.feat_counts_df = pd.DataFrame(index=self.feature_classes.keys())
        self.lib_stats_df = pd.DataFrame(index=LibraryStats.summary_categories)
        self.rule_counts_df = pd.DataFrame()
        self.chrom_misses = Counter()
        self.nt_len_matrices = {}
        self.warnings = []

        if self.report_diags:
            self.aln_diags = pd.DataFrame(columns=Diagnostics.aln_diag_categories)
            self.selection_diags = {}

    def write_report_files(self, alias: dict) -> None:
        self.print_warnings()
        if self.report_diags:
            Diagnostics.write_summary(self.out_prefix, self.aln_diags, self.selection_diags)

        self.write_feat_counts(alias, self.out_prefix)
        self.write_rule_counts(self.out_prefix)
        self.write_alignment_statistics(self.out_prefix)
        self.write_pipeline_statistics(self.out_prefix)
        self.write_nt_len_matrices(self.out_prefix)

    def write_alignment_statistics(self, prefix: str) -> None:
        # Sort columns by title and round all counts to 2 decimal places
        self.df_to_csv(self.lib_stats_df, "Alignment Statistics", prefix, 'alignment_stats')

    def write_pipeline_statistics(self, prefix: str) -> None:
        if self.report_summary_statistics:
            # Sort columns by title and round all counts to 2 decimal places
            self.df_to_csv(self.pipeline_stats_df, "Summary Statistics", prefix, "summary_stats")

    def write_feat_counts(self, alias: dict, prefix: str) -> None:
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

        summary = self.feat_counts_df
        # Sort columns by title and round all counts to 2 decimal places
        summary = self.sort_cols_and_round(summary)
        # Add Feature Name column, which is the feature alias (default is Feature ID if no alias exists)
        summary.insert(0, "Feature Name", summary.index.map(lambda feat: ', '.join(alias.get(feat, ''))))
        # Add Classes column for classes associated with the given feature
        feat_class_map = lambda feat: ', '.join(self.feature_classes[feat])
        summary.insert(1, "Feature Class", summary.index.map(feat_class_map))
        # Sort by index, make index its own column, and rename it to Feature ID
        summary = summary.sort_index().reset_index().rename(columns={"index": "Feature ID"})

        summary.to_csv(prefix + '_feature_counts.csv', index=False)

    def write_nt_len_matrices(self, prefix: str) -> None:
        """Writes each library's 5' end nucleotide / length matrices to their own files."""

        for lib_name, matrices in self.nt_len_matrices.items():
            sanitized_lib_name = lib_name.replace('/', '_')

            mapped_nt_len_df = pd.DataFrame(matrices[0]).sort_index().fillna(0)
            mapped_nt_len_df.to_csv(f'{prefix}_{sanitized_lib_name}_mapped_nt_len_dist.csv')

            assigned_nt_len_df = pd.DataFrame(matrices[1]).sort_index().round(decimals=2).fillna(0)
            assigned_nt_len_df.to_csv(f'{prefix}_{sanitized_lib_name}_assigned_nt_len_dist.csv')

    def write_rule_counts(self, prefix: str) -> None:
        # Reread the Features Sheet since FeatureSelector.rules_table is processed and less readable
        with open(self.features_csv, 'r', encoding='utf-8-sig') as f:
            # Convert each CSV row to a string with values labeled by column headers
            rules = ['; '.join([': '.join([c,v]) for c,v in row.items()]) for row in csv.DictReader(f)]

        self.rule_counts_df = self.rule_counts_df.join(pd.Series(rules, name="Rule"), how="outer")
        self.rule_counts_df = self.sort_cols_and_round(self.rule_counts_df).set_index('Rule').fillna(0)
        self.rule_counts_df.loc["Mapped Reads"] = mapped_reads = self.pipeline_stats_df.loc["Mapped Reads"]

        if mapped_reads.empty:
            # Row will be empty if Summary Stats were not gathered
            # Todo: Mapped Reads can be calculated without pipeline outputs
            self.rule_counts_df.drop("Mapped Reads", axis=0, inplace=True)

        self.df_to_csv(self.rule_counts_df, "Rule", prefix, 'counts_by_rule')

    def add_library(self, other: LibraryStats) -> None:
        name = other.library["Name"]
        # Index varies per library -> join
        self.rule_counts_df = self.rule_counts_df.join(pd.Series(other.ident_counts, name=name), how='outer')
        # Index is consistent across libraries -> index.map
        self.feat_counts_df[name] = self.feat_counts_df.index.map(other.feat_counts)
        self.lib_stats_df[name] = self.lib_stats_df.index.map(other.library_stats)
        self.nt_len_matrices[name] = [other.mapped_nt_len, other.assigned_nt_len]
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
                            ['Reads Assigned to Multiple Features',
                             'Reads Assigned to Single Feature',
                             'Total Unassigned Reads']])
        assigned_reads = other.library_stats["Total Assigned Reads"]
        total_reads, retained_reads = self.get_fastp_stats(other)
        unique_seqs = self.get_collapser_stats(other)

        other_summary = {
            "Total Reads": total_reads,
            "Retained Reads": retained_reads,
            "Unique Sequences": unique_seqs,
            "Mapped Sequences": mapped_seqs,
            "Mapped Reads": mapped_reads,
            "Assigned Reads": assigned_reads
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
            with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                from_pos = mm.rfind(b">") + 1
                to_pos = mm.rfind(b"_count=")
                count = mm[from_pos:to_pos]

        try:
            # Zero-based count
            return int(count) + 1
        except ValueError:
            self.warnings.append(f"Unable to parse {other.library['collapsed']} for Summary Statistics.")
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
        read = aln['seq'] \
            if aln['strand'] == '+' \
            else aln['seq'][::-1].translate(self.complement)

        # sequence, cor_counts, strand, start, end, feat1;feat2;feat3
        self.alignments.append((read, bundle['corr_count'], aln['strand'], aln['start'], aln['end'],
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
            if aln['strand'] == '+':
                self.alignment_diags['Uncounted alignments (+)'] += 1
            else:
                self.alignment_diags['Uncounted alignments (-)'] += 1
            if n_candidates == 0:
                self.alignment_diags['No feature counts'] += bundle['corr_count']
            else:
                self.alignment_diags['Eliminated counts'] += bundle['corr_count']

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
                out.append('\t' + ','.join(feat_class))
                for stat in sorted(selection_summary[lib][feat_class].keys()):
                    out.append("\t\t%s: %d" % (stat, selection_summary[lib][feat_class][stat]))

        selection_summary_filename = make_filename([out_prefix, 'selection_diags'], ext='.txt')
        with open(selection_summary_filename, 'w') as f:
            f.write('\n'.join(out))