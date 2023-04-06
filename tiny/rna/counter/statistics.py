import pandas as pd
import gzip
import json
import mmap
import csv
import sys
import os
import re

from abc import abstractmethod, ABC
from typing import Tuple, Optional, Union
from collections import Counter, defaultdict
from glob import glob

from ..util import make_filename


class LibraryStats:

    summary_categories = ['Total Assigned Reads', 'Total Unassigned Reads',
                          'Total Assigned Sequences', 'Total Unassigned Sequences',
                          'Assigned Single-Mapping Reads', 'Assigned Multi-Mapping Reads',
                          'Reads Assigned to Single Feature', 'Sequences Assigned to Single Feature',
                          'Reads Assigned to Multiple Features', 'Sequences Assigned to Multiple Features']

    def __init__(self, out_prefix: str = None, report_diags: bool = False, **prefs):
        self.library = {'Name': 'Unassigned', 'File': 'Unassigned', 'Norm': '1'}
        self.out_prefix = out_prefix
        self.diags = Diagnostics(out_prefix) if report_diags else None
        self.norm_gh = prefs.get('normalize_by_genomic_hits', True)
        self.norm_fh = prefs.get('normalize_by_feature_hits', True)

        self.feat_counts = Counter()
        self.rule_counts = Counter()
        self.chrom_misses = Counter()
        self.mapped_nt_len = defaultdict(Counter, {nt: Counter() for nt in ['A', 'T', 'G', 'C']})
        self.assigned_nt_len = defaultdict(Counter, {nt: Counter() for nt in ['A', 'T', 'G', 'C']})
        self.library_stats = {stat: 0 for stat in LibraryStats.summary_categories}

    def assign_library(self, library: dict):
        self.library = library

    def count_bundle(self, aln_bundle: iter, read_counts: int) -> dict:
        """Called for each multiple-alignment bundle before it is processed"""

        bundle_read = aln_bundle[0]
        loci_counts = len(aln_bundle)
        corr_counts = read_counts / loci_counts if self.norm_gh else read_counts
        nt5, seqlen = bundle_read['nt5end'], bundle_read['Length']

        # Fill in 5p nt/length matrix
        self.mapped_nt_len[nt5][seqlen] += read_counts

        return {
            'read_count': read_counts,
            'corr_count': corr_counts,
            'assigned_ftags': set(),
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

        asgn_count = len(assignments)
        corr_count = bundle['corr_count']

        if asgn_count == 0:
            self.library_stats['Total Unassigned Reads'] += corr_count
        else:
            fcorr_count = corr_count / asgn_count if self.norm_fh else corr_count
            bundle['assigned_reads'] += fcorr_count * asgn_count
            bundle['assigned_ftags'] |= assignments.keys()

            for ftag, matched_rules in assignments.items():
                self.feat_counts[ftag] += fcorr_count
                rcorr_count = fcorr_count / len(matched_rules)
                for rule in matched_rules:
                    self.rule_counts[rule] += rcorr_count

        if self.diags is not None:
            self.diags.record_diagnostics(assignments.keys(), n_candidates, aln, bundle)
            self.diags.record_alignment_details(aln, bundle, assignments.keys())

    def finalize_bundle(self, bundle: dict) -> None:
        """Called at the conclusion of processing each multiple-alignment bundle"""

        assigned_feat_count = len({feat[0] for feat in bundle['assigned_ftags']})

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


class MergedStat(ABC):
    prefix = None

    @abstractmethod
    def add_library(self, other: LibraryStats): ...

    @abstractmethod
    def finalize(self): ...

    @abstractmethod
    def write_output_logfile(self): ...

    @staticmethod
    def add_warning(msg):
        MergedStatsManager.warnings.append(msg)

    @staticmethod
    def sort_cols_and_round(df: pd.DataFrame, axis="columns") -> pd.DataFrame:
        """Convenience function to sort columns by title and round all values to 2 decimal places"""
        return df.round(decimals=2).sort_index(axis=axis)

    @staticmethod
    def df_to_csv(df: pd.DataFrame, idx_name: str, prefix: str, postfix: str = None, sort_axis="columns"):
        if postfix is None:
            postfix = '_'.join(map(str.lower, idx_name.split(' ')))

        if sort_axis is not None:
            out_df = MergedStat.sort_cols_and_round(df, axis=sort_axis)
        else:
            out_df = df.round(decimals=2)

        out_df.index.name = idx_name
        file_name = make_filename([prefix, postfix])
        out_df.to_csv(file_name)


class MergedStatsManager:
    chrom_misses = Counter()
    warnings = []

    def __init__(self, Features_obj, features_csv, prefs):
        MergedStat.prefix = prefs.get('out_prefix', '')

        self.merged_stats = [
            FeatureCounts(Features_obj), RuleCounts(features_csv),
            AlignmentStats(), SummaryStats(),
            NtLenMatrices()
        ]

        if prefs.get('report_diags', False):
            self.merged_stats.append(MergedDiags())

    def write_report_files(self) -> None:
        try:
            for in_process in self.merged_stats:
                in_process.finalize()

            # Validate stats now that all libraries have been processed
            StatisticsValidator(self.merged_stats).validate()
        finally:
            self.print_warnings()

        for output_stat in self.merged_stats:
            output_stat.write_output_logfile()
        self.print_warnings()

    def add_library(self, other: LibraryStats) -> None:
        for merger in self.merged_stats:
            merger.add_library(other)
        self.chrom_misses.update(other.chrom_misses)

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


class FeatureCounts(MergedStat):
    def __init__(self, Features_obj):
        self.feat_counts_df = pd.DataFrame(index=set.union(*Features_obj.classes.values()))
        self.aliases = Features_obj.aliases
        self.finalized = False
        self.norm_prefs = {}

    def add_library(self, other: LibraryStats) -> None:
        assert not self.finalized, "Cannot add libraries after FeatureCounts object has been finalized."
        self.feat_counts_df[other.library["Name"]] = self.feat_counts_df.index.map(other.feat_counts)
        self.norm_prefs[other.library["Name"]] = other.library['Norm']

    def write_output_logfile(self) -> None:
        self.write_feature_counts()
        self.write_norm_counts()

    def finalize(self):
        """Called once all libraries have been counted and added to the FeatureCounts object.

        This method will:
        - Name the multi-index columns: "Feature ID" and "Classifier"
        - Add a Feature Name column containing the values associated with user-defined alias tags
        - Sort columns by name and round decimal values to 2 places.

        Aliases are concatenated by ", " if multiple exist for a feature. Alias values are not
        added for the "ID" tag as the Feature ID column already contains this information."""

        # Sort columns by title and round all counts to 2 decimal places
        summary = self.sort_cols_and_round(self.feat_counts_df)
        # Add Feature Name column, which is the feature alias (default is blank if no alias exists)
        summary.insert(0, "Feature Name", summary.index.map(lambda feat: ', '.join(self.aliases.get(feat[0], ''))))
        # Sort by index, make index its own column, and rename it to Feature ID
        summary = summary.sort_index().reset_index().rename(columns={"level_0": "Feature ID", "level_1": "Classifier"})

        self.feat_counts_df = summary
        self.finalized = True

    def write_feature_counts(self):
        """Writes selected features and their associated counts to {prefix}_feature_counts.csv
        Should be called after finalize()"""

        assert self.finalized, "FeatureCounts object must be finalized before writing output."
        self.feat_counts_df.to_csv(self.prefix + '_feature_counts.csv', index=False)

    def write_norm_counts(self):
        """Writes normalized feature counts to {prefix}_norm_counts.csv
        Normalization method is determined by Samples Sheet's Normalization column.
        Should be called after finalize()"""

        assert self.finalized, "FeatureCounts object must be finalized before writing output."

        if all([norm in ['1', "", None] for norm in self.norm_prefs.values()]):
            return

        counts = self.feat_counts_df.copy()
        for library, norm in self.norm_prefs.items():
            if norm in ['1', "", None]:
                continue
            elif re.match(r'^\s*\d+\s*$', norm):
                factor = int(norm)
            elif re.match(r'^\s*[\d.]+\s*$', norm):
                factor = float(norm)
            elif norm.upper().strip() == "RPM":
                mapped_reads = SummaryStats.pipeline_stats_df.loc['Mapped Reads', library]
                factor = mapped_reads / 1_000_000
            else:
                raise ValueError(f'Invalid normalization value for {library}: "{norm}"\n'
                                 "Please check your Samples Sheet.")

            counts[library] = counts[library] / factor

        counts.to_csv(self.prefix + '_norm_counts.csv', index=False)


class RuleCounts(MergedStat):
    def __init__(self, features_csv):
        self.rules = self.read_features_sheet(features_csv)
        self.rule_counts_df = pd.DataFrame()
        self.finalized = False

    def add_library(self, other: LibraryStats):
        assert not self.finalized, "Cannot add libraries after RuleCounts object has been finalized."
        counts = pd.Series(other.rule_counts, name=other.library['Name'])
        self.rule_counts_df = self.rule_counts_df.join(counts, how='outer')

    def finalize(self):
        rules = self.get_rule_strings()
        self.rule_counts_df = self.sort_cols_and_round(self.rule_counts_df)

        # Add Rule String column to the left of the counts
        order = ["Rule String"] + self.rule_counts_df.columns.to_list()
        self.rule_counts_df = self.rule_counts_df.join(pd.Series(rules, name="Rule String"), how="outer")[order].fillna(0)

        # Add Mapped Reads row below the counts
        self.rule_counts_df.loc["Mapped Reads"] = SummaryStats.pipeline_stats_df.loc["Mapped Reads"]
        self.finalized = True

    def write_output_logfile(self):
        assert self.finalized, "RuleCounts object must be finalized before writing output."
        self.df_to_csv(self.rule_counts_df, "Rule Index", self.prefix, 'counts_by_rule', sort_axis=None)

    def read_features_sheet(self, features_csv):
        """Reads the Features Sheet and returns a list of rules in
        the same order as the FeatureSelector.rules_table. We don't use the
        CSVReader class here because it returns rows with shortened column names,
        and we want to preserve the full column names for the Rule String column
        of the output table."""

        with open(features_csv, 'r', encoding='utf-8-sig') as f:
            return [row for row in csv.DictReader(f)]

    def get_inverted_classifiers(self):
        """Returns a dictionary of classifiers and the indices of the rules
        where they are defined. Note the use of the user-facing column name
        rather than the internal shortname for "Classify as..." since the
        CSVReader class wasn't used."""

        inverted_classifiers = defaultdict(list)
        for i, rule in enumerate(self.rules):
            inverted_classifiers[rule['Classify as...']].append(i)

        return dict(inverted_classifiers)

    def get_rule_strings(self):
        """Returns a list of formatted strings representing
        each rule in the Features Sheet."""

        rows = []
        for rule in self.rules:
            row = [': '.join([col.replace('...', ''), defn])
                   for col, defn in rule.items()]
            rows.append('; '.join(row))
        return rows


class NtLenMatrices(MergedStat):
    def __init__(self):
        self.nt_len_matrices = {}
        self.finalized = False

    def add_library(self, other: LibraryStats):
        assert not self.finalized, "Cannot add libraries after NtLenMatrices object has been finalized."

        name = other.library['Name']
        self.nt_len_matrices[name] = (
            other.mapped_nt_len,
            other.assigned_nt_len
        )

    def finalize(self):
        self.nt_len_matrices = dict(sorted(self.nt_len_matrices.items()))
        for lib_name, (mapped, assigned) in self.nt_len_matrices.items():
            mapped_nt_len_df = self.finalize_mapped(mapped)
            assigned_nt_len_df = self.finalize_assigned(assigned)
            self.nt_len_matrices[lib_name] = (mapped_nt_len_df, assigned_nt_len_df)

        self.finalized = True

    def finalize_mapped(self, mapped: DefaultDict[str, Counter]) -> pd.DataFrame:
        # Whole number counts because number of reads mapped is always integer
        return pd.DataFrame(mapped).sort_index().fillna(0)

    def finalize_assigned(self, assigned: DefaultDict[str, Counter]) -> pd.DataFrame:
        # Fractional counts due to (loci count) and/or (assigned feature count) normalization
        assigned_nt_len_df = pd.DataFrame(assigned).sort_index().round(decimals=2)

        # Drop non-nucleotide columns if they don't contain counts
        assigned_nt_len_df.drop([
            col for col, values in assigned_nt_len_df.iteritems()
            if col not in ['A', 'T', 'G', 'C'] and values.isna().all()
        ], axis='columns', inplace=True)

        # Assign default count of 0 for remaining cases
        return assigned_nt_len_df.fillna(0)

    def write_output_logfile(self):
        """Writes each library's 5' end nucleotide / length matrices to separate files."""

        assert self.finalized, "NtLenMatrices object must be finalized before writing output."

        for lib_name, (mapped, assigned) in self.nt_len_matrices.items():
            sanitized_lib_name = lib_name.replace('/', '_')
            mapped.to_csv(f'{self.prefix}_{sanitized_lib_name}_mapped_nt_len_dist.csv')
            assigned.to_csv(f'{self.prefix}_{sanitized_lib_name}_assigned_nt_len_dist.csv')


class AlignmentStats(MergedStat):
    def __init__(self):
        self.alignment_stats_df = pd.DataFrame(index=LibraryStats.summary_categories)
        self.finalized = False

    def add_library(self, other: LibraryStats):
        assert not self.finalized, "Cannot add libraries after AlignmentStats object has been finalized."
        library, stats = other.library['Name'], other.library_stats
        self.alignment_stats_df[library] = self.alignment_stats_df.index.map(stats)

    def finalize(self):
        self.finalized = True  # Nothing else to do

    def write_output_logfile(self):
        assert self.finalized, "AlignmentStats object must be finalized before writing output."
        self.df_to_csv(self.alignment_stats_df, "Alignment Statistics", self.prefix, 'alignment_stats')


class SummaryStats(MergedStat):

    constant_categories     = ["Mapped Sequences", "Mapped Reads", "Assigned Reads"]
    conditional_categories  = ["Total Reads", "Retained Reads", "Unique Sequences"]
    summary_categories      = conditional_categories + constant_categories

    pipeline_stats_df = pd.DataFrame(index=summary_categories)

    def __init__(self):
        self.missing_fastp_outputs = []
        self.missing_collapser_outputs = []
        self.finalized = False

    def add_library(self, other: LibraryStats):
        """Add incoming summary stats as new column in the master table"""

        assert not self.finalized, "Cannot add libraries after SummaryStats object has been finalized."

        other.library['basename'] = self.get_library_basename(other)
        other_summary = {
            "Mapped Sequences": self.get_mapped_seqs(other),
            "Mapped Reads": self.get_mapped_reads(other),
            "Assigned Reads": other.library_stats["Total Assigned Reads"]
        }

        if self.library_has_fastp_outputs(other):
            total_reads, retained_reads = self.get_fastp_stats(other)
            other_summary.update({
                "Total Reads": total_reads,
                "Retained Reads": retained_reads
            })

        if self.library_has_collapser_outputs(other):
            unique_seqs = self.get_collapser_stats(other)
            other_summary["Unique Sequences"] = unique_seqs

        self.pipeline_stats_df[other.library["Name"]] = self.pipeline_stats_df.index.map(other_summary)

    @staticmethod
    def get_mapped_seqs(other: LibraryStats):
        return sum([other.library_stats[stat] for stat in [
            "Total Assigned Sequences",
            "Total Unassigned Sequences"
        ]])

    @staticmethod
    def get_mapped_reads(other: LibraryStats):
        return sum([other.library_stats[stat] for stat in [
            'Reads Assigned to Multiple Features',
            'Reads Assigned to Single Feature',
            'Total Unassigned Reads'
        ]])

    def finalize(self):
        if len(self.missing_fastp_outputs):
            missing = '\n\t'.join(self.missing_fastp_outputs)
            self.add_warning("Total Reads and Retained Reads could not be determined for the following libraries "
                             "because their fastp outputs were not found in the working directory:\n\t" + missing)
        if len(self.missing_collapser_outputs):
            missing = '\n\t'.join(self.missing_collapser_outputs)
            self.add_warning("The Unique Sequences stat could not be determined for the following libraries because "
                             "their tiny-collapse outputs were not found in the working directory:\n\t" + missing)

        # Only display conditional categories if they were collected for at least one library
        empty_rows = self.pipeline_stats_df.loc[self.conditional_categories].isna().all(axis='columns')
        self.pipeline_stats_df.drop(empty_rows.index[empty_rows], inplace=True)
        self.pipeline_stats_df = self.sort_cols_and_round(self.pipeline_stats_df)
        self.finalized = True

    def write_output_logfile(self):
        assert self.finalized, "SummaryStats object must be finalized before writing output."
        self.df_to_csv(self.pipeline_stats_df, "Summary Statistics", self.prefix, "summary_stats")

    def library_has_collapser_outputs(self, other: LibraryStats) -> bool:
        # Collapser outputs may have been gzipped. Accept either filename.
        collapsed_fa = glob(other.library['basename'] + "_collapsed.fa*")
        if len(collapsed_fa):
            other.library['collapsed'] = collapsed_fa[0]
            return True
        else:
            self.missing_collapser_outputs.append(other.library['basename'])
            return False

    def library_has_fastp_outputs(self, other: LibraryStats) -> bool:
        fastp_logfile = other.library['basename'] + "_qc.json"
        if os.path.isfile(fastp_logfile):
            other.library['fastp_log'] = fastp_logfile
            return True
        else:
            self.missing_fastp_outputs.append(other.library['basename'])
            return False

    def get_fastp_stats(self, other: LibraryStats) -> Union[Tuple[int, int], Tuple[None, None]]:
        """Determine the total number of reads for this library, and the total number retained by fastp"""

        try:
            with open(other.library['fastp_log'], 'r') as f:
                fastp_summary = json.load(f)['summary']

            total_reads = fastp_summary['before_filtering']['total_reads']
            retained_reads = fastp_summary['after_filtering']['total_reads']
        except (KeyError, json.JSONDecodeError):
            self.add_warning(f"Unable to parse {other.library['fastp_log']} for Summary Statistics.")
            return None, None

        return total_reads, retained_reads

    def get_collapser_stats(self, other: LibraryStats) -> Optional[int]:
        """Determine the total number of unique sequences (after quality filtering) in this library"""

        collapsed_fa = other.library['collapsed']

        if collapsed_fa.endswith('.gz'):
            # Random access is much harder with gzip
            # Instead, get last header in last 250 bytes
            with gzip.open(collapsed_fa) as g:
                g.seek(-250, os.SEEK_END)
                lines = g.readlines()
                while lines[-1] == "": lines.pop(-1)
                last_header = lines[-2].decode()
                count = last_header.split('_count=')[0][1:]
        else:
            with open(collapsed_fa, 'r') as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                from_pos = mm.rfind(b">") + 1
                to_pos = mm.rfind(b"_count=")
                count = mm[from_pos:to_pos]

        try:
            # Zero-based count
            return int(count) + 1
        except ValueError:
            self.add_warning(f"Unable to parse {other.library['collapsed']} for Summary Statistics.")
            return None

    @staticmethod
    def get_library_basename(other: LibraryStats) -> str:
        sam_basename = os.path.splitext(os.path.basename(other.library['File']))[0]
        lib_basename = sam_basename.replace("_aligned_seqs", "")
        return lib_basename


class Diagnostics:

    aln_diag_categories = ['Eliminated counts', 'No feature counts',
                           'Uncounted alignments (+)', 'Uncounted alignments (-)']

    complement = bytes.maketrans(b'ACGTacgt', b'TGCAtgca')

    def __init__(self, out_prefix: str):
        self.prefix = out_prefix
        self.alignment_diags = {stat: 0 for stat in Diagnostics.aln_diag_categories}
        self.selection_diags = defaultdict(Counter)
        self.alignments = []

    def record_alignment_details(self, aln, bundle, assignments):
        """Record detailed alignment info if user elects to save diagnostics info with the run

        This is called once per locus per read (every alignment) when the user elects to save
        diagnostics. The recorded information is later written to {library['Name']}_aln_table.txt
        after the entire SAM file has been processed."""

        # Perform reverse complement for anti-sense reads
        read = aln['Seq'] \
            if aln['Strand'] is True \
            else aln['Seq'][::-1].translate(self.complement)

        # sequence, cor_counts, strand, start, end, feat1;feat2;feat3
        self.alignments.append((read, bundle['corr_count'], aln['Strand'], aln['Start'], aln['End'],
                                ';'.join(assignments)))

    def record_diagnostics(self, assignments, n_candidates, aln, bundle):
        """Records basic diagnostic info"""

        if len(assignments) == 0:
            if aln['Strand'] is True:
                self.alignment_diags['Uncounted alignments (+)'] += 1
            else:
                self.alignment_diags['Uncounted alignments (-)'] += 1
            if n_candidates == 0:
                self.alignment_diags['No feature counts'] += bundle['corr_count']
            else:
                self.alignment_diags['Eliminated counts'] += bundle['corr_count']


class MergedDiags(MergedStat):
    def __init__(self):
        self.aln_diags = pd.DataFrame(columns=Diagnostics.aln_diag_categories)
        self.selection_diags = {}
        self.alignment_tables = {}

    def add_library(self, other: LibraryStats):
        other_lib = other.library['Name']
        self.alignment_tables[other_lib] = other.diags.alignments
        self.selection_diags[other_lib] = other.diags.selection_diags
        self.aln_diags.loc[other_lib] = other.diags.alignment_diags

    def write_output_logfile(self):
        self.write_alignment_diags()
        # self.write_selection_diags()  Not currently collected
        self.write_alignment_tables()

    def write_alignment_diags(self):
        MergedStat.df_to_csv(self.aln_diags, "Sample", self.prefix, "alignment_diags", sort_axis="index")

    def write_alignment_tables(self):
        for library_name, table in self.alignment_tables.items():
            outfile = make_filename([self.prefix, library_name, 'aln_table'], ext='.txt')
            with open(outfile, 'w') as imf:
                imf.writelines(
                    # sequence, cor_counts, strand, start, end, feat1a/feat1b;feat2;...
                    map(lambda rec: "%s\t%f\t%c\t%d\t%d\t%s\n" % rec, table)
                )

    def write_selection_diags(self):
        out = []
        for lib in sorted(self.selection_diags.keys()):
            out.append(lib)
            for feat_class in sorted(self.selection_diags[lib].keys()):
                out.append('\t' + ','.join(feat_class))
                for stat in sorted(self.selection_diags[lib][feat_class].keys()):
                    out.append("\t\t%s: %d" % (stat, self.selection_diags[lib][feat_class][stat]))

        selection_summary_filename = make_filename([self.prefix, 'selection_diags'], ext='.txt')
        with open(selection_summary_filename, 'w') as f:
            f.write('\n'.join(out))
