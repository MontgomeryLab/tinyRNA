import pandas as pd
import traceback
import gzip
import json
import mmap
import csv
import sys
import os
import re

from abc import abstractmethod, ABC
from typing import Tuple, Optional, Union, List, DefaultDict
from collections import Counter, defaultdict
from glob import glob

from ..util import get_csv_dialect, CSVReader, make_filename, report_execution_time


class LibraryStats:

    summary_categories = ['Mapped Reads Basis',
                          'Total Assigned Reads', 'Total Unassigned Reads',
                          'Total Assigned Sequences', 'Total Unassigned Sequences',
                          'Assigned Single-Mapping Reads', 'Assigned Multi-Mapping Reads',
                          'Reads Assigned to Single Feature', 'Sequences Assigned to Single Feature',
                          'Reads Assigned to Multiple Features', 'Sequences Assigned to Multiple Features']

    def __init__(self, features, **prefs):
        self.library = {'Name': 'Unassigned', 'File': 'Unassigned', 'Norm': '1'}
        self.diags = Diagnostics(features) if prefs.get('report_diags') else None
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
        self.library_stats['Mapped Reads Basis'] += read_counts

        return {
            'read_count': read_counts,
            'corr_count': corr_counts,
            'loci_count': loci_counts,
            'assigned_ftags': set(),
            'assigned_reads': 0.0,
            'unassigned_reads': 0.0,
            'nt5_assigned': self.assigned_nt_len[nt5],
            'nt5_mapped': self.mapped_nt_len[nt5],
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
            bundle['unassigned_reads'] += corr_count
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
            self.diags.record_assignments(assignments.keys(), aln, bundle, n_candidates)

    def finalize_bundle(self, bundle: dict) -> None:
        """Called at the conclusion of processing each multiple-alignment bundle"""

        assigned_feat_count = len({feat[0] for feat in bundle['assigned_ftags']})
        unassigned_reads = bundle['unassigned_reads']
        assigned_reads = bundle['assigned_reads']

        self.library_stats['Total Unassigned Reads'] += unassigned_reads
        bundle['nt5_mapped'][bundle['seq_length']] += assigned_reads + unassigned_reads

        if assigned_feat_count == 0:
            self.library_stats['Total Unassigned Sequences'] += 1
        else:
            self.library_stats['Total Assigned Sequences'] += 1
            self.library_stats['Total Assigned Reads'] += assigned_reads
            self.library_stats[bundle['mapping_stat']] += assigned_reads
            bundle['nt5_assigned'][bundle['seq_length']] += assigned_reads

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
    def finalize(self):
        """Called once all libraries have been added to the MergedStatsManager.
        The implementation of this method should perform any preparation steps required for
        stats validation, such as building dataframes or sorting columns so the appropriate
        comparisons can be made, but it should not perform any operations that reduce decimal
        precision, such as rounding or conversion to int. Validation depends on this precision
        for internal consistency checks. Rounding should instead be done in write_output_logfile()"""
        ...

    @abstractmethod
    def write_output_logfile(self):
        """Called once all stats classes have been validated.
        The implementation of this method should round counts to the appropriate decimal
        precision or convert float values to int if necessary, then write output CSV files."""
        ...

    @staticmethod
    def add_warning(msg):
        MergedStatsManager.warnings.append(msg)

    @staticmethod
    def df_to_csv(df: pd.DataFrame, prefix: str, postfix: str = None, sort_axis: Optional[str] = "columns"):
        """Rounds counts, optionally sorts an axis, and writes the dataframe to CSV with the appropriate filename.

        Args:
            df: The dataframe to write to CSV
            prefix: The output filename prefix (required)
            postfix: The optional filename postfix. If not provided, a postfix is created from
                the index name, which will be split on space, joined on "_", and made lowercase.
            sort_axis: The axis to sort. Default: columns. If None, neither axis is sorted.
        """

        if postfix is None:
            postfix = '_'.join(map(str.lower, df.index.name.split(' ')))

        if sort_axis is not None:
            out_df = df.round(decimals=2).sort_index(axis=sort_axis)
        else:
            out_df = df.round(decimals=2)

        file_name = make_filename([prefix, postfix])
        out_df.to_csv(file_name)


class MergedStatsManager:
    chrom_misses = Counter()
    warnings = []

    def __init__(self, Features_obj, features_csv, prefs):
        MergedStat.prefix = prefs.get('out_prefix', '')
        self.verify_stats = prefs.get('verify_stats')
        self.prefs = prefs

        self.merged_stats = [
            FeatureCounts(Features_obj), RuleCounts(features_csv),
            AlignmentStats(), SummaryStats(prefs),
            NtLenMatrices(prefs)
        ]

        if prefs.get('report_diags', False):
            self.merged_stats.append(MergedDiags())

    @report_execution_time('Validating and writing report files')
    def write_report_files(self) -> None:
        try:
            for in_process in self.merged_stats:
                in_process.finalize()

            if self.verify_stats:
                # Validate stats now that all libraries have been processed
                StatisticsValidator(self.merged_stats).validate()
        finally:
            self.print_warnings()

        for output_stat in self.merged_stats:
            output_stat.write_output_logfile()

    def add_library_stats(self, other: LibraryStats) -> None:
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
            for chrom in sorted(self.chrom_misses.keys()):
                self.warnings.append("\t" + chrom + ": " + str(self.chrom_misses[chrom]))
            self.warnings.append("\n")

        for warning in self.warnings:
            print(warning, file=sys.stderr)


class FeatureCounts(MergedStat):
    def __init__(self, Features_obj):
        self.feat_counts_df = pd.DataFrame(index=pd.MultiIndex.from_tuples(set.union(*Features_obj.classes.values())))
        self.feat_counts_df.index.names = ["Feature ID", "Classifier"]
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

        # Order columns alphabetically by header (sample_rep_n)
        summary = self.feat_counts_df.sort_index(axis="columns")
        # Add Feature Name column, which is the feature alias (default is blank if no alias exists)
        summary.insert(0, "Feature Name", summary.index.map(lambda feat: ', '.join(self.aliases.get(feat[0], ''))))

        self.feat_counts_df = summary
        self.finalized = True

    def write_feature_counts(self):
        """Writes selected features and their associated counts to {prefix}_feature_counts.csv
        Should be called after finalize()"""

        assert self.finalized, "FeatureCounts object must be finalized before writing output."
        self.df_to_csv(self.feat_counts_df, self.prefix, "feature_counts", sort_axis="index")

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

        self.df_to_csv(counts, self.prefix, "norm_counts", sort_axis="index")


class RuleCounts(MergedStat):
    def __init__(self, features_csv):
        self.rules = self.read_features_sheet(features_csv)
        self.rule_counts_df = pd.DataFrame(index=range(len(self.rules))).rename_axis("Rule Index")
        self.finalized = False

    def add_library(self, other: LibraryStats):
        assert not self.finalized, "Cannot add libraries after RuleCounts object has been finalized."
        self.rule_counts_df[other.library["Name"]] = self.rule_counts_df.index.map(other.rule_counts)

    def finalize(self):
        self.rule_counts_df = self.rule_counts_df.sort_index(axis="columns").fillna(0)

        # Add Rule String column
        rules = self.get_rule_strings()
        self.rule_counts_df.insert(0, "Rule String", rules)

        # Add Mapped Reads row below the counts
        self.rule_counts_df.loc["Mapped Reads"] = SummaryStats.pipeline_stats_df.loc["Mapped Reads"]
        self.finalized = True

    def write_output_logfile(self):
        assert self.finalized, "RuleCounts object must be finalized before writing output."
        self.df_to_csv(self.rule_counts_df, self.prefix, 'counts_by_rule', sort_axis=None)

    def read_features_sheet(self, features_csv):
        """Reads the Features Sheet and returns a list of rules in
        the same order as the FeatureSelector.rules_table. We don't use the
        CSVReader class here because it returns rows with shortened column names,
        and we want to preserve the full column names for the Rule String column
        of the output table."""

        with open(features_csv, 'r', encoding='utf-8-sig') as f:
            return [row for row in csv.DictReader(f, dialect=get_csv_dialect(f))]

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
    def __init__(self, prefs):
        self.nt_len_matrices = {}
        self.norm_gh = prefs.get('normalize_by_genomic_hits', True)
        self.norm_fh = prefs.get('normalize_by_feature_hits', True)
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
            mapped_nt_len_df = self._finalize_mapped(mapped)
            assigned_nt_len_df = self._finalize_assigned(assigned)
            self.nt_len_matrices[lib_name] = (mapped_nt_len_df, assigned_nt_len_df)

        self.finalized = True

    def _finalize_mapped(self, mapped: DefaultDict[str, Counter]) -> pd.DataFrame:
        return (pd.DataFrame(mapped, dtype='float64')
                .rename_axis("Length")
                .sort_index()
                .fillna(0))

    def _finalize_assigned(self, assigned: DefaultDict[str, Counter]) -> pd.DataFrame:
        # Fractional counts due to (loci count) and/or (assigned feature count) normalization
        assigned_nt_len_df = (pd.DataFrame(assigned, dtype='float64')
                              .rename_axis("Length")
                              .sort_index())

        # Drop non-nucleotide columns if they don't contain counts
        assigned_nt_len_df.drop([
            col for col, values in assigned_nt_len_df.items()
            if col not in ['A', 'T', 'G', 'C'] and values.isna().all()
        ], axis='columns', inplace=True)

        # Assign default count of 0 for remaining cases
        return assigned_nt_len_df.fillna(0)

    def write_output_logfile(self):
        """Writes each library's 5' end nucleotide / length matrices to separate files."""

        assert self.finalized, "NtLenMatrices object must be finalized before writing output."

        for lib_name, (mapped, assigned) in self.nt_len_matrices.items():
            sanitized_lib_name = lib_name.replace('/', '_')
            self._write_mapped(mapped, sanitized_lib_name)
            self._write_assigned(assigned, sanitized_lib_name)

    def _write_mapped(self, mapped: pd.DataFrame, lib_name: str):
        """Writes the mapped matrix to a file.

        Conceptually, mapped reads should always be a whole number. If either
        normalization step is disabled, these counts will be fractional. However,
        even when both normalization steps are enabled, we still get fractional
        counts that are very close to a whole number due to the way these counts
        are tallied. Floating point error is at the root of this issue. Up to
        this point, they are fractional and rounded only by the limitations of
        float representation for the purpose of validation."""

        postfix = f"{lib_name}_mapped_nt_len_dist"

        if self.norm_gh ^ self.norm_fh:
            self.df_to_csv(mapped, self.prefix, postfix, sort_axis=None)
        else:
            mapped_whole_numbers = mapped.round(decimals=2).astype('int64')
            self.df_to_csv(mapped_whole_numbers, self.prefix, postfix, sort_axis=None)

    def _write_assigned(self, assigned: pd.DataFrame, lib_name: str):
        """Writes the assigned matrix to a file (no further work required)."""

        self.df_to_csv(assigned, self.prefix, f"{lib_name}_assigned_nt_len_dist", sort_axis=None)


class AlignmentStats(MergedStat):
    def __init__(self):
        self.alignment_stats_df = (pd.DataFrame(index=LibraryStats.summary_categories)
                                   .rename_axis("Alignment Statistics"))
        self.finalized = False

    def add_library(self, other: LibraryStats):
        assert not self.finalized, "Cannot add libraries after AlignmentStats object has been finalized."
        library, stats = other.library['Name'], other.library_stats
        self.alignment_stats_df[library] = self.alignment_stats_df.index.map(stats)

    def finalize(self):
        self.alignment_stats_df.drop(index="Mapped Reads Basis", inplace=True)  # Handled by SummaryStats
        self.alignment_stats_df.sort_index(axis="columns", inplace=True)
        self.finalized = True

    def write_output_logfile(self):
        assert self.finalized, "AlignmentStats object must be finalized before writing output."
        self.df_to_csv(self.alignment_stats_df, self.prefix, 'alignment_stats')


class SummaryStats(MergedStat):

    conditional_categories  = ["Total Reads", "Retained Reads", "Unique Sequences"]
    counted_categories      = ["Mapped Sequences", "Normalized Mapped Reads", "Mapped Reads", "Assigned Reads"]
    summary_categories      = conditional_categories + counted_categories

    pipeline_stats_df = pd.DataFrame(index=summary_categories).rename_axis("Summary Statistics")

    def __init__(self, prefs):
        self.norm_gh = prefs.get('normalize_by_genomic_hits', True)
        self.norm_fh = prefs.get('normalize_by_feature_hits', True)
        self.missing_collapser_outputs = []
        self.missing_fastp_outputs = []
        self.finalized = False

    def add_library(self, other: LibraryStats):
        """Add incoming summary stats as new column in the master table"""

        assert not self.finalized, "Cannot add libraries after SummaryStats object has been finalized."

        other.library['basename'] = self.get_library_basename(other)
        other_summary = {
            "Mapped Sequences": self.calculate_mapped_seqs(other),
            "Mapped Reads": self.calculate_mapped_reads(other),
            "Assigned Reads": other.library_stats["Total Assigned Reads"]
        }

        if not (self.norm_gh and self.norm_fh):
            norm_mapped = other.library_stats['Mapped Reads Basis']
            other_summary['Normalized Mapped Reads'] = norm_mapped

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
    def calculate_mapped_seqs(other: LibraryStats):
        return sum([other.library_stats[stat] for stat in [
            "Total Assigned Sequences",
            "Total Unassigned Sequences"
        ]])

    @staticmethod
    def calculate_mapped_reads(other: LibraryStats):
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
        self.pipeline_stats_df.sort_index(axis="columns", inplace=True)
        self.finalized = True

    def write_output_logfile(self):
        assert self.finalized, "SummaryStats object must be finalized before writing output."

        if self.norm_gh and self.norm_fh:
            # No need to differentiate when default normalization is on
            irrelevant = "Normalized Mapped Reads"
            out_df = self.pipeline_stats_df.drop(index=irrelevant)
        else:
            # With normalization steps turned off, the Mapped Reads stat might be
            # inflated but necessary for calculating proportions in tiny-plot.
            differentiate = {'Mapped Reads': "Non-normalized Mapped Reads"}
            out_df = self.pipeline_stats_df.rename(index=differentiate)

        self.df_to_csv(out_df, self.prefix, "summary_stats")

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

    summary_categories = ['Eliminated counts', 'No feature counts',
                          'Uncounted alignments (+)', 'Uncounted alignments (-)']

    alignment_columns =  ["Sequence", "Raw Count", "Normalized Count", "Genomic Hits",
                          "Chrom", "Strand", "Start", "End", "Mismatches",
                          "Candidates", "Feature Hits", "Feature Aliases"]

    complement = bytes.maketrans(b'ACGTacgt', b'TGCAtgca')
    map_strand = {True: '+', False: '-', None: '.'}

    def __init__(self, Features_obj):
        self.assignment_diags = {stat: 0 for stat in Diagnostics.summary_categories}
        self.selection_diags = defaultdict(Counter)
        self.aliases = Features_obj.aliases
        self.alignments = []

    def record_assignments(self, assignments, alignment, bundle, n_candidates):
        self.record_alignment_details(assignments, alignment, bundle, n_candidates)
        self.record_summary_diagnostics(assignments, alignment, bundle, n_candidates)

    def record_alignment_details(self, assignments, aln, bundle, n_candidates):
        """Record detailed alignment info if user elects to save diagnostics info with the run

        This is called once per locus per read (every alignment) when the user elects to save
        diagnostics. The recorded information is later written to {library['Name']}_aln_table.txt
        after the entire alignment file has been processed."""

        # Map internal strand representation to +/-/.
        strand = self.map_strand[aln['Strand']]

        # Perform reverse complement for anti-sense reads
        seq = aln['Seq'] if strand == '+' \
            else aln['Seq'][::-1].translate(self.complement)

        # For easy parsing, report as: (id, classifier); ...
        feature_hits = '; '.join(f"({fid}, {tag})" for fid, tag in assignments)

        # For easy parsing, report as: (alias1, alias2, ...); ...
        feat_aliases = [', '.join(self.aliases.get(fid, '')) for fid, _ in assignments]
        feat_aliases = '; '.join(f"({aliases})" for aliases in feat_aliases)

        counts = (bundle['read_count'], bundle['corr_count'], bundle['loci_count'])
        pos = (aln['Chrom'], strand, aln['Start'] + 1, aln['End'], aln['NM'])
        row = (seq, *counts, *pos, n_candidates, feature_hits, feat_aliases)

        self.alignments.append(row)

    def record_summary_diagnostics(self, assignments, aln, bundle, n_candidates):
        """Records basic diagnostic info"""

        if len(assignments) == 0:
            if aln['Strand'] is True:
                self.assignment_diags['Uncounted alignments (+)'] += 1
            else:
                self.assignment_diags['Uncounted alignments (-)'] += 1
            if n_candidates == 0:
                self.assignment_diags['No feature counts'] += bundle['corr_count']
            else:
                self.assignment_diags['Eliminated counts'] += bundle['corr_count']


class MergedDiags(MergedStat):
    def __init__(self):
        self.assignment_diags = pd.DataFrame(index=Diagnostics.summary_categories).rename_axis("Assignment Diags")
        self.selection_diags = {}
        self.alignment_tables = {}

    def add_library(self, other: LibraryStats):
        other_lib = other.library['Name']
        self.alignment_tables[other_lib] = other.diags.alignments
        # self.selection_diags[other_lib] = other.diags.selection_diags  Not currently collected
        self.assignment_diags[other_lib] = self.assignment_diags.index.map(other.diags.assignment_diags)

    def finalize(self): pass  # Nothing to do

    def write_output_logfile(self):
        self.write_assignment_diags()
        # self.write_selection_diags()  Not currently collected
        self.write_alignment_tables()

    def write_assignment_diags(self):
        self.df_to_csv(self.assignment_diags, self.prefix)

    def write_alignment_tables(self):
        header = Diagnostics.alignment_columns
        for library_name, table in self.alignment_tables.items():
            outfile = make_filename([self.prefix, library_name, 'alignment_table'], ext='.csv')
            with open(outfile, 'w', newline='') as ao:
                csv_writer = csv.writer(ao)
                csv_writer.writerow(header)
                csv_writer.writerows(table)

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


class StatisticsValidator:

    def __init__(self, merged_stats: List[MergedStat]):
        self.stats = {type(stat).__name__: stat for stat in merged_stats}
        self.prefix = MergedStat.prefix

    def validate(self):
        """Check the counts reported in all MergedStat objects for internal consistency.
        This step should be allowed to fail without preventing outputs from being written."""

        try:
            self.validate_stats(self.stats['AlignmentStats'], self.stats['SummaryStats'])
            self.validate_counts(self.stats['FeatureCounts'], self.stats['RuleCounts'], self.stats['SummaryStats'])
            self.validate_nt_len_matrices(self.stats['NtLenMatrices'], self.stats['SummaryStats'])
        except Exception as e:
            print(traceback.format_exc(), file=sys.stderr)
            MergedStat.add_warning("Error validating statistics: " + str(e))

    def validate_stats(self, aln_stats, sum_stats):
        # SummaryStats categories
        MS, MR, AR, NMR = "Mapped Sequences", "Mapped Reads", "Assigned Reads", "Normalized Mapped Reads"

        # AlignmentStats categories
        TAR, TUR = "Total Assigned Reads", "Total Unassigned Reads"
        TAS, TUS = "Total Assigned Sequences", "Total Unassigned Sequences"
        ASMR, AMMR = "Assigned Single-Mapping Reads", "Assigned Multi-Mapping Reads"
        RASF, RAMF = "Reads Assigned to Single Feature", "Reads Assigned to Multiple Features"
        SASF, SAMF = "Sequences Assigned to Single Feature", "Sequences Assigned to Multiple Features"

        a_df = aln_stats.alignment_stats_df
        s_df = sum_stats.pipeline_stats_df
        table = pd.concat([s_df.loc[[MR, MS], :], a_df])

        res = pd.DataFrame(columns=a_df.columns)
        res.loc["TAS - (SASF + SAMF)"] = table.loc[TAS] - (table.loc[SASF] + table.loc[SAMF])
        res.loc["TAR - (ASMR + AMMR)"] = table.loc[TAR] - (table.loc[ASMR] + table.loc[AMMR])
        res.loc["TAR - (RASF + RAMF)"] = table.loc[TAR] - (table.loc[RASF] + table.loc[RAMF])
        res.loc["MS - (TUS + TAS)"] = table.loc[MS] - (table.loc[TUS] + table.loc[TAS])
        res.loc["MS - (TUS + SASF + SAMF)"] = table.loc[MS] - (table.loc[TUS] + table.loc[SASF] + table.loc[SAMF])
        res.loc["MR - (TUR + TAR)"] = table.loc[MR] - (table.loc[TUR] + table.loc[TAR])
        res.loc["MR - (TUR + RASF + RAMF)"] = table.loc[MR] - (table.loc[TUR] + table.loc[RASF] + table.loc[RAMF])
        res.loc["MR - (TUR + ASMR + AMMR)"] = table.loc[MR] - (table.loc[TUR] + table.loc[ASMR] + table.loc[AMMR])
        res.loc["Sum"] = res.sum(axis=0)
        res = res.round(decimals=2)

        if not self.approx_equal(a_df.loc[TAR], s_df.loc[AR]):
            MergedStat.add_warning("Alignment stats and summary stats disagree on the number of assigned reads.")
            MergedStat.add_warning(self.indent_table(a_df.loc[TAR], s_df.loc[AR], index=("Alignment Stats", "Summary Stats")))

        if not self.approx_equal(s_df.loc[NMR], s_df.loc[MR]) and s_df.loc[NMR].gt(s_df.loc[MR]).any():
            MergedStat.add_warning("Summary stats reports normalized mapped reads > non-normalized mapped reads.")
            MergedStat.add_warning(self.indent_table(s_df.loc[NMR], s_df.loc[MR], index=("Normalized", "Non-normalized")))

        if not self.approx_equal(res.loc["Sum"].sum(), 0):
            outfile_name = "stats_check"
            MergedStat.add_warning("Unexpected disagreement between alignment and summary statistics.")
            MergedStat.add_warning(f"See the checksum table ({outfile_name}.csv)")
            MergedStat.df_to_csv(res, self.prefix, outfile_name)

    def validate_counts(self, feat_counts, rule_counts, sum_stats):
        f_df = feat_counts.feat_counts_df
        r_df = rule_counts.rule_counts_df
        s_df = sum_stats.pipeline_stats_df
        s_sum = s_df.loc["Assigned Reads"]

        # Get sum of assigned counts in the RuleCounts table
        r_rows = r_df.index.difference(["Mapped Reads"], sort=False)
        r_cols = r_df.columns.difference(["Rule Index", "Rule String"], sort=False)
        r_df_sum = r_df.loc[r_rows, r_cols].sum(axis=0)
        r_df_mapped = r_df.loc["Mapped Reads", r_cols]

        # Get sum of assigned counts in the FeatureCounts table
        f_cols = f_df.columns.difference(["Feature Name"], sort=False)
        f_df_sum = f_df[f_cols].sum(axis=0)

        # Compare assigned counts in FeatureCounts, RuleCounts, and SummaryStats
        if not self.approx_equal(r_df_sum, f_df_sum):
            MergedStat.add_warning("Rule counts and feature counts disagree on the number of assigned reads:")
            MergedStat.add_warning(self.indent_table(r_df_sum, f_df_sum, index=("Rule Counts", "Feature Counts")))
        if not self.approx_equal(r_df_sum, s_sum):
            MergedStat.add_warning("Rule counts and summary stats disagree on the number of assigned reads:")
            MergedStat.add_warning(self.indent_table(r_df_sum, s_sum, index=("Rule Counts", "Summary Stats")))
        if r_df_sum.gt(r_df_mapped).any():
            MergedStat.add_warning("Rule counts reports assigned reads > mapped reads.")
            MergedStat.add_warning(self.indent_table(r_df_sum, r_df_mapped, index=("Assigned", "Mapped")))

        # Compare counts per classifier in FeatureCounts to corresponding rules in RuleCounts
        f_cls_sums = f_df[f_cols].groupby(level="Classifier").sum()
        r_cls_idxs = rule_counts.get_inverted_classifiers()
        for cls, counts in f_cls_sums.iterrows():
            f_cls_sum = f_cls_sums.loc[cls].sum()
            r_cls_sum = sum(r_df.loc[rule, r_cols].sum() for rule in r_cls_idxs[cls])
            if not self.approx_equal(f_cls_sum, r_cls_sum):
                MergedStat.add_warning(f'Feature counts and rule counts disagree on counts for the "{cls}" classifier.')
                MergedStat.add_warning(self.indent_vs(f_cls_sum, r_cls_sum))

    def validate_nt_len_matrices(self, matrices, sum_stats):
        s_df = sum_stats.pipeline_stats_df
        w_mapped = "Length matrix for {} disagrees with summary stats on the number of mapped reads."
        w_assigned = "Length matrix for {} disagrees with summary stats on the number of assigned reads."

        for lib_name, (mapped, assigned) in matrices.nt_len_matrices.items():
            m_sum = mapped.sum().sum()
            a_sum = assigned.sum().sum()
            if not self.approx_equal(m_sum, s_df.loc["Mapped Reads", lib_name]):
                MergedStat.add_warning(w_mapped.format(lib_name))
                MergedStat.add_warning(self.indent_vs(m_sum, s_df.loc['Mapped Reads', lib_name]))
            if not self.approx_equal(a_sum, s_df.loc["Assigned Reads", lib_name]):
                MergedStat.add_warning(w_assigned.format(lib_name))
                MergedStat.add_warning(self.indent_vs(a_sum, s_df.loc['Assigned Reads', lib_name]))

    @staticmethod
    def approx_equal(x: Union[pd.Series, int], y: Union[pd.Series, int], tolerance=1.0):
        """Tolerate small differences in floating point numbers due to floating point error.
        The error tends to be greater for stats that are incremented many times by small
        amounts. The default tolerance, while conceptually appropriate, is excessive."""

        approx = abs(x - y) < tolerance
        if isinstance(approx, pd.Series):
            approx = approx.all()

        return approx

    @staticmethod
    def indent_table(*series: pd.Series, index: tuple, indent=1):
        """Joins the series into a table and adds its string
        representation to the warning log with each line indented."""

        table_str = str(pd.DataFrame(tuple(series), index=index))
        table_ind = '\n'.join('\t' * indent + line for line in table_str.splitlines())
        return table_ind

    @staticmethod
    def indent_vs(stat_a: Union[int, float], stat_b: Union[int, float], indent=1):
        """Simply indents the string representation of two numbers being compared
        and formats them to two decimal places."""

        return '\t' * indent + f"{stat_a:.2f} vs {stat_b:.2f}"
