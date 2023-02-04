import functools
import subprocess
import sys
import os

from collections import Counter, defaultdict
from typing import List, Dict

from tiny.rna.counter.hts_parsing import parse_gff, ReferenceTables
from tiny.rna.counter.features import FeatureSelector
from tiny.rna.util import sorted_natural, gzip_open


class ReportFormatter:
    error_header   = "========[ ERROR ]=============================================="
    warning_header = "========[ WARNING ]============================================"

    def __init__(self, key_mapper):
        self.key_mapper = key_mapper
        self.warnings = []
        self.errors = []

    def print_report(self):
        for error_section in self.errors:
            print(self.error_header)
            print(error_section)
            print()

        for warning_section in self.warnings:
            print(self.warning_header)
            print(warning_section)
            print()

    def add_error_section(self, header, body=None):
        self.add_section(header, body, self.errors)

    def add_warning_section(self, header, body=None):
        self.add_section(header, body, self.warnings)

    def add_section(self, header, body, dest):
        if isinstance(body, dict):
            body = self.nested_dict(body)
        elif isinstance(body, list):
            body = '\n'.join(body)
        if body:
            dest.append('\n'.join([header, body]))
        else:
            dest.append(header)

    def nested_dict(self, summary, indent='\t'):
        report_lines = self.recursive_indent(summary, indent)
        return '\n'.join(report_lines)

    def recursive_indent(self, mapping: dict, indent: str):
        """Converts a nested dictionary into a properly indented list of lines

        Args:
            mapping: the nested dictionary to be formatted
            indent: the indent string to prepend to lines on the current level
        """

        lines = []
        for key, val in mapping.items():
            if not val: continue
            key_header = f"{indent}{self.key_mapper.get(key, key)}: "
            if isinstance(val, dict):
                lines.append(key_header)
                lines.extend(self.recursive_indent(val, indent + '\t'))
            elif isinstance(val, (list, set)):
                lines.append(key_header)
                lines.extend([indent + '\t' + line for line in map(str, val)])
            else:
                lines.append(key_header + str(val))
        return lines


class GFFValidator:
    """Validates GFF files based on their contents and the contents of sequencing files to which
    the GFF files are expected to be applied."""

    targets = {
        "ID attribute": "Features missing a valid identifier attribute",
        "sam files": "Potentially incompatible SAM alignment files",
        "seq chromosomes": "Chromosomes present in sequence files",
        "gff chromosomes": "Chromosomes present in GFF files",
        "strand": "Features missing strand information",
    }

    def __init__(self, gff_files, rules, ebwt=None, genomes=None, alignments=None):
        self.ReferenceTables = ReferenceTables(gff_files)
        self.column_filters = self.build_column_filters(rules)
        self.report = ReportFormatter(self.targets)
        self.chrom_set = set()

        self.seq_files = [ebwt, genomes, alignments]
        self.gff_files = gff_files

    def validate(self):
        print("Validating annotation files...")
        self.parse_and_validate_gffs(self.gff_files)
        self.validate_chroms(*self.seq_files)
        self.report.print_report()
        if self.report.errors:
            sys.exit(1)

    def parse_and_validate_gffs(self, file_set):
        gff_infractions = defaultdict(Counter)
        for file, *_ in file_set.items():
            row_fn = functools.partial(self.validate_gff_row, issues=gff_infractions, file=file)
            parse_gff(file, row_fn=row_fn)

        if gff_infractions:
            self.generate_gff_report(gff_infractions)

    def validate_gff_row(self, row, issues, file):
        # Skip definitions of whole chromosomes and obey source/type filters
        if row.type.lower() == "chromosome": return
        if not self.ReferenceTables.column_filter_match(row, self.column_filters): return

        if row.iv.strand not in ('+', '-'):
            issues['strand'][file] += 1

        try:
            self.ReferenceTables.get_feature_id(row)
        except:
            issues['ID attribute'][file] += 1

        self.chrom_set.add(row.iv.chrom)

    def generate_gff_report(self, infractions):
        header = "The following issues were found in the GFF files provided. "

        if "strand" in infractions:
            strand_issues = {"strand":
                sorted_natural([
                    f"{count} missing in {file}"
                    for file, count in infractions['strand'].items()
                ], reverse=True)
            }

            ext_header = 'Unstranded features are allowed, but they can lead to potentially unexpected results.\n' \
                         'These features will match "sense", "antisense", and "both" strand selectors. 5\'/3\' anchored\n' \
                         "overlap selectors for these features will evaluate for termini shared with the alignment,\n" \
                         "but will not distinguish between the alignment's 5' and 3' ends."

            self.report.add_warning_section('\n'.join([header, ext_header]), strand_issues)

        if "ID attribute" in infractions:
            idattr_issues = {"ID attribute":
                sorted_natural([
                    f"{count} missing in {file}"
                    for file, count in infractions['ID attribute'].items()
                ], reverse=True)
            }

            self.report.add_error_section(header, idattr_issues)

    def validate_chroms(self, ebwt=None, genomes=None, alignments=None):
        # First search bowtie indexes if they are available
        if ebwt is not None:
            try:
                shared, chroms = self.chroms_shared_with_ebwt(ebwt)
                self.generate_chrom_report(shared, chroms)
                return
            except Exception:
                pass # Fallback to other input options

        # Next search the genome fasta(s) if available
        if genomes is not None:
            shared, chroms = self.chroms_shared_with_genomes(genomes)
            if chroms:
                self.generate_chrom_report(shared, chroms)
                return
            # Continue search if chroms weren't gathered from genomes

        # Preferred inputs aren't available; continue testing with heuristic options
        if alignments is not None:
            self.validate_chroms_heuristic(alignments)
        else:
            self.report.add_warning_section("Shared chromosome identifiers could not be validated.")

    def validate_chroms_heuristic(self, alignments):
        suspect_files = self.alignment_chroms_mismatch_heuristic(alignments)
        self.generate_chrom_heuristics_report(suspect_files)

    def chroms_shared_with_ebwt(self, index_prefix):
        """Determines the complete set of chromosomes in the bowtie index, then finds
        their set intersection with the chromosomes parsed from GFF files.
        Returns:
            a tuple of (shared chromosomes, index chromosomes)"""

        summary = subprocess.check_output(
            ['bowtie-inspect', index_prefix, '-s'],
            stderr=subprocess.DEVNULL
        ).decode('latin1').splitlines()

        ebwt_chroms = set()
        for line in summary:
            if line.startswith("Sequence-"):
                try:
                    ebwt_chroms.add(line.split('\t')[1])
                except IndexError:
                    pass

        shared = ebwt_chroms & self.chrom_set
        return shared, ebwt_chroms

    def chroms_shared_with_genomes(self, genome_fastas):
        """Determines the complete set of chromosomes defined in each genome file,
        then finds their set intersection with the chromosomes parsed from GFF files.
        Returns:
            a tuple of (shared chromosomes, genome chromosomes)"""

        genome_chroms = set()
        for fasta in genome_fastas:
            if not os.path.isfile(fasta):
                continue
            elif fasta.endswith('.gz'):
                file_if = gzip_open
            else:
                file_if = open

            with file_if(fasta, 'rb') as f:
                for line in f:
                    if line[0] == ord('>'):
                        genome_chroms.add(line[1:].strip().decode())

        shared = genome_chroms & self.chrom_set
        return shared, genome_chroms

    def alignment_chroms_mismatch_heuristic(self, sam_files: List[str], subset_size=50000) -> Dict[str, set]:
        """Since alignment files can be very large, we only check that there's at least one shared
        chromosome identifier and only the first subset_size lines are read from each file. The
        returned dictionary contains only the SAM files whose sampled chromosome set failed to
        intersect with the chromosomes parsed from GFF files.
        Returns:
            a dictionary of {SAM filename: SAM chromosomes sampled}"""

        files_wo_overlap = {}

        for file in sam_files:
            file_chroms = set()
            with open(file, 'rb') as f:
                for i, line in zip(range(subset_size), f):
                    if line[0] == ord('@'): continue
                    file_chroms.add(line.split(b'\t')[2].strip().decode())
                    if i % 10000 == 0 and len(file_chroms & self.chrom_set): break

            if not len(file_chroms & self.chrom_set):
                files_wo_overlap[file] = file_chroms

        return files_wo_overlap

    def generate_chrom_report(self, shared, chroms):
        if shared: return
        header = "GFF files and sequence files don't share any chromosome identifiers."
        summary = {
            "gff chromosomes": sorted(self.chrom_set),
            "seq chromosomes": sorted(chroms)
        }

        self.report.add_error_section(header, summary)

    def generate_chrom_heuristics_report(self, suspect_files):
        if not suspect_files: return

        header = "GFF files and sequence files might not contain the same chromosome identifiers.\n" \
                 "This is determined from a subset of each sequence file, so false positives may be reported."
        chroms = {file: [f"Chromosomes sampled: {', '.join(sorted(chroms))}"]
                  for file, chroms in suspect_files.items()}

        summary = {
            "sam files": chroms,
            "gff chromosomes": sorted(self.chrom_set)
        }

        self.report.add_warning_section(header, summary)

    @staticmethod
    def build_column_filters(rules):
        """Builds a "rules table" containing only GFF column selectors.
        All listed filters are gathered and collapsed into a single rule"""

        selector_defs = {
            selector: ','.join([row[selector] for row in rules if row[selector]])
            for selector in ["Filter_s", "Filter_t"]
        }

        return FeatureSelector.build_selectors([selector_defs])[0]