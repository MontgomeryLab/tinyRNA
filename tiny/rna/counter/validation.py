import functools
import subprocess
import sys

from collections import Counter, defaultdict

from rna.counter.hts_parsing import parse_gff, ReferenceTables


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

    def recursive_indent(self, mapping, indent):
        lines = []
        for key, val in mapping.items():
            if not val: return lines
            key_header = f"{indent}{self.key_mapper.get(key, key)}: "
            if isinstance(val, dict):
                lines.append(key_header)
                lines.extend(self.recursive_indent(val, indent + '\t'))
            elif isinstance(val, list):
                lines.append(key_header)
                lines.extend([indent + '\t' + line for line in val])
            else:
                lines.append(key_header + str(val))
        return lines

    def indent(self, lines, level=1, sep='\n'):
        ind_token = '\t' * level
        out = ind_token
        out += (sep + ind_token).join(lines)
        return out


class GFFValidator:
    """Validates GFF files based on their contents and the contents of sequencing files to which
    the GFF files are expected to be applied."""

    targets = {
        "ID attribute": "Features missing a valid identifier attribute",
        "seq chromosomes": "Chromosomes present in sequence files",
        "gff chromosomes": "Chromosomes present in GFF files",
        "strand": "Features missing strand information",
    }

    def __init__(self, gff_files, prefs, ebwt=None, genomes=None, alignments=None):
        self.ReferenceTables = ReferenceTables(gff_files, None, **prefs)
        self.report = ReportFormatter(self.targets)
        self.chrom_set = set()

        self.seq_files = [ebwt, genomes, alignments]
        self.gff_files = gff_files
        self.prefs = prefs

    def validate(self):
        self.parse_and_validate_gffs(self.gff_files)
        self.validate_chroms(*self.seq_files)
        self.report.print_report()
        if self.report.errors:
            sys.exit(1)

    def parse_and_validate_gffs(self, file_set):
        gff_infractions = defaultdict(Counter)
        for file, *_ in file_set.items():
            row_fn = functools.partial(self.validate_gff_row, report=gff_infractions[file])
            parse_gff(file, row_fn=row_fn)

        if len(gff_infractions.values()):
            self.generate_gff_report(gff_infractions)

    def validate_gff_row(self, row, report):
        # Check for reasons to normally skip row
        if row.type.lower() == "chromosome": return             # Skip definitions of whole chromosomes regardless
        if not self.ReferenceTables.filter_match(row): return   # Obey source/type filters before validation

        if row.iv.strand not in ('+', '-'):
            report["strand"] += 1

        try:
            self.ReferenceTables.get_feature_id(row)
        except:
            report['ID attribute'] += 1

        self.chrom_set.add(row.iv.chrom)

    def generate_gff_report(self, infractions):
        header = "The following issues were found in the GFF files provided:"
        self.report.add_error_section(header, infractions)

    def validate_chroms(self, ebwt=None, genomes=None, alignments=None):
        # First search bowtie indexes if they are available
        if ebwt is not None:
            try:
                chroms, shared = self.chroms_shared_with_ebwt(ebwt)
                self.generate_chrom_report(chroms, shared)
                return
            except Exception:
                pass # Fallback to other input options

        # Next search the genome fasta(s) if available
        if genomes is not None:
            chroms, shared = self.chroms_shared_with_genomes(genomes)
            self.generate_chrom_report(chroms, shared)
            return

        # Preferred inputs aren't available; continue testing with heuristic options
        if alignments is not None:
            self.validate_chroms_heuristic(alignments)
        else:
            self.report.add_warning_section("Shared chromosome identifiers could not be validated.")

    def validate_chroms_heuristic(self, alignments):
        suspect_files = self.alignment_chroms_mismatch_heuristic(alignments)
        self.generate_chrom_heuristics_report(suspect_files)

    def chroms_shared_with_ebwt(self, index_prefix):
        """Returns the set intersection between parsed GFF chromosomes and those in the bowtie index"""

        summary = subprocess.check_output(['bowtie-inspect', index_prefix, '-s']).decode('latin1').splitlines()
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
        """Returns the set intersection between parsed GFF chromosomes and those in the bowtie index"""

        genome_chroms = set()
        for fasta in genome_fastas:
            with open(fasta, 'rb') as f:
                for line in f:
                    if line[0] == ord('>'):
                        genome_chroms.add(line[1:].strip().decode())

        shared = genome_chroms & self.chrom_set
        return shared, genome_chroms

    def alignment_chroms_mismatch_heuristic(self, sam_files):
        """Since alignment files can be very large, we only check that there's at least one shared
        chromosome identifier and only the first subset_size lines are read from each file."""

        subset_size = 5000
        files_wo_overlap = []

        for file in sam_files:
            file_chroms = set()
            with open(file, 'rb') as f:
                for line, i in zip(f, range(subset_size)):
                    if line[0] == b"@": continue
                    file_chroms.add(line.split(b'\t')[2])
                    if i % 1000 == 0 and len(file_chroms & self.chrom_set): break

            if not len(file_chroms & self.chrom_set):
                files_wo_overlap.append(file)

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
        summary = {
            "The following sequence files might be incompatible": sorted(suspect_files),
            "The following chromosomes are present in GFF files": sorted(self.chrom_set)
        }

        self.report.add_warning_section(header, summary)