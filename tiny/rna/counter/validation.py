import functools
import subprocess
import pysam
import sys
import os

from collections import Counter, defaultdict
from typing import List, Dict

from tiny.rna.counter.hts_parsing import parse_gff, ReferenceFeatures, SAM_reader
from tiny.rna.counter.features import FeatureSelector
from tiny.rna.util import ReportFormatter, sorted_natural, gzip_open


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
        self.ReferenceFeatures = ReferenceFeatures(gff_files)
        self.column_filters = self.build_column_filters(rules)
        self.report = ReportFormatter(self.targets)
        self.chrom_set = set()

        self.seq_files = [ebwt, genomes, alignments]
        self.gff_files = gff_files

    def validate(self):
        print("Validating annotation files... ", end='')
        self.parse_and_validate_gffs(self.gff_files)
        self.validate_chroms(*self.seq_files)
        print("done.")

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
        if not self.ReferenceFeatures.column_filter_match(row, self.column_filters): return

        if row.iv.strand not in ('+', '-'):
            issues['strand'][file] += 1

        try:
            self.ReferenceFeatures.get_feature_id(row)
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


class SamSqValidator:
    """Validates @SQ headers for tiny-count's sequence-based counting mode"""

    targets = {
        "inter sq": "Sequence identifiers with inconsistent lengths",
        "intra sq": "SAM files with repeated sequence identifiers",
        "incomplete sq": "SAM files with incomplete @SQ headers",
        "missing sq": "SAM files that lack @SQ headers"
    }

    def __init__(self, sam_files):
        self.report = ReportFormatter(self.targets)
        self.sam_files = sam_files
        self.reference_seqs = {}
        self.sq_headers = {}

    def validate(self):
        print("Validating sequence identifiers in SAM files... ", end='')
        self.read_sq_headers()
        self.validate_sq_headers()
        print("done.")

        self.report.print_report()
        if self.report.errors:
            sys.exit(1)

    def validate_sq_headers(self):
        """Performs validation tests in the required order
        It is important to check for and return upon syntax error, otherwise
        the remaining tests are likely to raise an exception before the report
        has a chance to be printed"""

        # First verify @SQ header syntax
        missing = self.get_missing_headers()
        incomplete = self.get_incomplete_sq_headers()
        if missing or incomplete:
            self.generate_header_syntax_report(missing, incomplete)
            return

        # Next verify identifier definitions
        duplicate = self.get_duplicate_identifiers()
        ambiguous = self.get_ambiguous_lengths()
        if duplicate or ambiguous:
            self.generate_identifier_report(duplicate, ambiguous)

    def get_ambiguous_lengths(self) -> List[str]:
        """Returns a list of sequence IDs that have inconsistent length definitions.
        Sequences with consistent length definitions are added to self.reference_seqs"""

        seq_lengths = defaultdict(set)
        for sam in self.sam_files:
            for sq in self.sq_headers[sam]:
                seq_id = sq['SN']
                seq_len = int(sq['LN'])
                seq_lengths[seq_id].add(seq_len)

        bad_seqs = []
        for seq_id, lengths in seq_lengths.items():
            if len(lengths) == 1:
                lengths = lengths.pop()
                self.reference_seqs[seq_id] = lengths
            else:
                bad_seqs.append(seq_id)

        return bad_seqs

    def get_duplicate_identifiers(self) -> Dict[str, List[str]]:
        """Returns a dictionary of SAM files that contain duplicate sequence identifiers"""

        bad_files = {}
        for file in self.sam_files:
            sq = self.sq_headers[file]
            id_count = Counter(seq['SN'] for seq in sq)
            duplicates = [seq_id for seq_id, count in id_count.items() if count > 1]
            if duplicates: bad_files[file] = duplicates

        return bad_files

    def get_incomplete_sq_headers(self) -> List[str]:
        """Returns a list of SAM files that have incomplete @SQ headers"""

        return [file for file, sqs in self.sq_headers.items()
                if not all("SN" in sq and "LN" in sq for sq in sqs)]

    def get_missing_headers(self) -> List[str]:
        """Returns a list of SAM files that lack @SQ headers"""

        return [file for file, sqs in self.sq_headers.items()
                if len(sqs) == 0]

    def generate_header_syntax_report(self, missing, incomplete):
        report = {}
        if missing:
            report['missing sq'] = sorted_natural(missing)
        if incomplete:
            report['incomplete sq'] = sorted_natural(incomplete)

        header = "Every SAM file must have complete @SQ headers with SN and LN\n" \
                 "fields when performing sequence-based counting.\n"
        self.report.add_error_section(header, report)

    def generate_identifier_report(self, duplicate, ambiguous):
        report = {}
        if duplicate:
            report['intra sq'] = sorted_natural(duplicate)
        if ambiguous:
            report['inter sq'] = sorted_natural(ambiguous)

        header = "Sequence identifiers must be unique and have consistent length definitions.\n"
        self.report.add_error_section(header, report)

    def read_sq_headers(self):
        for file in self.sam_files:
            sr = SAM_reader()
            pr = pysam.AlignmentFile(file, check_sq=True)
            sr._gather_metadata(pr)

            self.sq_headers[file] = sr._header_dict.get('SQ', [])