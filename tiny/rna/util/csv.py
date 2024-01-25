import csv
import os

from collections import OrderedDict
from typing import TextIO, Type

def get_csv_dialect(file: TextIO) -> Type[csv.Dialect]:
    """Returns the default CSV dialect with delimiter determined by heuristics
    It is necessary to do this before parsing the CSV because the delimiter is
    ultimately determined by locale where last edited, and locales that use a
    comma for decimal separator might instead use a semicolon delimiter."""

    file_sample = file.read(16 * 1024)  # 16 kb
    heuristic = csv.Sniffer().sniff(file_sample)

    # All other dialect attributes should be default values to ensure compatibility
    # with our own outputs. Currently, the quote heuristics produce incorrect results
    # during testing (e.g. single row Features Sheet with embedded quotes in one field)

    dialect = csv.excel  # the default dialect
    dialect.delimiter = heuristic.delimiter
    file.seek(0)
    return dialect


class CSVReader(csv.DictReader):
    """A simple wrapper class for csv.DictReader

    This makes field labels consistent across the project, simplifies the code, and
    allows for validation and reordering of column names. We also keep track of the
    row number for diagnostic outputs; the base class offers the line_num attribute,
    but line_num != row_num if a record spans multiple lines in the csv.
    """

    # user-facing name -> internal short name
    tinyrna_sheet_fields = {
        "Features Sheet": OrderedDict({
           "Select for...":     "Key",
           "with value...":     "Value",
           "Classify as...":    "Class",
           "Source Filter":     "Filter_s",
           "Type Filter":       "Filter_t",
           "Hierarchy":         "Hierarchy",
           "Strand":            "Strand",
           "5' End Nucleotide": "nt5end",
           "Length":            "Length",
           "Overlap":           "Overlap",
           "Mismatches":        "Mismatch"
        }),
        "Samples Sheet": OrderedDict({
            "Input Files":       "File",
            "Sample/Group Name": "Group",
            "Replicate Number":  "Replicate",
            "Control":           "Control",
            "Normalization":     "Normalization"
        })
    }

    def __init__(self, filename: str, doctype: str = None):
        self.doctype = doctype
        self.tinyrna_file = filename
        self.row_num = 0
        try:
            self.tinyrna_fields = tuple(CSVReader.tinyrna_sheet_fields[doctype].values())
        except KeyError:
            raise ValueError(f'No existing configuration for doctype "{doctype}"')

    def rows(self):
        self.replace_excel_ellipses()
        with open(os.path.expanduser(self.tinyrna_file), 'r', encoding='utf-8-sig', newline='') as f:
            super().__init__(f, fieldnames=self.tinyrna_fields, dialect=get_csv_dialect(f))
            header = next(self)

            # Compatibility check. Column headers are still often changed at this stage
            # and it doesn't make sense to infer column identity
            self.validate_csv_header(header)

            def is_empty_row(row):
                return not any(c.strip() for c in row.values() if c is not None)

            for row in self:
                self.row_num += 1
                if is_empty_row(row): continue
                if self.restkey in row.keys() or self.restval in row.values():
                    raise csv.Error(f"Inconsistent column count (row {self.row_num} of {self.doctype})")
                yield row

        if self.row_num == 0:
            raise csv.Error(f"{self.doctype} is empty or contains only whitespace")

    def replace_excel_ellipses(self):
        """Excel has an autocorrect setting that converts three periods to an ellipsis character.
        The resulting ellipsis is not UTF-8 compatible and causes the decoder to fail. Switching
        to the correct encoding is still problematic because the character, while visually correct,
        doesn't compare equal to "..." so header validation fails. Just replace it."""

        try:
            with open(self.tinyrna_file, 'r', encoding='utf-8') as f:
                next(f)
        except UnicodeDecodeError as e:
            if e.reason == "invalid continuation byte" and e.object[e.start] == 0xc9:
                with open(self.tinyrna_file, 'r', encoding='latin-1') as f:
                    lines = f.readlines()
                    lines[0] = lines[0].replace("\xc9", '...')
                    lines[0].encode('latin-1').decode('utf-8-sig')  # sanity check
                # Encode with utf-8 rather than utf-8-sig for Excel compatibility
                with open(self.tinyrna_file, 'w', encoding='utf-8') as f:
                    f.writelines(lines)
                print(f"Replaced invalid character in {self.doctype} header.")
            else:
                raise

    def validate_csv_header(self, header: OrderedDict):
        # The expected header values
        doc_reference = self.tinyrna_sheet_fields[self.doctype]
        expected = {key.lower() for key in doc_reference.keys()}

        # The header values that were read
        read_vals = {val.lower() for key, val in header.items() if val is not None}
        self.check_backward_compatibility(read_vals)

        # Find differences between actual and expected headers
        unknown = read_vals - expected
        missing = expected - read_vals

        if len(missing):
            raise csv.Error('\n\t'.join([f"The following columns are missing from your {self.doctype}:", *missing]))
        if len(unknown):
            raise csv.Error('\n\t'.join([f"The following columns in your {self.doctype} are unrecognized:", *unknown]))

        doc_ref_lowercase = {key.lower(): value for key, value in doc_reference.items()}
        header_lowercase = {key: value.lower() for key, value in header.items()}

        if tuple(header_lowercase.values()) != tuple(doc_ref_lowercase.keys()):
            # Remap column order to match client's
            self.fieldnames = tuple(doc_ref_lowercase[key] for key in header_lowercase.values())

    def check_backward_compatibility(self, header_vals):
        """Catch column differences from old versions so a helpful error can be provided"""

        header_vals_lc = {val.lower() for val in header_vals}

        compat_errors = []
        if self.doctype == "Features Sheet":
            if len(header_vals_lc & {'alias by...', 'feature source'}):
                compat_errors.append('\n'.join([
                    "It looks like you're using a Features Sheet from an earlier version of",
                    "tinyRNA. Feature aliases and GFF files are now defined in the Paths File.",
                    "Please review the Paths File documentation in Configuration.md, update your",
                    'Paths File, and remove the "Alias by..." and "Feature Source" columns from',
                    "your Features Sheet to avoid this error."
                ]))

            if len(header_vals_lc & {'source filter', 'type filter'}) != 2:
                compat_errors.append('\n'.join([
                    "It looks like you're using a Features Sheet from an earlier version of",
                    "tinyRNA. Source and type filters are now defined in the Features Sheet.",
                    "They are no longer defined in the Run Config. Please review the Stage 1",
                    "section in tiny-count's documentation, then add the new columns",
                    '"Source Filter" and "Type Filter" to your Features Sheet to avoid this error.'
                ]))

            if 'tag' in header_vals_lc:
                compat_errors.append('\n'.join([
                    "It looks like you're using a Features Sheet from a version of tinyRNA",
                    'that offered "tagged counting". The "Tag" header has been repurposed as a feature',
                    "classifier and its meaning within the pipeline has changed. Additionally, feature",
                    "class is no longer determined by the Class= attribute. Please review the Stage 1",
                    'section in tiny-count\'s documentation, then rename the "Tag" column to',
                    '"Classify as..." to avoid this error.'
                ]))

            if 'mismatches' not in header_vals_lc:
                compat_errors.append('\n'.join([
                    "It looks like you're using a Features Sheet from an earlier version of",
                    'tinyRNA. An additional column, "Mismatches", is now expected. Please review',
                    "the Stage 2 section in tiny-count's documentation for more info, then add",
                    "the new column to your Features Sheet to avoid this error."
                ]))

        if self.doctype == "Samples Sheet":
            if 'fastq/sam files' in header_vals_lc:
                compat_errors.append('\n'.join([
                    "It looks like you're using a Samples Sheet from an earlier version of",
                    'tinyRNA. The "FASTQ/SAM files" column has been renamed to "Input Files"',
                    'due to the addition of BAM file support. Please rename the column in',
                    "your Samples Sheet to avoid this error."
                ]))

        if compat_errors: raise ValueError('\n\n'.join(compat_errors))
