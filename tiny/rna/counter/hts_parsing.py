import numpy as np
import HTSeq
import sys
import re
import os

from collections import Counter, defaultdict
from typing import Tuple, List, Dict
from tiny.rna.util import report_execution_time

# For parse_GFF_attribute_string()
_re_attr_main = re.compile(r"\s*([^\s=]+)[\s=]+(.*)")
_re_attr_empty = re.compile(r"^\s*$")


class Alignment:
    """The data structure in which parsed SAM alignments are stored.

    Strand-non-specific 5' end nucleotide is stored for efficient lookup. This
    allows us to skip performing full reverse complement of sequences aligned
    to the antisense strand.
    """

    complement = {ord('A'): 'T', ord('T'): 'A', ord('G'): 'C', ord('C'): 'G'}

    class Sequence:
        def __init__(self, name, seq, nt5):
            self.name = name
            self.len = len(seq)
            self.seq = seq
            self.nt5 = nt5

        def __repr__(self):
            return f"<Sequence Object: '{self.name}', {self.seq} ({self.len} bases)>"

        def __len__(self):
            return self.len

    def __init__(self, iv, name, seq):
        nt5 = self.complement[seq[-1]] if iv.strand == '-' else chr(seq[0])
        self.read = self.Sequence(name, seq, nt5)
        self.iv = iv

    def __repr__(self):
        return f"<Alignment Object: Read '{self.read.name}' aligned to {self.iv}>"


def read_SAM(file):
    """A minimal SAM parser which only handles data relevant to the workflow, for performance."""

    with open(file, 'rb') as f:
        line = f.readline()

        # Skip headers
        while line[0] == ord('@'):
            line = f.readline()

        while line:
            cols = line.split(b'\t')

            # Note: we assume sRNA sequencing data is NOT reversely stranded
            strand = "+" if (int(cols[1]) & 16) >> 4 == 0 else "-"
            chrom = cols[2].decode('utf-8')
            name = cols[0].decode('utf-8')
            start = int(cols[3]) - 1
            seq = cols[9]

            iv = HTSeq.GenomicInterval(chrom, start, start + len(seq), strand)
            yield Alignment(iv, name, seq)

            line = f.readline()


def infer_strandedness(sam_file: str, intervals: dict) -> str:
    """Infers strandedness from a sample SAM file and intervals from a parsed GFF file

    Credit: this technique is an adaptation of those in RSeQC's infer_experiment.py.
    It has been modified to accept a GFF reference file rather than a BED file,
    and to use HTSeq rather than bx-python.
    """

    unstranded = HTSeq.GenomicArrayOfSets("auto", stranded=False)

    for orig_iv in intervals.values():
        iv_convert = orig_iv.copy()
        iv_convert.strand = '.'
        unstranded[iv_convert] = orig_iv.strand

    # Assumes read_SAM() defaults to non-reverse strandedness
    sample_read = read_SAM(sam_file)
    gff_sam_map = Counter()
    for count in range(1, 20000):
        try:
            rec = next(sample_read)
            rec.iv.strand = '.'
            gff_strand = ':'.join(unstranded[rec.iv].get_steps())
            sam_strand = rec.iv.strand
            gff_sam_map[sam_strand + gff_strand] += 1
        except StopIteration:
            break

    non_rev = (gff_sam_map['++'] + gff_sam_map['--']) / sum(gff_sam_map.values())
    reverse = (gff_sam_map['+-'] + gff_sam_map['-+']) / sum(gff_sam_map.values())
    unknown = 1 - reverse - non_rev

    if reverse > non_rev: return "reverse"
    else: return "non-reverse"


def parse_GFF_attribute_string(attrStr, extra_return_first_value=False):
    """Parses a GFF attribute string and returns it as a dictionary.

    This is a slight modification of the same method found in HTSeq.features.
    It has been adapted to allow features to have comma separated attribute
    value lists. For downstream compatibility with membership operations
    (e.g. (Attribute Key, Attribute Value) rules when matching feature
    candidates) non-list attribute values are also recorded as tuples.

    Per the original HTSeq docstring:
        "If 'extra_return_first_value' is set, a pair is returned: the dictionary
        and the value of the first attribute. This might be useful if this is the
        ID."
    """

    attribute_dict = {}
    first_val = "_unnamed_"
    for i, attr in enumerate(HTSeq._HTSeq.quotesafe_split(attrStr.rstrip().encode())):
        attr = attr.decode()
        if _re_attr_empty.match(attr):
            continue
        if attr.count('"') not in (0, 2):
            raise ValueError(
                "The attribute string seems to contain mismatched  quotes.")
        mo = _re_attr_main.match(attr)
        if not mo:
            raise ValueError("Failure parsing GFF attribute line")
        key = mo.group(1)
        val = mo.group(2)
        if val.startswith('"') and val.endswith('"'):
            val = val[1:-1]
        # Modification: allow for comma separated attribute values
        attribute_dict[sys.intern(key)] = (sys.intern(val),) \
            if ',' not in val \
            else tuple(c.strip() for c in val.split(','))
        if extra_return_first_value and i == 0:
            first_val = val
    if extra_return_first_value:
        return attribute_dict, first_val
    else:
        return attribute_dict


class ReferenceTables:
    """A GFF parser which builds feature, attribute, and alias tables, with intelligent appends

    Features may be defined by multiple GFF files. If multiple files offer different attributes for
    the same feature, the unique among those attributes are appended to the record. If multiple aliases
    (or Name Attributes, per the Features Sheet) are defined for a feature, the unique among those
    names are appended. Each GFF file defined in the Features Sheet is parsed only once regardless of
    the number of Name Attributes associated with it. Each feature ID may be defined for only one
    interval; if multiple interval definitions are supplied, then these feature IDs are renamed on
    the basis of their source GFF file.
    """

    source_filter = []
    type_filter = []

    def __init__(self, gff_files: Dict[str, list], rules: List[dict], **kwargs):
        self.gff_files = gff_files
        self.feats = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        self.attrs, self.alias, self.intervals = {}, {}, defaultdict(list)

        # Obtain an ordered list of unique attributes of interest from selection rules
        self.attrs_of_interest = list(np.unique(["Class"] + [rule['Identity'][0] for rule in rules]))

        # Patch the GFF attribute parser to support comma separated attribute value lists
        setattr(HTSeq.features, 'parse_GFF_attribute_string', parse_GFF_attribute_string)

        self.set_filters(**kwargs)

    @report_execution_time("GFF parsing")
    def get(self) -> Tuple['HTSeq.GenomicArrayOfSets', Dict[str, list], dict, dict]:
        """Initiates GFF parsing and returns the resulting reference tables"""

        for file, alias_keys in self.gff_files.items():
            gff = HTSeq.GFF_Reader(file)
            for row in gff:
                if row.iv.strand == ".":
                    raise ValueError(f"Feature {row.name} in {file} has no strand information.")
                if not self.filter_match(row):
                    continue

                try:
                    feature_id = row.attr["ID"][0]
                    # Add feature_id <-> feature_interval records
                    root_id = self.add_feature_iv(feature_id, row)
                    # Append alias to feat if unique
                    self.add_alias(alias_keys, root_id, row.attr)
                    # Select attributes of interest from row
                    row_attrs = self.get_interesting_attrs(feature_id, row.attr)
                except KeyError as ke:
                    raise ValueError(f"Feature {row.name} does not contain a {ke} attribute in {file}")

                # Add feature_id -> feature_attributes record
                self.incorporate_attributes(root_id, row_attrs)

        return self.feats, self.attrs, self.alias, self.intervals

    @classmethod
    def filter_match(cls, row):
        select = True
        if len(cls.source_filter):
            select &= row.source in cls.source_filter
        if len(cls.type_filter):
            select &= row.type in cls.type_filter
        return select

    def get_root_feature(self, feature_id: str, feature_attrs: List[Tuple[str, tuple]]) -> str:
        """Returns the ID of the feature's root parent if one exists. Otherwise the original ID is returned."""

        if "Parent" not in feature_attrs:
            return feature_id

        try:
            feature_id = self.get_parent(feature_attrs)
            # Keep checking parent of parents until the root parent is found
            while any(attr for attr in self.attrs[feature_id] if attr[0] == 'Parent'):
                feature_id = self.get_parent(self.attrs[feature_id])

            return feature_id
        except KeyError as ke:
            raise ValueError(f"Feature ID {ke} is referenced as a parent before being defined. Please "
                            "ensure that it occurs before children in your GFF file.")

    def add_feature_iv(self, feature_id: str, row) -> str:
        """Adds the new feature and its intervals to reference tables, then returns its attributes of interest"""

        root_id = self.get_root_feature(feature_id, row.attr)

        self.feats[row.iv] += root_id
        if row.iv not in self.intervals[root_id]:
            self.intervals[root_id].append(row.iv)

        return root_id

    def add_alias(self, alias_keys, root_id, row_attr):
        """"""

        curr_alias = self.alias.get(root_id, ())
        for key in alias_keys:
            # Append to feature's aliases if it does not already contain
            if row_attr[key] not in curr_alias:
                # Add feature_id -> feature_alias_tuple record
                curr_alias += row_attr[key]
        self.alias[root_id] = curr_alias

    def get_interesting_attrs(self, local_id, row_attrs):
        # Copy only the attributes of interest
        interests = [(interest, row_attrs[interest]) for interest in self.attrs_of_interest]
        if "Parent" in row_attrs:
            self.attrs[local_id] = interests + [('Parent', tuple(row_attrs['Parent']))]

        return interests

    def incorporate_attributes(self, root_id, row_attrs):
        if root_id in self.attrs and row_attrs != self.attrs[root_id]:
            # If an attribute record already exists for this feature, and this row provides new attributes,
            #  append the new attribute values to the existing values
            cur_attrs = self.attrs[root_id]
            new_attrs = []
            for cur, new in zip(cur_attrs, row_attrs):
                attr_key = cur[0]
                updated_vals = set(cur[1] + new[1])
                new_attrs.append((attr_key, tuple(updated_vals)))

            self.attrs[root_id] = new_attrs
        else:
            self.attrs[root_id] = row_attrs

    @staticmethod
    def get_parent(feature_attrs):
        parent = []

        # For GFF row attrs
        if type(feature_attrs) is dict:
            parent = feature_attrs["Parent"]
        # For processed attrs (ancestor lookup)
        elif type(feature_attrs) is list:
            for attr in feature_attrs:
                if attr[0] == "Parent":
                    parent = attr[1]

        if len(parent) > 1:
            raise ValueError(f"{feature_attrs['ID']} defines multiple parents which is unsupported at this time.")

        return parent[0]

    @classmethod
    def set_filters(cls, **kwargs):
        for pref, val in kwargs:
            setattr(cls, pref, val)
