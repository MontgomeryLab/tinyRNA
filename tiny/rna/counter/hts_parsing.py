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


@report_execution_time("GFF parsing")
def build_reference_tables(gff_files: Dict[str, list], rules: List[dict]) \
        -> Tuple['HTSeq.GenomicArrayOfSets', Dict[str, list], dict, dict]:
    """A GFF parser which builds feature, attribute, and alias tables, with intelligent appends

    Features may be defined by multiple GFF files. If multiple files offer different attributes for
    the same feature, the unique among those attributes are appended to the record. If multiple aliases
    (or Name Attributes, per the Features Sheet) are defined for a feature, the unique among those
    names are appended. Each GFF file defined in the Features Sheet is parsed only once regardless of
    the number of Name Attributes associated with it. Each feature ID may be defined for only one
    interval; if multiple interval definitions are supplied, then these feature IDs are renamed on
    the basis of their source GFF file.
    """

    # Patch the GFF attribute parser to support comma separated attribute value lists
    setattr(HTSeq.features, 'parse_GFF_attribute_string', parse_GFF_attribute_string)

    # Obtain an ordered list of unique attributes of interest from selection rules
    attrs_of_interest = list(np.unique(["Class"] + [attr['Identity'][0] for attr in rules]))

    feats = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    attrs, alias, intervals = {}, {}, defaultdict(list)

    def check_ancestors(feature_id, row):
        feature_id = row.attr.get("Parent", (feature_id,))[0]

        # Climb generational tree to find root parent
        while any(attr for attr in attrs.get(feature_id, []) if len(attr) and attr[0] == 'Parent'):
            feature_id = attrs[feature_id]['Parent'][0]

        return feature_id

    def add_feature_iv(feature_id, row):
        feats[row.iv] += feature_id
        if row.iv not in intervals[feature_id]:
            intervals[feature_id].append(row.iv)

        # Copy only the attributes of interest
        interests = [(interest, row.attr[interest]) for interest in attrs_of_interest]
        if "Parent" in row.attr:
            interests.append(('Parent', tuple(row.attr['Parent'])))

        return interests

    def add_alias(feature_id, row_attr):
        curr_alias = alias.get(feature_id, ())
        for pref_id in preferred_ids:
            # Append to feature's aliases if it does not already contain
            if row_attr[pref_id] not in curr_alias:
                # Add feature_id -> feature_alias_tuple record
                curr_alias += row_attr[pref_id]
        alias[feature_id] = curr_alias

    def incorporate_attributes(feature_id, row_attrs):
        if feature_id in attrs and row_attrs != attrs[feature_id]:
            # If an attribute record already exists for this feature, and this row provides new attributes,
            #  append the new attribute values to the existing values
            cur_attrs = attrs[feature_id]
            new_attrs = []
            for cur, new in zip(cur_attrs, row_attrs):
                attr_key = cur[0]
                updated_vals = set(cur[1] + new[1])
                new_attrs.append((attr_key, tuple(updated_vals)))

            attrs[feature_id] = new_attrs
        else:
            attrs[feature_id] = row_attrs

    # BEGIN main routine
    for file, preferred_ids in gff_files.items():
        gff = HTSeq.GFF_Reader(file)
        for row in gff:
            if row.iv.strand == ".":
                raise ValueError(f"Feature {row.name} in {file} has no strand information.")

            try:
                feature_id = row.attr["ID"][0]
                # Climb ancestral tree to find root feature ID
                feature_id = check_ancestors(feature_id, row)
                # Add feature_id <-> feature_interval records
                row_attrs = add_feature_iv(feature_id, row)
                # Append alias to feat if unique
                add_alias(feature_id, row.attr)
            except KeyError as ke:
                raise ValueError(f"Feature {row.name} does not contain a {ke} attribute in {file}")

            # Add feature_id -> feature_attributes record
            incorporate_attributes(feature_id, row_attrs)

    return feats, attrs, alias, intervals
