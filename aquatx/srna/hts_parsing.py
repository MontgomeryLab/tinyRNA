import numpy as np
import HTSeq
import time
import sys
import re

from typing import Tuple, Set, List, Dict

# For parse_GFF_attribute_string()
# Todo: I believe _re_attr_main may fail if user GFFs have escape characters, which are valid per GFF3 specification
_re_attr_main = re.compile(r"\s*([^\s=]+)[\s=]+(.*)")
_re_attr_empty = re.compile(r"^\s*$")

# Type aliases for human readability
Features = HTSeq.GenomicArrayOfSets  # interval -> set of associated features
Attributes = Dict[str, list]  # feature -> feature attributes
FeatureSources = Set[Tuple[str, str]]
SelectionRules = List[dict]
Alias = dict


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
            return f"<Sequence Object: '{self.name}', {self.seq} ({self.len} bases)"

        def __len__(self):
            return self.len

    def __init__(self, iv, name, seq):
        nt5 = self.complement[seq[-1]] if iv.strand == '-' else chr(seq[0])
        self.read = self.Sequence(name, seq, nt5)
        self.iv = iv

    def __repr__(self):
        return f"<Alignment Object: Read '{self.read.name}' aligned to {self.iv}"


def read_SAM(file):
    """A minimal SAM parser which only handles data relevant to the workflow, for performance."""

    with open(file, 'rb') as f:
        line = f.readline()

        # Skip headers
        while line[0] == ord('@'):
            line = f.readline()

        while line:
            cols = line.split(b'\t')

            strand = "+" if (int(cols[1]) & 16) >> 4 == 0 else "-"
            chrom = cols[2].decode('utf-8')
            name = cols[0].decode('utf-8')
            start = int(cols[3]) - 1
            seq = cols[9]

            iv = HTSeq.GenomicInterval(chrom, start, start + len(seq), strand)
            yield Alignment(iv, name, seq)

            line = f.readline()


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


def build_reference_tables(gff_files: FeatureSources, rules: SelectionRules) -> Tuple[Features, Attributes, Alias]:
    """A GFF parser which builds feature, attribute, and alias tables, with intelligent appends

    Features may be defined by multiple GFF files. If multiple files offer different attributes for
    the same feature, the unique among those attributes are appended to the record. If multiple aliases
    (or Name Attributes, per the Features Sheet) are defined for a feature, the unique among those
    names are appended. Each GFF file defined in the Features Sheet is parsed only once regardless of
    the number of Name Attributes associated with it. Each feature ID may be defined for only one
    interval; if multiple interval definitions are supplied, then these feature IDs are renamed on
    the basis of their source GFF file.
    """

    start_time = time.time()

    # Patch the GFF attribute parser to support comma separated attribute value lists
    setattr(HTSeq.features, 'parse_GFF_attribute_string', parse_GFF_attribute_string)

    # Obtain an ordered list of unique attributes of interest from selection rules
    attrs_of_interest = list(np.unique(["Class"] + [attr['Identity'][0] for attr in rules]))

    feats = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    attrs, alias = {}, {}

    for file, preferred_id in gff_files:
        gff = HTSeq.GFF_Reader(file)
        for row in gff:
            if row.iv.strand == ".":
                raise ValueError(f"Feature {row.name} in {file} has no strand information.")

            try:
                # Add feature_id -> feature_interval record
                feature_id = row.attr["ID"][0]
                feats[row.iv] += feature_id
                row_attrs = [(interest, row.attr[interest]) for interest in attrs_of_interest]

                if preferred_id != "ID":
                    # Add feature_id -> feature_alias_tuple record
                    # If an alias already exists for this feature, append to feature's aliases
                    alias[feature_id] = alias.get(feature_id, ()) + row.attr[preferred_id]
            except KeyError as ke:
                raise ValueError(f"Feature {row.name} does not contain a {ke} attribute in {file}")

            # Todo: ensure the second part of this condition is sound (since list, shouldn't check "in" instead of ==?)...
            if feature_id in attrs and row_attrs != attrs[feature_id]:
                # If an attribute record already exists for this feature, and this row provides new attributes,
                #  append the new attribute values to the existing values
                cur_attrs = attrs[feature_id]
                row_attrs = [(cur[0], cur[1] + new[1]) for cur, new in zip(cur_attrs, row_attrs)]

            # Add feature_id -> feature_attributes record
            attrs[feature_id] = row_attrs

    print("GFF parsing took %.2f seconds" % (time.time() - start_time))
    return feats, attrs, alias