import numpy as np
import HTSeq
import sys
import re

from typing import Tuple, List, Dict
from aquatx.srna.util import report_execution_time

# For parse_GFF_attribute_string()
# Todo: I believe _re_attr_main may fail if user GFFs have escape characters, which are valid per GFF3 specification
_re_attr_main = re.compile(r"\s*([^\s=]+)[\s=]+(.*)")
_re_attr_empty = re.compile(r"^\s*$")

# Type aliases for human readability
Features = HTSeq.GenomicArrayOfSets  # interval -> set of associated features
Attributes = Dict[str, list]  # feature -> feature attributes
FeatureSources = Dict[str, list]
SelectionRules = List[dict]
Alias = dict


class Alignment:
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
    It has been adapted to allow features to have multiple classes, which are
    split and stored as a tuple rather than a comma separated string. For
    downstream compatibility with membership operations
    (e.g. (attr_key, attr_val) in feature_candidate), non-list attribute values
    are recorded as tuples.

    If 'extra_return_first_value' is set, a pair is returned: the dictionary
    and the value of the first attribute. This might be useful if this is the
    ID.
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
def build_reference_tables(gff_files: FeatureSources, rules: SelectionRules) -> Tuple[Features, Attributes, Alias]:
    """A simplified and slightly modified version of HTSeq.create_genomicarrayofsets

    This modification changes the information stored in an interval's step vector
    within the features table. This stores (feature ID, feature type) tuples rather
    than the feature ID alone. It also allows for cataloguing features by any attribute,
    not just by ID, on a per GFF file basis.

    Note: at this time if the same feature is defined in multiple GFF files using
    different ID attributes, the feature will
    """

    # Patch the GFF attribute parser to support comma separated attribute value lists
    setattr(HTSeq.features, 'parse_GFF_attribute_string', parse_GFF_attribute_string)

    # Obtain an ordered list of unique attributes of interest from selection rules
    attrs_of_interest = list(np.unique(["Class"] + [attr['Identity'][0] for attr in rules]))

    feats = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    attrs, alias, intervals = {}, {}, {}

    for file, preferred_ids in gff_files.items():
        gff = HTSeq.GFF_Reader(file)
        for row in gff:
            if row.iv.strand == ".":
                raise ValueError(f"Feature {row.name} in {file} has no strand information.")

            try:
                feature_id = row.attr["ID"][0]

                # If this feature_id has already been seen, make sure their intervals match
                if feature_id in intervals and row.iv != intervals[feature_id]:
                    raise ValueError(f"{feature_id} is being redefined under a different interval in {file}.")

                # Add feature_id <-> feature_interval records
                feats[row.iv] += feature_id
                intervals[feature_id] = row.iv
                row_attrs = [(interest, row.attr[interest]) for interest in attrs_of_interest]

                curr_alias = alias.get(feature_id, ())
                for pref_id in preferred_ids:
                    if pref_id != "ID" and row.attr[pref_id] not in curr_alias:
                        # Add feature_id -> feature_alias_tuple record
                        # Append to feature's aliases if it does not already contain
                        alias[feature_id] = curr_alias + row.attr[pref_id]
            except KeyError as ke:
                raise ValueError(f"Feature {row.name} does not contain a {ke} attribute in {file}")

            if feature_id in attrs and row_attrs != attrs[feature_id]:
                # If an attribute record already exists for this feature, and this row provides new attributes,
                #  append the new attribute values to the existing values
                cur_attrs = attrs[feature_id]
                row_attrs = [(cur[0], cur[1] + new[1]) for cur, new in zip(cur_attrs, row_attrs)]

            # Add feature_id -> feature_attributes record
            attrs[feature_id] = row_attrs

    return feats, attrs, alias