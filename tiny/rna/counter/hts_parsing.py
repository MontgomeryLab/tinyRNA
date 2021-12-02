import numpy as np
import HTSeq
import sys
import os
import re

from line_profiler_pycharm import profile

from collections import Counter, defaultdict
from typing import Tuple, List, Dict, Union, Set

from tiny.rna.util import report_execution_time

# For parse_GFF_attribute_string()
_re_attr_main = re.compile(r"\s*([^\s=]+)[\s=]+(.*)")
_re_attr_empty = re.compile(r"^\s*$")

# For Alignment
complement = {ord('A'): 'T', ord('T'): 'A', ord('G'): 'C', ord('C'): 'G'}

class Alignment:
    """The data structure in which parsed SAM alignments are stored.

    Strand-non-specific 5' end nucleotide is stored for efficient lookup. This
    allows us to skip performing full reverse complement of sequences aligned
    to the antisense strand.
    """

    def __init__(self, name, chrom, start, seq, strand):
        nt5 = complement[seq[-1]] if strand == '-' else chr(seq[0])
        length = len(seq)
        self.name = name
        self.len = length
        self.seq = seq
        self.nt5 = nt5

        self.chrom = chrom
        self.start = start
        self.end = start + length
        self.strand = strand

    def __repr__(self):
        return f"<Alignment Object: Read '{self.name}' {self.seq} ({self.len} bases)>"

    def __len__(self):
        return self.len


def read_SAM(file):
    """A minimal SAM parser which only handles data relevant to the workflow, for performance."""

    with open(file, 'rb') as f:
        line = f.readline()

        # Skip headers
        while line[0] == ord('@'):
            line = f.readline()

        while line:
            cols = line.split(b'\t')
            line = f.readline()

            seq = cols[9]
            start = int(cols[3]) - 1
            length = len(seq)

            # Note: we assume sRNA sequencing data is NOT reversely stranded
            if (int(cols[1]) & 16):
                strand = '-'
                nt5 = complement[seq[-1]]
            else:
                strand = '+'
                nt5 = chr(seq[0])

            yield {
                "name": cols[0].decode('utf-8'),
                "len": length,
                "seq": seq,
                "nt5": nt5,
                "chrom": cols[2].decode('utf-8'),
                "start": start,
                "end": start + length,
                "strand": strand
            }


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
    It has been adapted to parse comma separated attribute values as separate values.
    Values are stored in a set for ease of handling in ReferenceTables and because
    duplicate values don't make sense in this context.

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
        attribute_dict[sys.intern(key)] = {sys.intern(val)} \
            if ',' not in val \
            else set(c.strip() for c in val.split(','))
        if extra_return_first_value and i == 0:
            first_val = val
    if extra_return_first_value:
        return attribute_dict, first_val
    else:
        return attribute_dict


class ReferenceTables:
    """A GFF parser which builds feature, identity, and alias tables, with intelligent appends

    Features may be defined by multiple GFF files. If multiple files offer different identities for
    the same feature, the unique among those attribute keys and values are merged with the record.
    If multiple aliases (or Name Attributes, per the Features Sheet) are defined for a feature, the
    unique among those names are appended. Each GFF file defined in the Features Sheet is parsed only
    once regardless of the number of Name Attributes associated with it.

    Features which define a Parent or share an ID attribute are treated as discontinuous features.
    In these cases the root ancestor feature receives merged identities, intervals, and aliases.
    Children of the root ancestor are otherwise not stored in the reference tables.

    Source and type filters allow the user to define acceptable values for columns 2 and 3 of the
    GFF, respectively. These filters are inclusive (only rows with matching values are parsed),
    and behave as a logical AND if both are defined. Empty filter lists allow all matches.
    Feature lineage is preserved even for filtered features; for these cases, the root ancestor
    is considered to be the highest unfiltered feature in the lineage.
    """

    source_filter = []
    type_filter = []

    def __init__(self, gff_files: Dict[str, list], rules: Dict[str, list], inv_ident, **kwargs):
        self.gff_files = gff_files
        self._set_filters(**kwargs)
        self.all_features = kwargs.get('all_features', False)
        # self.attrs_of_interest = defaultdict(set)
        # for rule in rules:
        #     self.attrs_of_interest[rule['Identity'][0]].add(rule['Identity'][1])
        self.matchies = defaultdict(set)
        self.inv_ident = inv_ident
        self.rules = rules

        import CyStep._stepvector as StepVector
        setattr(HTSeq._HTSeq, "StepVector", StepVector)

        # self.feats = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        self.feats = HTSeq.GenomicArray("auto", typecode='O', stranded=True)
        self.ident, self.alias = defaultdict(lambda: defaultdict(set)), defaultdict(set)
        self.classes, self.parents, self.filtered = {}, {}, set()
        self.intervals = defaultdict(list)
        self.classes = defaultdict(set)

        # Patch the GFF attribute parser to support comma separated attribute value lists
        setattr(HTSeq.features, 'parse_GFF_attribute_string', parse_GFF_attribute_string)

    @report_execution_time("GFF parsing")
    def get(self) -> \
            Tuple['HTSeq.GenomicArrayOfSets', Dict[str, list], Dict[str, tuple], Dict[str, list], Dict[str, list]]:
        """Initiates GFF parsing and returns the resulting reference tables"""

        for file, alias_keys in self.gff_files.items():
            gff = HTSeq.GFF_Reader(file)
            try:
                for row in gff:
                    if row.iv.strand == ".":
                        raise ValueError(f"Feature {row.name} in {file} has no strand information.")
                    if not self.filter_match(row):
                        self.exclude_row(row)
                        continue
                    try:
                        # Grab the primary key for this feature
                        feature_id = next(iter(row.attr["ID"]))
                        # Select only identities (key-val pairs) of interest
                        idxs, classes = self.get_interesting_idents(feature_id, row.attr)
                        # Only add features with identity matches if all_features is False
                        if not self.all_features and not len(idxs):
                            self.exclude_row(row)
                            continue
                        # Add feature_id <-> feature_interval records
                        root_id = self.add_feature_iv(feature_id, row, idxs, classes)
                        # Append alias to root feature if it is unique
                        self.add_alias(alias_keys, root_id, row.attr)
                        # Add feature_id -> feature_idents record
                        #self.incorporate_idents(root_id, idents)
                    except KeyError as ke:
                        raise ValueError(f"Feature {row.name} does not contain a {ke} attribute.")
            except Exception as e:
                # Append to error message while preserving exception provenance and traceback
                e.args = (str(e.args[0]) + "\nError occurred on line %d of %s" % (gff.line_no, file),)
                raise e.with_traceback(sys.exc_info()[2]) from e

        self.finalize_tables()
        return self.feats, self.ident, self.alias, dict(self.intervals), self.classes

    def get_root_feature(self, feature_id: str, row_attrs: Dict[str, set]) -> str:
        """Returns the ID of the feature's root parent if one exists. Otherwise the original ID is returned."""

        if "Parent" not in row_attrs:
            return feature_id

        parent_id = self.get_row_parent(feature_id, row_attrs)
        tree = [feature_id, parent_id]

        # Climb ancestral tree until the root parent is found
        while parent_id in self.parents:
            parent_id = self.parents[parent_id]
            tree.append(parent_id)

        # Descend tree until the descendent is found in the identites table
        # This is because ancestor feature(s) may have been filtered
        for ancestor in tree[::-1]:
            if ancestor in self.ident or ancestor == feature_id:
                return ancestor

    def add_feature_iv(self, feature_id: str, row, idxs, classes) -> str:
        """Adds the feature and its intervals to corresponding reference tables"""

        root_id = self.get_root_feature(feature_id, row.attr)
        self.classes[root_id] |= classes
        self.matchies[root_id] |= idxs

        # self.feats[row.iv] += root_id
        self.ident.setdefault(root_id, self.ident.default_factory())
        if row.iv not in self.intervals[root_id]:
            self.intervals[root_id].append(row.iv)

        return root_id

    def add_alias(self, alias_keys, root_id, row_attr):
        """Merge unique aliases with the root feature's"""

        for alias_key in alias_keys:
            for row_val in row_attr[alias_key]:
                self.alias[root_id].add(row_val)

    def get_interesting_idents(self, feat_id: str, row_attrs: Dict[str, set]) -> Set[int]:
        """Returns only the identities of interest from the row's identities"""

        interest_matches = set()
        # for interest in self.attrs_of_interest:
        #     match = row_attrs[interest] & self.attrs_of_interest[interest]
        #     if match: interest_matches[interest] = match
        for ident, rule_idxs in self.inv_ident.items():
            match = row_attrs.get(ident[0], None)
            if match is not None and ident[1] in match:
                interest_matches.update(r for r in rule_idxs)

        # {(rule, rank, strict), ...}
        interest_matches = {(r, self.rules[r]['Hierarchy'], self.rules[r]['Strict']) for r in interest_matches}
        classes = {c for c in row_attrs["Class"]}

        return interest_matches, classes

    def incorporate_idents(self, root_id, row_attrs):
        """Add unique identities (keys and values) to root feature's identities"""

        for key in row_attrs:
            self.ident[root_id][key] |= row_attrs[key]

    def get_row_parent(self, feature_id: str, row_attrs: Dict[str, set]) -> str:
        """Get the current feature's parent while cooperating with filtered features"""

        parent_attr = row_attrs.get("Parent", {None})
        parent = next(iter(parent_attr))

        if len(parent_attr) > 1:
            raise ValueError(f"{feature_id} defines multiple parents which is unsupported at this time.")
        if len(parent_attr) == 0 or parent is None:
            return feature_id
        if (parent not in self.ident                # If parent is not a root feature
                and parent not in self.parents      # If parent doesn't have a parent itself
                and parent not in self.filtered):   # If parent was not a filtered root feature
            raise ValueError(f"Feature ID {parent} is referenced as a parent before being defined. Please "
                             "ensure that it occurs before children in your GFF file.")

        self.parents[feature_id] = parent
        return parent

    def exclude_row(self, row):
        """The current row was filtered, but we still need to account for it and its parent

        We don't want to add filtered features to the identities table because we later call
        get_keys() on the table to obtain a list of all features considered for counting.
        We still need to keep track of parents for situations where an ancestral tree has
        a gap due to filtering; in this case we still want to merge descendents with the
        highest considered feature in the tree."""

        feature_id = row.attr['ID'].pop()
        self.filtered.add(feature_id)
        if "Parent" in row.attr:
            self.parents[feature_id] = self.get_row_parent(feature_id, row.attr)
        if row.iv.chrom not in self.feats.chrom_vectors:
            self.feats.add_chrom(row.iv.chrom)

    def finalize_tables(self):
        """Internally these are sets for ease, but Counter expects tuples for speed and hashability"""

        self.classes = {feat: tuple(sorted(classes)) for feat, classes in self.classes.items()}
        self.alias = {feat: tuple(sorted(aliases)) for feat, aliases in self.alias.items()}
        self.ident = {feat: [(select_for, with_value)
                             for select_for in self.ident[feat]
                             for with_value in self.ident[feat][select_for]]
                      for feat in self.ident}

        # (root_id, [(rule, rank, strict), ...])
        for root_id, ivs in self.intervals.items():
            matchies = tuple(sorted(self.matchies[root_id], key=lambda x: x[1]))
            for iv in ivs:
                self.feats[iv] += {(root_id, matchies)}

    @classmethod
    def _set_filters(cls, **kwargs):
        """Assigns inclusive filter values"""

        for filt in ["source_filter", "type_filter"]:
            setattr(cls, filt, kwargs.get(filt, []))

    @classmethod
    def filter_match(cls, row):
        """Checks if the GFF row passes the inclusive filter(s)

        If both filters are defined then the must both evaluate to true for a match"""

        select = True
        if len(cls.source_filter):
            select &= row.source in cls.source_filter
        if len(cls.type_filter):
            select &= row.type in cls.type_filter
        return select