import HTSeq
import sys
import re

from collections import Counter, defaultdict
from typing import Tuple, List, Dict, Iterator

from tiny.rna.util import report_execution_time

# For parse_GFF_attribute_string()
_re_attr_main = re.compile(r"\s*([^\s=]+)[\s=]+(.*)")
_re_attr_empty = re.compile(r"^\s*$")

# For Alignment
complement = {ord('A'): 'T', ord('T'): 'A', ord('G'): 'C', ord('C'): 'G'}


def read_SAM(file) -> Iterator[dict]:
    """A minimal SAM reader that bundles multiple-alignments and only parses data relevant to the workflow"""

    with open(file, 'rb') as f:
        line = f.readline()

        # Skip headers
        while line[0] == ord('@'):
            line = f.readline()

        # Bundle multiple alignments by name
        aln_iter = iter(parse_alignments(f, line))
        bundle = [next(aln_iter)]
        for aln in aln_iter:
            if aln['name'] != bundle[0]['name']:
                yield bundle
                bundle = [aln]
            else:
                bundle.append(aln)
        yield bundle


def parse_alignments(f, line):

    while line:
        cols = line.split(b'\t')
        line = f.readline()

        start = int(cols[3]) - 1
        seq = cols[9]
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

    Args:
        sam_file: the path of the SAM file to evaluate
        intervals: the intervals instance attribute of ReferenceTables, populated by .get()
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
            strandless = HTSeq.GenomicInterval(rec['chrom'], rec['start'], rec['end'])
            sam_strand = rec['strand']
            gff_strand = ':'.join(unstranded[strandless].get_steps())
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
    """A GFF parser which builds feature, alias, and class reference tables

    Discontinuous features are supported, and comma-separated attribute values (GFF3 column 9)
    are treated as separate values. Multiple GFF3 files may define the same feature; only the
    unique and relevant definitions are retained. Each GFF file defined in the Features Sheet
    is parsed only once regardless of the number of Name Attributes associated with it.

    Features which define a Parent or share an ID attribute are treated as discontinuous features.
    In these cases the root ancestor feature receives merged match-tuples, classes, and aliases.
    Children of the root ancestor are otherwise not stored in the reference tables.

    Match-tuples are created for each Features Sheet rule which matches a feature's attributes.
    They are structured as (rank, rule, strict). "Rank" is the heirarchy value of the matching
    rule, "rule" is the index of that rule in FeatureSelector's rules table, and "strict" is a
    boolean representing whether a strict alignment overlap is required for a match.

    Source and type filters allow the user to define acceptable values for columns 2 and 3 of the
    GFF, respectively. These filters are inclusive (only rows with matching values are parsed),
    and behave as a logical AND if both are defined. Empty filter lists allow all matches.
    Feature lineage is preserved even for filtered features; for these cases, the root ancestor
    is considered to be the highest unfiltered feature in the lineage.
    """

    source_filter = []
    type_filter = []

    def __init__(self, gff_files: Dict[str, list], feature_selector, **kwargs):
        self.all_features = kwargs.get('all_features', False)
        self.inv_ident = feature_selector.inv_ident
        self.rules = feature_selector.rules_table
        self._set_filters(**kwargs)
        self.gff_files = gff_files

        self.feats = HTSeq.GenomicArrayOfSets("auto", stranded=False)
        self.parents, self.filtered = {}, set()
        self.intervals = defaultdict(list)
        self.matches = defaultdict(set)
        self.classes = defaultdict(set)
        self.alias = defaultdict(set)

        # Patch the GFF attribute parser to support comma separated attribute value lists
        setattr(HTSeq.features, 'parse_GFF_attribute_string', parse_GFF_attribute_string)

    @report_execution_time("GFF parsing")
    def get(self) -> Tuple['HTSeq.GenomicArray', Dict[str, tuple], Dict[str, tuple]]:
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
                        feature_id = row.attr["ID"][0]
                        # Get feature's classes and identity match tuples
                        matches, classes = self.get_matches_and_classes(row.attr)
                        # Only add features with identity matches if all_features is False
                        if not self.all_features and not len(matches):
                            self.exclude_row(row)
                            continue
                        # Add feature data to root ancestor in the reference tables
                        root_id = self.add_feature(feature_id, row, matches, classes)
                        # Add alias to root ancestor if it is unique
                        self.add_alias(root_id, alias_keys, row.attr)
                    except KeyError as ke:
                        raise ValueError(f"Feature {row.name} does not contain a {ke} attribute.")
            except Exception as e:
                # Append to error message while preserving exception provenance and traceback
                e.args = (str(e.args[0]) + "\nError occurred on line %d of %s" % (gff.line_no, file),)
                raise e.with_traceback(sys.exc_info()[2]) from e

        return self.finalize_tables()

    def get_root_feature(self, feature_id: str, row_attrs: Dict[str, tuple]) -> str:
        """Returns the ID of the feature's root parent if one exists. Otherwise the original ID is returned."""

        if "Parent" not in row_attrs:
            return feature_id

        parent_id = self.get_row_parent(feature_id, row_attrs)
        tree = [feature_id, parent_id]

        # Climb ancestral tree until the root parent is found
        while parent_id in self.parents:
            parent_id = self.parents[parent_id]
            tree.append(parent_id)

        # Descend tree until the descendent is found in the matches table
        # This is because ancestor feature(s) may have been filtered
        for ancestor in tree[::-1]:
            if ancestor in self.matches or ancestor == feature_id:
                return ancestor

    def add_feature(self, feature_id: str, row, matches: set, classes: set) -> str:
        """Adds the feature and its interval under the root ancestor's ID"""

        root_id = self.get_root_feature(feature_id, row.attr)
        self.classes[root_id] |= classes
        self.matches[root_id] |= matches

        if row.iv not in self.intervals[root_id]:
            self.intervals[root_id].append(row.iv)

        return root_id

    def add_alias(self, root_id: str, alias_keys: List[str], row_attr: Dict[str, tuple]) -> None:
        """Add feature's aliases to the root ancestor's alias set"""

        for alias_key in alias_keys:
            for row_val in row_attr[alias_key]:
                self.alias[root_id].add(row_val)

    def get_matches_and_classes(self, row_attrs: Dict[str, set]) -> Tuple[set, set]:
        """Grabs classes and match tuples from attributes that match identity rules"""

        classes = {c for c in row_attrs["Class"]}

        identity_matches = set()
        for ident, rule_indexes in self.inv_ident.items():
            match = row_attrs.get(ident[0], None)
            if match is not None and ident[1] in match:
                identity_matches.update(
                    (r, self.rules[r]['Hierarchy'], self.rules[r]['Strict'])
                    for r in rule_indexes
                )
        # -> identity_matches: {(rule, rank, strict), ...}
        return identity_matches, classes

    def get_row_parent(self, feature_id: str, row_attrs: Dict[str, tuple]) -> str:
        """Get the current feature's parent while cooperating with filtered features"""

        parent_attr = row_attrs.get("Parent", [None])
        parent = parent_attr[0]

        if len(parent_attr) > 1:
            raise ValueError(f"{feature_id} defines multiple parents which is unsupported at this time.")
        if len(parent_attr) == 0 or parent is None:
            return feature_id
        if (parent not in self.matches              # If parent is not a root feature
                and parent not in self.parents      # If parent doesn't have a parent itself
                and parent not in self.filtered):   # If parent was not a filtered root feature
            raise ValueError(f"Feature ID {parent} is referenced as a parent before being defined. Please "
                             "ensure that it occurs before children in your GFF file.")

        self.parents[feature_id] = parent
        return parent

    def exclude_row(self, row):
        """Performs residual accounting of features that were excluded from the tables

        Features may have been excluded by filters or for a lack of identity matches.
        We still need to keep track of parents for cases where these exclusions form
        gaps in the ancestral tree. This allows us to maintain a path to the tree's
        root."""

        feature_id = row.attr['ID'][0]
        self.filtered.add(feature_id)
        if "Parent" in row.attr:
            self.parents[feature_id] = self.get_row_parent(feature_id, row.attr)
        if row.iv.chrom not in self.feats.chrom_vectors:
            self.feats.add_chrom(row.iv.chrom)

    def finalize_tables(self) -> Tuple['HTSeq.GenomicArray', Dict[str, tuple], Dict[str, tuple]]:
        """Convert sets to sorted tuples for performance, hashability, and deterministic outputs"""

        self._finalize_classes()
        self._finalize_aliases()
        self._finalize_features()

        if self.get_feats_table_size() == 0 or len(self.classes) == 0:
            raise ValueError("No features or classes were retained while parsing your GFF file.\n"
                             "This may be due to a lack of features matching 'Select for...with value...'")

        return self.feats, self.alias, self.classes

    def _finalize_features(self):
        """Performs final accounting of discontinuous features and adds feature records to the StepVector"""

        for root_id, family_ivs in self.intervals.items():
            # Sort match tuples by rank for more efficient feature selection
            sorted_match_tuples = tuple(sorted(self.matches[root_id], key=lambda x: x[1]))

            # Sort intervals so that adjacencies are adjacent by index
            family_ivs.sort(key=lambda x: x.start)
            family_iter = iter(family_ivs)
            continuous = next(family_iter).copy()
            continuous_ivs = []

            # Merge adjacent intervals for discontinuous features
            for iv in family_iter:
                if continuous.overlaps(iv):
                    continuous.extend_to_include(iv)
                else:
                    continuous_ivs.append(continuous)
                    continuous = iv.copy()
            continuous_ivs.append(continuous)

            # Optimization opportunity: for match tuples with partial interval matching,
            # eliminate all but the highest ranking match in each feature family to reduce
            # loop count in phase 1 selection.

            for iv in continuous_ivs:
                self.feats[iv] += (root_id, iv.start, iv.end, iv.strand, sorted_match_tuples)

    def _finalize_aliases(self):
        self.alias = {feat: tuple(sorted(aliases)) for feat, aliases in self.alias.items()}

    def _finalize_classes(self):
        self.classes = {feat: tuple(sorted(classes)) for feat, classes in self.classes.items()}

    def get_feats_table_size(self) -> int:
        """Returns the sum of features across all chromosomes and strands"""

        total_feats = 0
        for chrom in self.feats.chrom_vectors:
            for strand in self.feats.chrom_vectors[chrom]:
                total_feats += self.feats.chrom_vectors[chrom][strand].array.num_steps()

        return total_feats

    @classmethod
    def _set_filters(cls, **kwargs):
        """Assigns inclusive filter values"""

        for filt in ["source_filter", "type_filter"]:
            setattr(cls, filt, kwargs.get(filt, []))

    @classmethod
    def filter_match(cls, row):
        """Checks if the GFF row passes the inclusive filter(s)

        If both filters are defined then they must both evaluate to true for a match"""

        select = True
        if len(cls.source_filter):
            select &= row.source in cls.source_filter
        if len(cls.type_filter):
            select &= row.type in cls.type_filter
        return select