import os.path
import HTSeq
import sys
import re

from collections import Counter, defaultdict
from typing import Tuple, List, Dict, Iterator, Optional, DefaultDict, Set, Union
from inspect import stack

from tiny.rna.counter.matching import Wildcard
from tiny.rna.util import report_execution_time, make_filename

# For parse_GFF_attribute_string()
_re_attr_main = re.compile(r"\s*([^\s=]+)[\s=]+(.*)")
_re_attr_empty = re.compile(r"^\s*$")

# For Alignment
complement = {ord('A'): 'T', ord('T'): 'A', ord('G'): 'C', ord('C'): 'G'}


class SAM_reader:
    """A minimal SAM reader that bundles multiple-alignments and only parses data relevant to the workflow"""

    def __init__(self, **prefs):
        self.decollapse = prefs.get("decollapse", False)
        self.out_prefix = prefs.get("out_prefix", None)
        self._decollapsed_filename = None
        self._decollapsed_reads = []
        self._headers = []
        self.file = None

    def bundle_multi_alignments(self, file) -> Iterator[List[dict]]:
        """Bundles multiple alignments by name"""

        self.file = file
        sam_parser = self._parse_alignments
        with open(file, 'rb') as f:
            aln_iter = iter(sam_parser(f))
            bundle = [next(aln_iter)]
            for aln in aln_iter:
                if aln['name'] != bundle[0]['name']:
                    yield bundle
                    bundle = [aln]
                else:
                    bundle.append(aln)
            yield bundle

        if self.decollapse and len(self._decollapsed_reads):
            self._write_decollapsed_sam()

    def _parse_alignments(self, file_obj) -> Iterator[dict]:
        """Parses and yields individual SAM alignments from the open file_obj"""

        line = self._read_thru_header(file_obj)
        line_no = len(self._headers)
        decollapse_sam = self.decollapse

        try:
            while line:
                line_no += 1
                cols = line.split(b'\t')

                if decollapse_sam:
                    self._decollapsed_reads.append((cols[0], line))
                    if len(self._decollapsed_reads) > 100000:
                        self._write_decollapsed_sam()

                line = file_obj.readline()  # Next line
                start = int(cols[3]) - 1
                seq = cols[9]
                length = len(seq)

                # Note: we assume sRNA sequencing data is NOT reversely stranded
                if (int(cols[1]) & 16):
                    strand = False  # -
                    try:
                        nt5 = complement[seq[-1]]
                    except KeyError:
                        nt5 = chr(seq[-1])
                else:
                    strand = True  # +
                    nt5 = chr(seq[0])

                yield {
                    "name": cols[0].decode(),
                    "Length": length,
                    "seq": seq,
                    "nt5end": nt5,
                    "chrom": cols[2].decode(),
                    "start": start,
                    "end": start + length,
                    "Strand": strand
                }
        except Exception as e:
            # Append to error message while preserving exception provenance and traceback
            e.args = (str(e.args[0]) + '\n' + f"Error occurred on line {line_no} of {self.file}",)
            raise e.with_traceback(sys.exc_info()[2]) from e

    def _get_decollapsed_filename(self):
        if self._decollapsed_filename is None:
            basename = os.path.splitext(os.path.basename(self.file))[0]
            self._decollapsed_filename = make_filename([basename, "decollapsed"], ext='.sam')
        return self._decollapsed_filename

    def _read_thru_header(self, file_obj):
        """Advance file_obj past the SAM header and return the first alignment unparsed"""

        line = file_obj.readline()
        while line[0] == ord('@'):
            self._headers.append(line.decode('utf-8'))
            line = file_obj.readline()

        if self.decollapse:
            # Write the same header data to the decollapsed file
            with open(self._get_decollapsed_filename(), 'w') as f:
                f.writelines(self._headers)

        return line

    def _write_decollapsed_sam(self):
        aln_out, prevname, seq_count = [], None, 0
        for name, line in self._decollapsed_reads:
            # Parse count just once per multi-alignment
            if name != prevname:
                seq_count = int(name.split(b"=")[1])

            aln_out.extend([line] * seq_count)
            prevname = name

        with open(self._get_decollapsed_filename(), 'ab') as sam_o:
            sam_o.writelines(aln_out)
            self._decollapsed_reads.clear()


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


def parse_GFF_attribute_string(attrStr, extra_return_first_value=False, gff_version=2):
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

    if attrStr.endswith("\n"):
        attrStr = attrStr[:-1]

    # Modification: store attributes in a dict subclass that allows case-insensitive ops
    attribute_dict = CaseInsensitiveAttrs()
    first_val = "_unnamed_"

    if gff_version == 2:
        iterator = HTSeq._HTSeq.quotesafe_split(attrStr.encode())
    else:
        # GFF3 does not care about quotes
        iterator = attrStr.encode().split(b';')

    for i, attr in enumerate(iterator):
        attr = attr.decode()
        if _re_attr_empty.match(attr):
            continue
        if (gff_version == 2) and attr.count('"') not in (0, 2):
            raise ValueError(
                "The attribute string seems to contain mismatched  quotes.")
        mo = _re_attr_main.match(attr)
        if not mo:
            raise ValueError("Failure parsing GFF attribute line")
        key = mo.group(1)
        val = mo.group(2)
        if (gff_version == 2) and val.startswith('"') and val.endswith('"'):
            val = val[1:-1]
        # Modification: allow for comma separated attribute values
        attribute_dict[key] = (val,) \
            if ',' not in val \
            else tuple(c.strip() for c in val.split(','))
        if extra_return_first_value and i == 0:
            first_val = val
    if extra_return_first_value:
        return attribute_dict, first_val
    else:
        return attribute_dict


class CaseInsensitiveAttrs(Dict[str, tuple]):
    """A dictionary subclass that allows for case-insensitive queries against feature attributes

    From a bird's eye view, this class holds the feature attribute's name as the key,
    and a tuple of values associated with that attribute name. The attribute's values did
    not contain any commas, this tuple will have a length of 1. Otherwise, the original value
    is tokenized on comma and each token is stored in a separate tuple index.

    Internally, each key is stored in lowercase form, and its associated value is a (nested) tuple
    that contains (at the following indices):
        [0]: The key (attribute name) in its original case
        [1]: A tuple of values in their original case
        [2]: A tuple of values in lowercase form

    Interactions with the Dict base class involve handling the "internal" nested tuple described
    above. Functions which call self[item] in turn call these methods, and therefore the handling of the
    internal tuple is abstracted away; these functions only deal with the key/values in their original form.
    """

    def __init__(self):
        super().__init__()

    def __setitem__(self, key: str, vals: tuple):
        lowercase_vals = tuple(map(str.lower, vals))
        super().__setitem__(key.lower(), (key, vals, lowercase_vals))

    def __getitem__(self, key: str):
        # Allows case-insensitive key lookups which return original case values
        # Ensure that KeyError contains the original key
        if key.lower() not in self: raise KeyError(key)
        return super().__getitem__(key.lower())[1]

    def __contains__(self, key: str):
        # Allows case-insensitive membership queries by key
        return super().__contains__(key.lower())

    def __str__(self):
        # Returns original case keys/values
        return str({v[0]: v[1] for v in super().values()})

    def __repr__(self):
        # Returns both original case and lowercase keys/values
        return str({f"{k}/{v[0]}": f"{v[2]}/{v[1]}" for k,v in super().items()})

    def setdefault(self, key: str, value: Optional[Tuple]=None):
        if key not in self:
            self[key] = value
            return value
        else:
            return self[key]

    def get(self, key: str, default=None):
        # Same as __getitem__() but with a default return value
        return self[key] if key in self else default

    def keys(self):
        # Roughly mimics a KeysView with original case
        for v in super(CaseInsensitiveAttrs, self).values():
            yield v[0]

    def values(self):
        # Roughly mimics a ValuesView with original case
        for v in super(CaseInsensitiveAttrs, self).values():
            yield v[1]

    def items(self):
        # Roughly mimics an ItemsView with original case
        yield from zip(self.keys(), self.values())

    # Non-overridden method
    def contains_ident(self, query: Tuple[Union[str, Wildcard], Union[str, Wildcard]]):
        """Checks if the identity tuple is present in the dictionary. Wildcards are supported."""

        key_type = type(query[0])
        val_type = type(query[1])

        key = query[0].lower() if key_type is str else query[0]
        val = query[1].lower() if val_type is str else query[1]

        if key_type is val_type is Wildcard:
            return True
        if key_type is val_type is str:
            # Allows case-insensitive membership queries by (key, value)
            return key in self and \
                   val in super(CaseInsensitiveAttrs, self).__getitem__(key)[2]
        if key_type is str and val_type is Wildcard:
            return key in self
        if key_type is Wildcard and val_type is str:
            for v in super(CaseInsensitiveAttrs, self).values():
                if val in v[2]: return True
            else: return False

    # Dict methods not implemented which are invalid if delegated to dict class
    def update(self, other, **kwargs): self.not_implemented()
    def fromkeys(self, it, val=None): self.not_implemented()
    def popitem(self): self.not_implemented()
    def copy(self): self.not_implemented()

    def not_implemented(self):
        raise NotImplementedError(f"CaseInsensitiveAttrs does not support {stack()[1].function}")


# Type aliases for human readability
ClassTable = AliasTable = DefaultDict[str, Tuple[str]]
StepVector = HTSeq.GenomicArrayOfSets
Tags = DefaultDict[str, Set[str]]


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
    They are structured as (rank, rule, overlap). "Rank" is the heirarchy value of the matching
    rule, "rule" is the index of that rule in FeatureSelector's rules table, and "overlap" is the
    IntervalSelector per the rule's overlap requirements.

    Source and type filters allow the user to define acceptable values for columns 2 and 3 of the
    GFF, respectively. These filters are inclusive (only rows with matching values are parsed),
    and behave as a logical AND if both are defined. Empty filter lists allow all matches.
    Feature lineage is preserved even for filtered features; for these cases, the root ancestor
    is considered to be the highest unfiltered feature in the lineage.
    """

    source_filter = []
    type_filter = []

    def __init__(self, gff_files: Dict[str, list], feature_selector, **prefs):
        self.all_features = prefs.get('all_features', False)
        self.selector = feature_selector
        self._set_filters(**prefs)
        self.gff_files = gff_files
        # ----------------------------------------------------------- Primary Key:
        if prefs['step_vector'] == 'Cython':
            from tiny.rna.counter.stepvector import StepVector
            setattr(HTSeq.StepVector, 'StepVector', StepVector)
            self.feats = HTSeq.GenomicArray("auto", stranded=False)         # Root Match ID
        else:
            self.feats = HTSeq.GenomicArrayOfSets("auto", stranded=False)   # Root Match ID

        self.parents, self.filtered = {}, set()                         # Original Feature ID
        self.intervals = defaultdict(list)                              # Root Feature ID
        self.matches = defaultdict(set)                                 # Root Match ID
        self.classes = defaultdict(set)                                 # Root Feature ID
        self.alias = defaultdict(set)                                   # Root Feature ID
        self.tags = defaultdict(set)                                    # Root Feature ID -> Root Match ID

        # Patch the GFF attribute parser to support comma separated attribute value lists
        setattr(HTSeq.features.GFF_Reader, 'parse_GFF_attribute_string', staticmethod(parse_GFF_attribute_string))

    @report_execution_time("GFF parsing")
    def get(self) -> Tuple[StepVector, AliasTable, ClassTable, dict]:
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

                    # Grab the primary key for this feature
                    feature_id = self.get_feature_id(row)
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
            except Exception as e:
                # Append to error message while preserving exception provenance and traceback
                e.args = (str(e.args[0]) + "\nError occurred on line %d of %s" % (gff.line_no, file),)
                raise e.with_traceback(sys.exc_info()[2]) from e

        return self.finalize_tables()

    def get_root_feature(self, lineage: list) -> str:
        """Returns the highest feature ID in the ancestral tree which passed stage 1 selection.
        The original feature ID is returned if there are no valid ancestors."""

        # Descend tree until the descendent is found in the matches table
        # This is because ancestor feature(s) may have been filtered
        for ancestor in lineage[::-1]:
            if self.was_matched(ancestor):
                return ancestor
        else:
            return lineage[0]  # Default: the original feature_id

    def get_feature_ancestors(self, feature_id: str, row_attrs: CaseInsensitiveAttrs):
        if "Parent" not in row_attrs:
            return [feature_id]

        parent_id = self.get_row_parent(feature_id, row_attrs)
        lineage = [feature_id, parent_id]

        # Climb ancestral tree until the root parent is found
        while parent_id in self.parents:
            parent_id = self.parents[parent_id]
            lineage.append(parent_id)

        return lineage

    def add_feature(self, feature_id: str, row, matches: defaultdict, classes: set) -> str:
        """Adds the feature to classes and intervals tables under its root ID, and to the matches table
        under its tagged ID. Note: matches are later assigned to intervals in _finalize_features()"""

        lineage = self.get_feature_ancestors(feature_id, row.attr)
        root_id = self.get_root_feature(lineage)

        self.classes[root_id] |= classes
        if row.iv not in self.intervals[root_id]:
            self.intervals[root_id].append(row.iv)

        if matches:
            for tag, matches in matches.items():
                tagged_id = (root_id, tag)
                self.tags[root_id].add(tagged_id)
                self.matches[tagged_id] |= matches
        else:
            # Features without identity matches are saved if self.all_features
            assert self.all_features is True, "Feat. without identity matches was saved but self.all_features is False"
            self.tags[root_id].add((root_id, ''))
            self.matches[(root_id, '')] |= set()

        return root_id

    def add_alias(self, root_id: str, alias_keys: List[str], row_attrs: CaseInsensitiveAttrs) -> None:
        """Add feature's aliases to the root ancestor's alias set"""

        for alias_key in alias_keys:
            for row_val in row_attrs.get(alias_key, ()):
                self.alias[root_id].add(row_val)

    def get_matches_and_classes(self, row_attrs: CaseInsensitiveAttrs) -> Tuple[DefaultDict, set]:
        """Grabs classes and match tuples from attributes that match identity rules"""

        row_attrs.setdefault("Class", ("_UNKNOWN_",))
        classes = {c for c in row_attrs["Class"]}

        identity_matches = defaultdict(set)
        for ident, rule_indexes in self.selector.inv_ident.items():
            if row_attrs.contains_ident(ident):
                for index in rule_indexes:
                    # Non-tagged matches are pooled under '' empty string
                    tag = self.selector.rules_table[index]['Tag']
                    identity_matches[tag].add((
                        index,
                        self.selector.rules_table[index]['Hierarchy'],
                        self.selector.rules_table[index]['Overlap']
                    ))
        # -> identity_matches: {tag: (rule, rank, overlap), ...}
        return identity_matches, classes

    def get_row_parent(self, feature_id: str, row_attrs: CaseInsensitiveAttrs) -> str:
        """Get the current feature's parent while cooperating with filtered features"""

        parent_attr = row_attrs.get("Parent", [None])
        parent = parent_attr[0]

        if len(parent_attr) > 1:
            raise ValueError(f"{feature_id} defines multiple parents which is unsupported at this time.")
        if len(parent_attr) == 0 or parent is None:
            return feature_id
        if (parent not in self.tags                 # If parent is not a root feature
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

        feature_id = self.get_feature_id(row)
        self.filtered.add(feature_id)
        if "Parent" in row.attr:
            self.parents[feature_id] = self.get_row_parent(feature_id, row.attr)
        self.chrom_vector_setdefault(row.iv.chrom)

    def finalize_tables(self) -> Tuple[StepVector, AliasTable, ClassTable, dict]:
        """Convert sets to sorted tuples for performance, hashability, and deterministic outputs
        Matches must also be assigned to intervals now that all intervals are known"""

        self._finalize_classes()
        self._finalize_aliases()
        self._finalize_features()

        if self.get_feats_table_size() == 0 and self.all_features is False:
            raise ValueError("No features were retained while parsing your GFF file.\n"
                             "This may be due to a lack of features matching 'Select for...with value...'")

        return self.feats, self.alias, self.classes, self.tags

    def _finalize_features(self):
        """Assigns matches to their corresponding intervals by populating StepVector with match tuples"""

        for root_id, unmerged_sub_ivs in self.intervals.items():
            merged_sub_ivs = self._merge_adjacent_subintervals(unmerged_sub_ivs)

            for tagged_id in self.tags[root_id]:
                # Optimization opportunity: for features whose interval selectors are all "partial",
                # eliminate all but the highest ranking match in each feature family to reduce
                # loop count in stage 2 selection.

                # Non-matching feature IDs are tagged if all_features = True
                if not self.matches[tagged_id]:
                    self.chrom_vector_setdefault(merged_sub_ivs[0].chrom)
                    continue

                # Sort match tuples by rank for more efficient feature selection
                sorted_matches = sorted(self.matches[tagged_id], key=lambda x: x[1])

                for sub_iv in merged_sub_ivs:
                    finalized_match_tuples = self.selector.build_interval_selectors(sub_iv, sorted_matches.copy())
                    self.feats[sub_iv] += (tagged_id, sub_iv.strand == '+', tuple(finalized_match_tuples))

    @staticmethod
    def _merge_adjacent_subintervals(unmerged_sub_ivs: List[HTSeq.GenomicInterval]) -> list:
        """Combines overlapping and adjacent elements in the provided list of subintervals"""

        # Sort intervals so that adjacencies are adjacent by index
        unmerged_sub_ivs = iter(sorted(unmerged_sub_ivs, key=lambda x: x.start))
        continuous_iv = next(unmerged_sub_ivs).copy()
        merged_ivs = []

        # Merge adjacent intervals for discontinuous features
        for sub_iv in unmerged_sub_ivs:
            if continuous_iv.overlaps(sub_iv):
                continuous_iv.extend_to_include(sub_iv)
            else:
                merged_ivs.append(continuous_iv)
                continuous_iv = sub_iv.copy()
        merged_ivs.append(continuous_iv)

        return merged_ivs

    def _finalize_aliases(self):
        self.alias = {feat: tuple(sorted(aliases, key=str.lower)) for feat, aliases in self.alias.items()}

    def _finalize_classes(self):
        self.classes = {feat: tuple(sorted(classes, key=str.lower)) for feat, classes in self.classes.items()}

    def get_feats_table_size(self) -> int:
        """Returns the sum of features across all chromosomes and strands"""

        total_feats = 0
        empty_size = 3  # steps: start, empty set, end
        for chrom in self.feats.chrom_vectors:
            for strand in self.feats.chrom_vectors[chrom]:
                total_feats += self.feats.chrom_vectors[chrom][strand].array.num_steps() - empty_size

        return total_feats

    def was_matched(self, untagged_id):
        """Checks if the feature ID previously matched on identity, regardless of whether
        the matching rule was tagged or untagged."""

        # any() will short circuit on first match when provided a generator function
        return any(tagged_id in self.matches for tagged_id in self.tags.get(untagged_id, ()))

    def chrom_vector_setdefault(self, chrom):
        """Behaves like dict.setdefault() for chrom_vectors. Even though chrom_vectors are
        dictionaries themselves, calling setdefault on them will break the GenomicArrayOfSets"""

        if chrom not in self.feats.chrom_vectors:
            self.feats.add_chrom(chrom)

    @staticmethod
    def get_feature_id(row):
        id_collection = row.attr.get("ID", None)

        if id_collection is None:
            raise ValueError(f"Feature {row.name} does not contain an ID attribute.")
        if len(id_collection) == 0:
            raise ValueError("A feature's ID attribute cannot be empty. This value is required.")
        if len(id_collection) > 1:
            err_msg = "A feature's ID attribute cannot contain multiple values. Only one ID per feature is allowed."
            raise ValueError(err_msg)

        return id_collection[0]

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
