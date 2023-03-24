import functools
import os.path
import HTSeq
import pysam
import sys
import re

from collections import Counter, defaultdict
from typing import Tuple, List, Dict, Iterator, Optional, DefaultDict, Set, Union, IO, Callable
from abc import ABC, abstractmethod
from urllib.parse import unquote
from inspect import stack

from tiny.rna.counter.matching import Wildcard
from tiny.rna.util import report_execution_time, make_filename, ReportFormatter, append_to_exception

# For parse_GFF_attribute_string()
_re_attr_main = re.compile(r"\s*([^\s=]+)[\s=]+(.*)")
_re_attr_empty = re.compile(r"^\s*$")

# For SAM_reader
AlignmentDict = Dict[str, Union[str, int, bytes]]
Bundle = Tuple[List[AlignmentDict], int]
_re_tiny = r"\d+_count=\d+"
_re_fastx = r"seq\d+_x(\d+)"


class SAM_reader:
    """A minimal SAM reader that bundles multiple-alignments and only parses data relevant to the workflow."""

    def __init__(self, **prefs):
        self.decollapse = prefs.get("decollapse", False)
        self.out_prefix = prefs.get("out_prefix", None)
        self.collapser_type = None
        self.file = None
        self._collapser_token = None
        self._decollapsed_filename = None
        self._decollapsed_reads = []
        self._header_dict = {}
        self._header = None

    def bundle_multi_alignments(self, file: str) -> Iterator[Bundle]:
        """Bundles multi-alignments (determined by a shared QNAME) and reports the associated read's count"""

        self.file = file
        pysam_reader = pysam.AlignmentFile(file)
        aln_iter = iter(self._parse_alignments(pysam_reader))
        bundle, read_count = self._new_bundle(next(aln_iter))

        for aln in aln_iter:
            if aln['Name'] != bundle[0]['Name']:
                yield bundle, read_count
                bundle, read_count = self._new_bundle(aln)
            else:
                bundle.append(aln)
        yield bundle, read_count

        if self.decollapse and len(self._decollapsed_reads):
            self._write_decollapsed_sam()

    def _new_bundle(self, aln: AlignmentDict) -> Bundle:
        """Wraps the provided alignment in a list and reports the read's count"""

        if self.collapser_type is not None:
            token = self.collapser_token
            count = int(aln['Name'].split(token)[-1])
        else:
            count = 1

        return [aln], count

    def _parse_alignments(self, reader: pysam.AlignmentFile) -> Iterator[AlignmentDict]:
        """Yields alignment dictionaries containing relevant info from each pysam.AlignedSegment"""

        self._gather_metadata(reader)
        first_line = len(str(self._header).splitlines()) + 1
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

        for line_no, aln in enumerate(reader, start=first_line):
            try:
                if self.decollapse:
                    self._decollapsed_reads.append(aln)
                    if len(self._decollapsed_reads) > 100_000:
                        self._write_decollapsed_sam()
                if aln.is_unmapped:
                    continue

                seq = aln.query_sequence
                start = aln.pos
                length = len(seq)
                strand = aln.is_forward

                if strand:
                    nt5 = seq[0]
                else:
                    try:
                        nt5 = complement[seq[-1]]
                    except KeyError:
                        nt5 = seq[-1]

                # Note: we assume sRNA sequencing data is NOT reversely stranded

                yield {
                    "Name": aln.query_name,
                    "Length": length,
                    "Seq": seq,
                    "nt5end": nt5,
                    "Chrom": aln.reference_name,
                    "Start": start,
                    "End": start + length,
                    "Strand": strand
                }
            except Exception as e:
                # Append to error message while preserving exception provenance and traceback
                msg = f"Error occurred on line {line_no} of {self.file}"
                append_to_exception(e, msg)
                raise

    def _gather_metadata(self, reader: pysam.AlignmentFile) -> None:
        """Saves header information, determines which collapser utility (if any)
        was used before alignment, and copies the input file's header to the
        decollapsed output file if necessary."""

        header = reader.header
        first_aln = next(reader.head(1))

        self._header = header
        self._header_dict = header.to_dict()  # performs validation
        self._determine_collapser_type(first_aln.query_name)

        if self.decollapse:
            self._write_header_for_decollapsed_sam(str(header))

    def _determine_collapser_type(self, first_aln_line: str) -> None:
        """Attempts to determine the collapsing utility that was used before producing the
        input alignment file, then checks basic requirements for the utility's outputs."""

        if re.match(_re_tiny, first_aln_line) is not None:
            self.collapser_type = "tiny-collapse"

        elif re.match(_re_fastx, first_aln_line) is not None:
            self.collapser_type = "fastx"

            sort_order = self._header_dict.get('HD', {}).get('SO', None)
            if sort_order is None or sort_order != "queryname":
                raise ValueError("SAM files from fastx collapsed outputs must be sorted by queryname\n"
                                 "(and the @HD [...] SO header must be set accordingly).")
        else:
            self.collapser_type = None

            if self.decollapse:
                self.decollapse = False
                print("Alignments do not appear to be derived from a supported collapser input. "
                      "Decollapsed SAM files will therefore not be produced.", file=sys.stderr)

    def _get_decollapsed_filename(self) -> str:
        """Returns the filename to be used for writing decollapsed outputs"""

        if self._decollapsed_filename is None:
            basename = os.path.splitext(os.path.basename(self.file))[0]
            self._decollapsed_filename = make_filename([basename, "decollapsed"], ext='.sam')
        return self._decollapsed_filename

    def _write_header_for_decollapsed_sam(self, header_str) -> None:
        """Writes the provided header to the decollapsed output file"""

        assert self.collapser_type is not None
        with open(self._get_decollapsed_filename(), 'w') as f:
            f.write(header_str)

    def _write_decollapsed_sam(self) -> None:
        """Expands the collapsed alignments in the _decollapsed_reads buffer
        and writes the result to the decollapsed output file"""

        assert self.collapser_type is not None
        aln_out, prev_name, seq_count = [], None, 0
        for aln in self._decollapsed_reads:
            # Parse count just once per multi-alignment
            name = aln.query_name
            if name != prev_name:
                seq_count = int(name.split(self.collapser_token)[-1])

            aln_out.extend([aln.to_string()] * seq_count)
            prev_name = name

        with open(self._get_decollapsed_filename(), 'ab') as sam_o:
            sam_o.writelines(aln_out)
            self._decollapsed_reads.clear()

    @property
    def collapser_token(self) -> bytes:
        """Returns the split token to be used for determining read count from the QNAME field"""

        if self._collapser_token is None:
            self._collapser_token = {
                "tiny-collapse": "=",
                "fastx": "_x",
                "BioSeqZip": ":"
            }[self.collapser_type]

        return self._collapser_token


def infer_strandedness(sam_file: str, intervals: dict) -> str:
    """Infers strandedness from a sample SAM file and intervals from a parsed GFF file

    Credit: this technique is an adaptation of those in RSeQC's infer_experiment.py.
    It has been modified to accept a GFF reference file rather than a BED file,
    and to use HTSeq rather than bx-python.

    Args:
        sam_file: the path of the SAM file to evaluate
        intervals: the intervals instance attribute of ReferenceFeatures, populated by .get()
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
            strandless = HTSeq.GenomicInterval(rec['Chrom'], rec['Start'], rec['End'])
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

    This slight modification of the same method found in HTSeq.features includes
    the following for improved compliance with the GFF format:
      - Attribute values containing commas are tokenized.
      - Attribute keys are URL decoded. Values are URL decoded after tokenization.

    Per the original HTSeq docstring:
        "If 'extra_return_first_value' is set, a pair is returned: the dictionary
        and the value of the first attribute. This might be useful if this is the
        ID."
    """

    # Modification: store attributes in a dict subclass that allows case-insensitive ops
    attribute_dict = CaseInsensitiveAttrs()
    first_val = "_unnamed_"

    if gff_version == 2:
        iterator = HTSeq._HTSeq.quotesafe_split(attrStr.rstrip().encode())
    else:
        # GFF3 does not care about quotes
        iterator = attrStr.rstrip().encode().split(b';')

    for i, attr in enumerate(iterator):
        attr = attr.decode()
        if _re_attr_empty.match(attr):
            continue
        if (gff_version == 2) and attr.count('"') not in (0, 2):
            raise ValueError(
                "The attribute string seems to contain mismatched quotes.")
        mo = _re_attr_main.match(attr)
        if not mo:
            raise ValueError("Failure parsing GFF attribute line")
        key = mo.group(1)
        val = mo.group(2)
        if (gff_version == 2) and val.startswith('"') and val.endswith('"'):
            val = val[1:-1]
        # Modification: allow for comma separated attribute values
        attribute_dict[unquote(key)] = (unquote(val),) \
            if ',' not in val \
            else tuple(unquote(c.strip()) for c in val.split(','))
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


class ReferenceBase(ABC):
    """The base class for reference parsers"""

    def __init__(self, **prefs):
        self.stepvector = prefs.get('stepvector', 'HTSeq')
        self.feats = self._init_genomic_array()

        # The selector is assigned whenever get() is called.
        # While it isn't the current use case, this would allow
        # for groups of GFF files to be processed with different
        # Stage 1 selection rules and pooled into the same tables
        self.selector = None

    @abstractmethod
    def get(self, feature_selector): pass

    def _init_genomic_array(self):
        """The Cython StepVector is more efficient but requires extra setup steps.
        If these fail, we want to fall back to HTSeq's StepVector and carry on."""

        if self.stepvector == 'Cython':
            try:
                from tiny.rna.counter.stepvector import StepVector
                setattr(HTSeq.StepVector, 'StepVector', StepVector)
                return HTSeq.GenomicArray("auto", stranded=False)
            except:
                self.stepvector = 'HTSeq'
                print("Failed to import Cython StepVector\n"
                      "Falling back to HTSeq's StepVector",
                      file=sys.stderr)

        if self.stepvector == 'HTSeq':
            return HTSeq.GenomicArrayOfSets("auto", stranded=False)

    def chrom_vector_setdefault(self, chrom):
        """Behaves like dict.setdefault() for chrom_vectors. Even though chrom_vectors are
        dictionaries themselves, calling setdefault on them will break the GenomicArrayOfSets"""

        if chrom not in self.feats.chrom_vectors:
            self.feats.add_chrom(chrom)

    def get_feats_table_size(self) -> int:
        """Returns the sum of features across all chromosomes and strands"""

        total_feats = 0
        empty_size = {"Cython": 1, "HTSeq": 3}[self.stepvector]
        for chrom in self.feats.chrom_vectors:
            for strand in self.feats.chrom_vectors[chrom]:
                total_feats += self.feats.chrom_vectors[chrom][strand].array.num_steps() - empty_size

        return total_feats

    def print_selector_warnings(self):
        """Warnings accumulate in FeatureSelector and are printed here using ReportFormatter"""

        if self.selector is None: return

        header = "Incompatible feature intervals were produced as a result of overlap shift parameters."
        bad_shift = "The following matches were omitted from selection because they result in"
        descriptions = {
            "null": f"{bad_shift} a zero-width interval",
            "inverted": f"{bad_shift} an inverted interval (start > end)",
            "negative start": f"{bad_shift} a negative start coordinate"
        }

        formatter = ReportFormatter(descriptions)
        formatter.add_warning_section(header, self.selector.warnings)
        formatter.print_report()

    @staticmethod
    def map_strand(htseq_value: str):
        """Maps HTSeq's strand representation (+/-/.) to
        tinyRNA's strand representation (True/False/None)"""

        return {
            '+': True,
            '-': False,
        }.get(htseq_value, None)


def parse_gff(file, row_fn: Callable, alias_keys=None):
    if alias_keys is not None:
        row_fn = functools.partial(row_fn, alias_keys=alias_keys)

    gff = HTSeq.GFF_Reader(file)
    try:
        for row in gff:
            row_fn(row)
    except Exception as e:
        # Append to error message while preserving exception provenance and traceback
        extended_msg = f"Error occurred on line {gff.line_no} of {file}"
        append_to_exception(e, extended_msg)
        raise


# Type aliases for human readability
AliasTable = DefaultDict[str, Tuple[str]]
TagTable = DefaultDict[str, Set[Tuple[str, str]]]
GenomicArray = HTSeq.GenomicArrayOfSets


class ReferenceFeatures(ReferenceBase):
    """A GFF parser which builds feature, alias, and class reference tables

    Discontinuous features are supported, and comma-separated attribute values (GFF3 column 9)
    are treated as separate values. Multiple GFF3 files may define the same feature; only the
    unique and relevant definitions are retained. Each GFF file defined in the Features Sheet
    is parsed only once regardless of the number of Name Attributes associated with it.

    Features which define a Parent or share an ID attribute are treated as discontinuous features.
    In these cases the root ancestor feature receives merged match-tuples, classes, and aliases.
    Children of the root ancestor are otherwise not stored in the reference tables.

    Match-tuples are created for each Features Sheet rule which matches a feature's attributes.
    They are structured as (rank, rule, overlap). "Rank" is the hierarchy value of the matching
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

    def __init__(self, gff_files: Dict[str, list], **prefs):
        super().__init__(**prefs)
        self.all_features = prefs.get('all_features', False)
        self.gff_files = gff_files
        # ----------------------------------------------- Primary Key:
        # self.feats                                        # Root Match ID
        self.parents, self.filtered = {}, set()             # Original Feature ID
        self.intervals = defaultdict(list)                  # Root Feature ID
        self.matches = defaultdict(set)                     # Root Match ID
        self.alias = defaultdict(set)                       # Root Feature ID
        self.tags = defaultdict(set)                        # Root Feature ID -> Root Match ID

        # Patch the GFF attribute parser to support comma separated attribute value lists
        setattr(HTSeq.features.GFF_Reader, 'parse_GFF_attribute_string', staticmethod(parse_GFF_attribute_string))

    @report_execution_time("GFF parsing")
    def get(self, feature_selector) -> Tuple[GenomicArray, AliasTable, TagTable]:
        """Initiates GFF parsing and returns complete feature, alias, and tag tables"""

        self.selector = feature_selector
        for file, alias_keys in self.gff_files.items():
            parse_gff(file, self.parse_row, alias_keys=alias_keys)

        self._finalize_aliases()
        self._finalize_features()

        if self.selector.warnings: self.print_selector_warnings()
        if self.get_feats_table_size() == 0 and self.all_features is False:
            raise ValueError("No features were retained while parsing your GFF file.\n"
                             "This may be due to a lack of features matching 'Select for...with value...'")

        return self.feats, self.alias, self.tags

    def parse_row(self, row, alias_keys=None):
        if row.type.lower() == "chromosome":
            self.exclude_row(row)
            return

        # Perform Stage 1 selection
        matches = self.get_matches(row)
        # Skip non-matching rows unless all_features=True
        if not (len(matches) or self.all_features):
            self.exclude_row(row)
            return
        # Grab the primary key for this feature
        feature_id = self.get_feature_id(row)
        # Add feature data to root ancestor in the reference tables
        root_id = self.add_feature(feature_id, row, matches)
        # Add alias to root ancestor if it is unique
        self.add_alias(root_id, alias_keys, row.attr)

    @staticmethod
    def get_feature_id(row):
        id_collection = row.attr.get('ID') \
                        or row.attr.get('gene_id') \
                        or row.attr.get('Parent')

        if id_collection is None:
            raise ValueError(f"Feature {row.name} does not have an ID attribute.")
        if len(id_collection) == 0:
            raise ValueError("A feature's ID attribute cannot be empty. This value is required.")
        if len(id_collection) > 1:
            return ','.join(id_collection)

        return id_collection[0]

    def get_matches(self, row: HTSeq.GenomicFeature) -> DefaultDict:
        """Performs Stage 1 selection and returns match tuples under their associated classifier"""

        identity_matches = defaultdict(set)
        for ident, rule_indexes in self.selector.inv_ident.items():
            if row.attr.contains_ident(ident):
                for index in rule_indexes:
                    rule = self.selector.rules_table[index]
                    if self.column_filter_match(row, rule):
                        # Unclassified matches are pooled under '' empty string
                        identity_matches[rule['Class']].add(
                            (index, rule['Hierarchy'], rule['Overlap'])
                        )

        # -> identity_matches: {classifier: {(rule, rank, overlap), ...}, ...}
        return identity_matches

    def add_feature(self, feature_id: str, row, matches: defaultdict) -> str:
        """Adds the feature to the intervals table under its root ID, and to the matches table
        under its tagged ID. Note: matches are later assigned to intervals in _finalize_features()"""

        lineage = self.get_feature_ancestors(feature_id, row.attr)
        root_id = self.get_root_feature(lineage)

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

    def get_root_feature(self, lineage: list) -> str:
        """Returns the highest feature ID in the ancestral tree which passed stage 1 selection.
        The original feature ID is returned if there are no valid ancestors."""

        # Descend tree until the descendant is found in the matches table
        # This is because ancestor feature(s) may have been filtered
        for ancestor in lineage[::-1]:
            if self.was_matched(ancestor):
                return ancestor
        else:
            return lineage[0]  # Default: the original feature_id

    def was_matched(self, untagged_id):
        """Checks if the feature ID previously matched on identity, regardless of whether
        the matching rule was tagged or untagged."""

        # any() will short circuit on first match when provided a generator function
        return any(tagged_id in self.matches for tagged_id in self.tags.get(untagged_id, ()))

    def get_row_parent(self, feature_id: str, row_attrs: CaseInsensitiveAttrs) -> str:
        """Get the current feature's parent while cooperating with filtered features"""

        parent_attr = row_attrs.get("Parent") or [None]
        parent = parent_attr[0]

        if len(parent_attr) > 1:
            raise ValueError(f"{feature_id} defines multiple parents which is unsupported at this time.")
        if parent in (None, feature_id):
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

    def add_alias(self, root_id: str, alias_keys: List[str], row_attrs: CaseInsensitiveAttrs) -> None:
        """Add feature's aliases to the root ancestor's alias set"""

        for alias_key in alias_keys:
            for row_val in row_attrs.get(alias_key, ()):
                if row_val: self.alias[root_id].add(row_val)

    def _finalize_aliases(self):
        self.alias = {feat: tuple(sorted(aliases, key=str.lower)) for feat, aliases in self.alias.items()}

    def _finalize_features(self):
        """Assigns matches to their corresponding intervals by populating GenomicArray with match tuples

        Each feature that is handled here has matched a rule, and its matches might be classified
        into subsets depending on the rules it matched. However, the feature's interval
        is the same regardless of the match classification. This interval might be
        discontinuous, so we start by merging any of its sub-intervals that are
        adjacent to reduce loop count in Stage 2 and 3 selection. Each group of
        classified matches is then added to the GenomicArray under this shared
        interval (but the interval is not necessarily shared by matches that
        have a shift parameter in their overlap selector...)
        """

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

                tagged_matches = self.matches[tagged_id]
                self._add_subinterval_matches(tagged_id, merged_sub_ivs, tagged_matches)

        # GenomicArray is built. Clear cache...
        self.selector.overlap_cache.clear()

    def _add_subinterval_matches(self, tagged_id: tuple, sub_ivs: list, matches: set):
        """Adds the classified group of matches to the GenomicArray under each of the
        feature's sub-intervals

        These sub-intervals might be further subset if the matches define a shift
        parameter. The shift operation has to be applied to each of the feature's
        sub-intervals before the given match can be added to the GenomicArray. The
        shifted interval must match the interval that the overlap selector expects.
        """

        for sub_iv in sub_ivs:
            # Build interval selectors for this
            matches_by_shifted_iv = self.selector.build_interval_selectors(sub_iv, matches, tagged_id)
            strand = self.map_strand(sub_iv.strand)

            for shifted_iv, built_matches in matches_by_shifted_iv.items():
                # Sort match tuples by rank for more efficient feature selection
                sorted_matches = sorted(built_matches, key=lambda x: x[1])
                self.feats[shifted_iv] += (tagged_id, strand, tuple(sorted_matches))

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

    @staticmethod
    def column_filter_match(row, rule):
        """Checks if the GFF row passes the inclusive filter(s)
        If both filters are defined then they must both evaluate to true for a match"""

        return row.source in rule['Filter_s'] and row.type in rule['Filter_t']


class ReferenceSeqs(ReferenceBase):

    def __init__(self, reference_seqs, **prefs):
        super().__init__(**prefs)
        self.tags = defaultdict(set)
        self.seqs = reference_seqs
        self.alias = {}

    def get(self, selector):
        self.selector = selector
        matches = self.get_matches()

        for seq_id, seq_len in self.seqs.items():
            self.add_reference_seq(seq_id, seq_len, matches)

        # GenomicArray is built. Clear cache...
        self.selector.overlap_cache.clear()

        # If the user provided shift parameters, there will likely be warnings
        if self.selector.warnings: self.print_selector_warnings()

        # Aliases are irrelevant for non-GFF references
        aliases = {seq_id: () for seq_id in self.tags}
        return self.feats, aliases, self.tags

    def get_matches(self) -> Dict[str, List[Tuple[int, int, str]]]:
        """Stage 1 selection is skipped in sequence-based counting.
        Want the reference sequences to enter Stage 2 as though they
        were features that matched every rule in Stage 1. Simply return
        match_tuples and their classifier for all rules. These are used
        uniformly in each reference sequence's feature_record_tuple"""

        matches_by_classifier = defaultdict(list)

        for idx, rule in enumerate(self.selector.rules_table):
            match_tuple = (idx, rule['Hierarchy'], rule['Overlap'])
            matches_by_classifier[rule['Class']].append(match_tuple)

        return matches_by_classifier

    def add_reference_seq(self, seq_id, seq_len, matches_by_classifier):
        for classifier, matches in matches_by_classifier.items():
            tagged_id = (seq_id, classifier)
            self.tags[seq_id].add(tagged_id)

            for strand in ('+', '-'):
                iv = HTSeq.GenomicInterval(seq_id, 0, seq_len, strand)
                matches_by_shifted_iv = self.selector.build_interval_selectors(iv, matches, tagged_id)
                strand = self.map_strand(strand)

                for shifted_iv, built_matches in matches_by_shifted_iv.items():
                    # Sort match tuples by rank for more efficient feature selection
                    sorted_matches = sorted(built_matches, key=lambda x: x[1])
                    self.feats[shifted_iv] += (tagged_id, strand, tuple(sorted_matches))
