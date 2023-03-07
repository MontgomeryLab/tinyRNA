import HTSeq
import sys

from collections import defaultdict
from typing import List, Tuple, Set, Dict, Mapping

from tiny.rna.counter.hts_parsing import SAM_reader, ReferenceFeatures, ReferenceSeqs
from .statistics import LibraryStats
from .matching import *
from ..util import append_to_exception

# Type aliases for human readability
match_tuple = Tuple[int, int, IntervalSelector]             # (rank, rule, IntervalSelector)
unbuilt_match_tuple = Tuple[int, int, str]                  # (rank, rule, interval selector keyword)
feature_record_tuple = Tuple[str, str, Tuple[match_tuple]]  # (feature ID, strand, match tuple)


class Features(metaclass=Singleton):
    chrom_vectors: HTSeq.ChromVector
    aliases: dict
    classes: dict

    def __init__(_, features: HTSeq.GenomicArrayOfSets, aliases: dict, classes: dict):
        Features.chrom_vectors = features.chrom_vectors  # For interval -> feature record tuple lookups
        Features.aliases = aliases                       # For feature ID -> preferred feature name lookups
        Features.classes = classes                       # For feature ID -> match IDs


class FeatureCounter:

    def __init__(self, references, selection_rules, **prefs):
        self.stats = LibraryStats(**prefs)
        self.sam_reader = SAM_reader(**prefs)
        self.selector = FeatureSelector(selection_rules, **prefs)

        if isinstance(references, ReferenceFeatures):
            self.mode = "by feature"
        elif isinstance(references, ReferenceSeqs):
            self.mode = "by sequence"
        else:
            raise TypeError("Expected ReferenceFeatures or ReferenceSeqs, got %s" % type(references))

        Features(*references.get(self.selector))
        self.prefs = prefs

    def count_reads(self, library: dict):
        """Collects statistics on features assigned to each alignment associated with each read"""

        read_seq = self.sam_reader.bundle_multi_alignments(library["File"])
        self.stats.assign_library(library)

        # For each sequence in the sam file...
        for bundle, read_count in read_seq:
            bstat = self.stats.count_bundle(bundle, read_count)

            # For each alignment of the given sequence...
            for alignment in bundle:
                hits, n_candidates = self.assign_features(alignment)
                self.stats.count_bundle_assignments(bstat, alignment, hits, n_candidates)

            self.stats.finalize_bundle(bstat)

        return self.stats

    def assign_features(self, al: dict) -> Tuple[dict, int]:
        """Determines features associated with the interval then performs rule-based feature selection"""

        try:
            feat_matches = set().union(
                            *Features.chrom_vectors[al['Chrom']]['.']  # GenomicArrayOfSets -> ChromVector
                                     .array[al['Start']:al['End']]     # ChromVector -> StepVector
                                     .get_steps(values_only=True))     # StepVector -> {features}
        except KeyError as ke:
            self.stats.chrom_misses[ke.args[0]] += 1
            return {}, 0

        # If features are associated with the alignment interval, perform selection
        assignment = self.selector.choose(feat_matches, al) if feat_matches else {}
        return assignment, len(feat_matches)


class FeatureSelector:
    """Performs hierarchical selection given a set of candidate features for a locus

    Two sources of data serve as targets for selection: feature attributes (sourced from
    input GFF files), and sequence attributes (sourced from input SAM files).
    All candidate features are assumed to have matched at least one Identity selector,
    as determined by hts_parsing.ReferenceFeatures.get_matches_and_classes()

    The first round of selection was performed during GFF parsing.

    The second round is performed against the hierarchy values and
    IntervalSelectors in each candidate's match-tuples.

    If more than one candidate remains, a final round of selection is performed
    against sequence attributes: strand, 5' end nucleotide, and length. Rules for
    5' end nucleotides support lists (e.g. C,G,U) and wildcards (e.g. "all").
    Rules for length support lists, wildcards, and ranges (e.g. 20-27) which
    may be intermixed in the same rule.
    """

    rules_table: List[dict]
    inv_ident: Dict[tuple, List[int]]

    def __init__(self, rules: List[dict], **kwargs):
        FeatureSelector.rules_table = self.build_selectors(rules)
        FeatureSelector.inv_ident = self.build_inverted_identities(FeatureSelector.rules_table)
        self.warnings = defaultdict(set)
        self.overlap_cache = {}

    @classmethod
    def choose(cls, candidates: Set[feature_record_tuple], alignment: dict) -> Mapping[str, set]:
        """Selects features according to the selection rules provided at construction

        The `candidates` argument is a set of nested tuples, each representing a
        feature whose interval overlapped the alignment interval by at least one base.

            It is structured as:
                { ( feature_record_tuple 1 ), ( feature_record_tuple 2 ), ... }

            Each feature_record_tuple is structured as:
                ( featureID, strand, ( match_tuple, ... ) )

            Each match_tuple represents a rule which matched the feature on identity:
                ( rule, rank, IntervalSelector )

        Args:
            candidates: a list of tuples, each representing features associated with
                an interval which overlapped the alignment interval by at least one base.
            alignment: the alignment to which features are being selected for assignment.

        Returns:
            selections: a set of features which passed selection
        """

        # Stage 2
        hits = [(rank, rule, feat, strand)
                for feat, strand, matches in candidates
                for rule, rank, iv_match in matches
                if alignment in iv_match]

        if not hits: return {}
        hits.sort(key=lambda x: x[0])

        # Stage 3
        min_rank = None
        selections = defaultdict(set)
        for rank, rule, feat, strand in hits:
            if min_rank is not None and rank != min_rank: break

            rule_def = FeatureSelector.rules_table[rule]
            if alignment['nt5end'] not in rule_def["nt5end"]: continue
            if alignment['Length'] not in rule_def["Length"]: continue
            if (alignment['Strand'], strand) not in rule_def["Strand"]: continue

            selections[feat].add(rule)
            min_rank = rank

        return selections

    @staticmethod
    def build_selectors(rules_table) -> List[dict]:
        """Builds single/list/range/wildcard membership-matching selectors.

        Applies to: strand, 5' end nucleotide, length, and identities containing wildcards

        This function replaces text-based selector definitions in the rules_table with
        their corresponding selector classes. Selector evaluation is then performed via
        the membership operator (keyword `in`) which is handled by the selector class'
        __contains__() method.

        Precondition: rules_table preserves original row order
        """

        selector_builders = {
            "Filter_s": GffSourceMatch,
            "Filter_t": GffTypeMatch,
            "Strand": StrandMatch,
            "nt5end": NtMatch,
            "Length": NumericalMatch,
            "Identity": lambda x:x
        }

        for i, row in enumerate(rules_table):
            try:
                for selector, defn in row.items():
                    if selector not in selector_builders:
                        continue
                    if type(defn) is str and defn.lower().strip() in Wildcard.kwds:
                        row[selector] = Wildcard()
                    elif type(defn) is tuple:
                        row[selector] = tuple(Wildcard() if x.lower().strip() in Wildcard.kwds else x for x in defn)
                    else:
                        row[selector] = selector_builders[selector](defn)
            except Exception as e:
                # Append to error message while preserving exception provenance and traceback
                msg = f"Error occurred while processing rule number {i + 1}"
                append_to_exception(e, msg)
                raise

        return rules_table

    def build_interval_selectors(self, iv: 'HTSeq.GenomicInterval', match_tuples: List[unbuilt_match_tuple], feat_id=None):
        """Builds partial/wildcard/nested/exact/5'anchored/3'anchored interval selectors

        Unlike build_selectors() and build_inverted_identities(), this function
        is not called at construction time. Instead, it is called when finalizing
        match-tuples in reference parsers. This is because the intervals of features
        passing Stage 1 selection, and the specific rules they matched, must be known.

        Index 2 of each match tuple is from the Overlap column of the Features Sheet.
        It defines the desired selector and, optionally, a shift parameter for shifting
        the 5' and/or 3' ends of the interval. Its syntax is:
            selector,M,N
                M = signed shift value for 5' end
                N = signed shift value for 3' end

        Args:
            iv: The interval of the feature from which each selector is built
            match_tuples: A list of tuples representing the feature's Stage 1 matches
            feat_id: The tagged feature ID associated with these matches (for error logging)

        Returns:
            matches_by_interval: a dictionary of shifted intervals and their associated
                match tuples, where each match tuple now contains a complete IntervalSelector
        """

        cache = self.overlap_cache
        selector_factory = {
            'exact': lambda x: IntervalExactMatch(x),
            'full': lambda x: IntervalNestedMatch(x),   # temporary backward compatibility
            'nested': lambda x: IntervalNestedMatch(x),
            'partial': lambda x: IntervalPartialMatch(x),
            'anchored': lambda x: IntervalAnchorMatch(x),
            "5'anchored": lambda x: Interval5pMatch(x) if iv.strand in ('+', '-') else IntervalAnchorMatch(x),
            "3'anchored": lambda x: Interval3pMatch(x) if iv.strand in ('+', '-') else IntervalAnchorMatch(x)}

        matches_by_interval = defaultdict(list)
        for match in match_tuples:
            # Split optional shift parameters from definition
            defn = re.split(r'\s*,\s*', match[2], 1)
            selector, shift = defn[0], defn[1:]
            selector = selector.replace(' ', '')

            if selector in Wildcard.kwds:
                selector_obj = Wildcard()
                match_iv = iv
            else:
                try:
                    # Shift the interval before constructing if shift parameters were provided
                    match_iv = IntervalSelector.get_shifted_interval(shift[0], iv) if shift else iv

                    # Cache instances to prevent duplicates for the same match type on the same iv
                    cache_key = (selector, match_iv.chrom, match_iv.start, match_iv.end)
                    selector_obj = cache.setdefault(cache_key, selector_factory[selector](match_iv))
                except KeyError:
                    raise ValueError(f'Invalid overlap selector: "{match[2]}"')
                except IllegalShiftError as e:
                    bad_shift_msg = f'{feat_id} and rule {match[0]} ({match[2]})'
                    self.warnings[e.subtype].add(bad_shift_msg)
                    continue

            # Replace the match tuple's definition with selector
            built_match_tuple = (match[0], match[1], selector_obj)
            matches_by_interval[match_iv].append(built_match_tuple)

        return matches_by_interval

    @staticmethod
    def build_inverted_identities(rules_table) -> Dict[Tuple[str, str], List[int]]:
        """Builds inverted identity rules for fast matching in phase 1 selection

        The resulting dictionary has (Attrib Key, Attrib Val) as key, [associated rule indexes] as val
        """

        inverted_identities = defaultdict(list)
        for i, rule in enumerate(rules_table):
            inverted_identities[rule['Identity']].append(i)

        return dict(inverted_identities)