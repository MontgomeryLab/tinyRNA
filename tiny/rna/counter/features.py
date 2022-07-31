import HTSeq
import sys

from collections import defaultdict
from typing import List, Tuple, Set, Dict, Mapping

from tiny.rna.counter.hts_parsing import ReferenceTables, SAM_reader
from .statistics import LibraryStats
from .matching import *

# Type aliases for human readability
match_tuple = Tuple[int, int, IntervalSelector]             # (rank, rule, interval selector)
feature_record_tuple = Tuple[str, str, Tuple[match_tuple]]  # (feature ID, strand, match tuple)


class Features(metaclass=Singleton):
    chrom_vectors: HTSeq.ChromVector
    aliases: dict
    classes: dict
    tags: dict

    def __init__(_, features: HTSeq.GenomicArrayOfSets, aliases: dict, classes: dict, tags: dict):
        Features.chrom_vectors = features.chrom_vectors  # For interval -> feature record tuple lookups
        Features.aliases = aliases                       # For feature ID -> preferred feature name lookups
        Features.classes = classes                       # For feature ID -> class lookups
        Features.tags = tags                             # For feature ID -> match IDs


class FeatureCounter:

    def __init__(self, gff_file_set, selection_rules, **prefs):
        self.stats = LibraryStats(**prefs)
        self.sam_reader = SAM_reader(**prefs)
        self.selector = FeatureSelector(selection_rules, self.stats, **prefs)

        reference_tables = ReferenceTables(gff_file_set, self.selector, **prefs)
        Features(*reference_tables.get())
        self.prefs = prefs

    def assign_features(self, al: dict) -> Tuple[dict, int]:
        """Determines features associated with the interval then performs rule-based feature selection"""

        feat_matches, assignment = set(), {}

        try:
            feat_matches = feat_matches.union(*(match for match in
                    (Features.chrom_vectors[al['chrom']]['.']       # GenomicArrayOfSets -> ChromVector
                             .array[al['start']:al['end']]          # ChromVector -> StepVector
                             .get_steps(values_only=True))          # StepVector -> {features}
                    # If an alignment does not map to a feature, an empty set is returned
                    if len(match) != 0))
        except KeyError as ke:
            self.stats.chrom_misses[ke.args[0]] += 1

        # If features are associated with the alignment interval, perform selection
        if len(feat_matches):
            assignment = self.selector.choose(feat_matches, al)

        return assignment, len(feat_matches)

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


class FeatureSelector:
    """Performs hierarchical selection given a set of candidate features for a locus

    Two sources of data serve as targets for selection: feature attributes (sourced from
    input GFF files), and sequence attributes (sourced from input SAM files).
    All candidate features are assumed to have matched at least one Identity selector,
    as determined by hts_parsing.ReferenceTables.get_matches_and_classes()

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

    def __init__(self, rules: List[dict], libstats: 'LibraryStats', report_diags=False, **kwargs):
        FeatureSelector.rules_table = self.build_selectors(rules)
        FeatureSelector.inv_ident = self.build_inverted_identities(FeatureSelector.rules_table)

        self.report_eliminations = report_diags
        if report_diags: self.elim_stats = libstats.diags.selection_diags

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

        hits, min_rank = [], sys.maxsize

        for feat, strand, matches in candidates:
            for rule, rank, iv_match in matches:
                if rank > min_rank: continue
                if alignment not in iv_match: continue
                if rank < min_rank: min_rank = rank
                hits.append((rank, rule, feat, strand))

        if not hits: return {}

        selections = defaultdict(set)
        for hit in hits:
            if hit[0] != min_rank: continue
            _, rule, feat, strand = hit

            strand = (alignment['strand'], strand)
            nt5end = alignment['nt5']
            length = alignment['len']

            rule_def = FeatureSelector.rules_table[rule]
            if strand not in rule_def["Strand"]: continue
            if nt5end not in rule_def["nt5end"]: continue
            if length not in rule_def["Length"]: continue
            selections[feat].add(rule)

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

        selector_builders = {"Strand": StrandMatch, "nt5end": NtMatch, "Length": NumericalMatch, "Identity": lambda x:x}

        for i, row in enumerate(rules_table):
            try:
                for selector, build_fn in selector_builders.items():
                    defn = row[selector]

                    if type(defn) is str and defn.lower().strip() in Wildcard.kwds:
                        row[selector] = Wildcard()
                    elif type(defn) is tuple:
                        row[selector] = tuple(Wildcard() if x.lower().strip() in Wildcard.kwds else x for x in defn)
                    else:
                        row[selector] = build_fn(defn)
            except Exception as e:
                # Append to error message while preserving exception provenance and traceback
                e.args = (str(e.args[0]) + '\n' + f"Error occurred while processing rule number {i + 2}",)
                raise e.with_traceback(sys.exc_info()[2]) from e

        return rules_table

    @staticmethod
    def build_interval_selectors(iv: 'HTSeq.GenomicInterval', match_tuples: List[Tuple]):
        """Builds partial/full/exact/3' anchored/5' anchored interval selectors

        Unlike build_selectors() and build_inverted_identities(), this function
        is not called at construction time. Instead, it is called when finalizing
        match-tuples in ReferenceTables. This is because interval selectors are
        created for each feature (requiring start/stop/strand to be known) for
        each of the feature's identity matches (each match-tuple).

        Args:
            iv: The interval of the feature from which each selector is built
            match_tuples: A list of tuples representing the feature's identity
                matches. Each tuple index 2 defines and is replaced by the selector.
        """

        built_selectors = {}
        selector_factory = {
            'full': lambda: IntervalFullMatch(iv),
            'exact': lambda: IntervalExactMatch(iv),
            'partial': lambda: IntervalPartialMatch(iv),
            "5' anchored": lambda: Interval5pMatch(iv),
            "3' anchored": lambda: Interval3pMatch(iv),
        }

        for i in range(len(match_tuples)):
            try:
                match = match_tuples[i]
                # Cache instances to prevent duplicates for the same match type on the same iv
                selector = built_selectors.setdefault(match[2], selector_factory[match[2]]())
                match_tuples[i] = (match[0], match[1], selector)
            except KeyError:
                raise ValueError(f'Invalid overlap selector: "{match_tuples[i][2]}"')

        return match_tuples

    @staticmethod
    def build_inverted_identities(rules_table) -> Dict[Tuple[str, str], List[int]]:
        """Builds inverted identity rules for fast matching in phase 1 selection

        The resulting dictionary has (Attrib Key, Attrib Val) as key, [associated rule indexes] as val
        """

        inverted_identities = defaultdict(list)
        for i, rule in enumerate(rules_table):
            inverted_identities[rule['Identity']].append(i)

        return dict(inverted_identities)