"""Selector classes which determine a match via __contains__()"""

import HTSeq
import re

from typing import Optional, Tuple

from tiny.rna.util import Singleton


class Wildcard(metaclass=Singleton):
    kwds = ('all', 'both', '*', '')

    @staticmethod
    def __contains__(*_): return True
    def __repr__(_): return "<all>"


class GffColumnMatch(frozenset):
    def __new__(cls, targets):
        tokenized_lowercase = map(str.strip, targets.lower().split(','))
        unique_non_empty = set(t for t in tokenized_lowercase if t)

        return super().__new__(cls, unique_non_empty)


class GffSourceMatch(GffColumnMatch):

    # def __new__(cls, sources):
    #     super().__new__(cls, sources)
    
    def __repr__(self):
        return f'<Sources (col. 2) allowed: "{", ".join(self)}">'


class GffTypeMatch(GffColumnMatch):

    # def __new__(cls, sources):
    #     super().__new__(cls, sources)

    def __repr__(self):
        return f'<Types (col. 3) allowed: "{", ".join(self)}">'


class StrandMatch:
    """Evaluates BOTH the alignment's strand and the feature's strand for a match

    If sense: alignment's strand == feature strand for a match
    If antisense: alignment's strand != feature strand for a match
    """
    def __init__(self, strand):
        strand = strand.lower().strip()
        self.validate_definition(strand)

        self.strand = strand
        self.select = (self.strand == 'sense')

    def __contains__(self, strand_relationship: Tuple[bool, Optional[bool]]):
        """If feature strand is None, a match is not possible for 'sense' or 'antisense'"""

        aln_strand, feat_strand = strand_relationship
        return feat_strand is None or aln_strand ^ feat_strand ^ self.select

    def __repr__(self): return str(self.strand)

    @staticmethod
    def validate_definition(defn: str):
        assert defn in ("sense", "antisense"), f'Invalid strand selector: "{defn}"'


class NtMatch(tuple):
    """For evaluating a single nucleotide against a list of desired bases

    FeatureSelector will determine a match by checking the tuple
    superclass's __contains__() function via the `in` keyword.
    Tuples are more lightweight for a small range of possible values.
    """

    def __new__(cls, nts):
        # Supports single or comma separated values
        nts = map(lambda x: x.strip().upper(), nts.split(','))

        unique_nts = set(nts)
        cls.validate_definition(unique_nts)

        # Call tuple superclass constructor
        return super().__new__(cls, unique_nts)

    @staticmethod
    def validate_definition(unique_nts: set):
        non_nt = unique_nts - {'A', 'T', 'G', 'C', 'N'}
        assert len(non_nt) == 0, f"""Invalid nucleotide selector: {
            ", ".join(f'"{x}"' for x in non_nt)}"""


class NumericalMatch(frozenset):
    """For evaluating sequence length against a list and/or range of desired values

    FeatureSelector will determine a match by checking the frozenset
    superclass's __contains__() function. Frozensets have
    marginally faster lookup for wider ranges of values.
    """

    def __new__(cls, lengths):
        cls.validate_definition(lengths)

        # Supports intermixed lists and ranges
        rule, lengths = lengths.split(','), []
        for piece in rule:
            if '-' in piece:
                for lo, hi in re.findall(r"(\d+)-(\d+)", piece):
                    lengths.extend([*range(int(lo), int(hi) + 1)])
            else:
                lengths.append(int(piece))

        # Call frozenset superclass constructor
        return super().__new__(cls, lengths)

    @staticmethod
    def validate_definition(defn: str):
        assert re.match(r'(^\d+$)|^(\d+)([\d\-,]|(, +))*(\d|\1)$', defn) is not None, \
            f'Invalid length selector: "{defn}"'


class IntervalSelector:
    __slots__ = ("start", "end")

    """IntervalSelector classes use __slots__ rather than class attributes to
    reduce memory footprint. As a result, each instance requires only 128 bytes
    (compare to: 184 bytes empty class without slots, or 248 bytes empty dict).
    Instances of these classes may be numerous depending on the user's GFF file;
    one IntervalSelector is created for each unique match definition, for each
    interval associated with each retained feature."""

    def __init__(self, iv: HTSeq.GenomicInterval):
        """Descendents only need to know the start and end coordinates of the target feature"""
        self.start = iv.start
        self.end = iv.end

    def __hash__(self):
        """Descendents must be hashable in order to be stored in a GenomicArrayOfSets"""
        return self.start ^ self.end

    def __eq__(self, other):
        return self.__class__.__name__ == other.__class__.__name__ and \
               self.start == other.start and \
               self.end == other.end


class IntervalPartialMatch(IntervalSelector):
    """This is a no-op filter

    Since StepVector only returns features that overlap each alignment by
    at least one base, no further evaluation is required when this filter
    is used by FeatureSelector.choose()"""

    __contains__ = Wildcard.__contains__

    def __init__(self, iv: HTSeq.GenomicInterval):
        super().__init__(iv)

    def __repr__(self):
        return f"<Any overlap [{self.start}, {self.end})>"


class IntervalFullMatch(IntervalSelector):

    def __init__(self, iv: HTSeq.GenomicInterval):
        super().__init__(iv)

    def __contains__(self, alignment):
        return self.start <= alignment['Start'] and alignment['End'] <= self.end

    def __repr__(self):
        return f"<Between [{self.start}, {self.end})>"


class IntervalExactMatch(IntervalSelector):

    def __init__(self, iv: HTSeq.GenomicInterval):
        super().__init__(iv)

    def __contains__(self, alignment):
        return self.start == alignment['Start'] and alignment['End'] == self.end

    def __repr__(self):
        return f"<Exactly [{self.start}, {self.end})>"


class IntervalAnchorMatch(IntervalSelector):
    """Evaluates whether an alignment's start matches the feature's start, and vice versa for end."""

    def __init__(self, iv: HTSeq.GenomicInterval):
        super().__init__(iv)
        assert iv.strand not in ('+', '-')

    def __contains__(self, alignment: dict):
        """The following diagram demonstrates unstranded anchored matching semantics.

                    Match  |   <------->|
                    Match  |<-----------|-->
                    Match  |<---------->|
                 No match  |  <---------|->
               <-----------|<==feat_A==>|----------->

        Args:
            alignment: An alignment dictionary containing strand, start, and end

        Returns: True if the alignment's 5' end is anchored to the strand-appropriate
            terminus of this feature's interval.
        """

        return alignment['Start'] == self.start or alignment['End'] == self.end

    def __repr__(self):
        return f"<Unstranded anchor from {self.start} to {self.end}>"

class Interval5pMatch(IntervalSelector):
    """Evaluates whether an alignment's 5' end is anchored to the corresponding terminus of the feature"""

    def __init__(self, iv: HTSeq.GenomicInterval):
        super().__init__(iv)

    def __contains__(self, alignment: dict):
        """The following diagram demonstrates 5' anchored matching semantics.
        Each "(no) match" label applies to BOTH feat_A AND feat_B.

                 No match  |   -------->|
                    Match  |------------|-->
                    Match  |----------->|
                    Match  |-------->   |
        (+) 5' ------------|==feat_A===>|-----------> 3'

        (-) 3' <-----------|<===feat_B==|------------ 5'
                           |   <--------|  Match
                           |<-----------|  Match
                        <--|------------|  Match
                           |<--------   |  No match
        Args:
            alignment: An alignment dictionary containing strand, start, and end

        Returns: True if the alignment's 5' end is anchored to the strand-appropriate
            terminus of this feature's interval.
        """

        if alignment['Strand'] is True:
            return alignment['Start'] == self.start
        else:
            return alignment['End'] == self.end

    def __repr__(self):
        return f"<5' anchored to {self.start} on (+) / {self.end} on (-)>"


class Interval3pMatch(IntervalSelector):
    """Evaluates whether an alignment's 5' end is anchored to the corresponding terminus of the feature"""

    def __init__(self, iv: HTSeq.GenomicInterval):
        super().__init__(iv)

    def __contains__(self, alignment):
        """The following diagram demonstrates 3' anchored matching semantics.
        Each "(no) match" label applies to BOTH feat_A AND feat_B.

               No match    |-------->   |
                  Match  --|----------->|
                  Match    |----------->|
                  Match    |   -------->|
        (+) 5' ------------|==feat_A===>|-----------> 3'

        (-) 3' <-----------|<===feat_B==|------------ 5'
                           |<--------   |    Match
                           |<-----------|    Match
                           |<-----------|--  Match
                           |   <--------|    No match
        Args:
            alignment: An alignment dictionary containing strand, start, and end

        Returns: True if the alignment's 3' end is anchored to the strand-appropriate
            terminus of this feature's interval.
        """

        if alignment['Strand'] is True:
            return alignment['End'] == self.end
        else:
            return alignment['Start'] == self.start

    def __repr__(self):
        return f"<3' anchored to {self.end} on (+) / {self.start} on (-)>"
