"""Selector classes which determine a match via __contains__()"""

import HTSeq
import re

from typing import Optional, Tuple

from tiny.rna.util import Singleton


class Wildcard(metaclass=Singleton):
    kwds = ('any', 'all', 'both', '*', '')

    @staticmethod
    def __contains__(*_): return True
    def __repr__(_): return "<all>"


class GffColumnMatch(frozenset):
    def __new__(cls, targets):
        tokenized_lowercase = map(str.strip, targets.lower().split(','))
        unique_non_empty = set(t for t in tokenized_lowercase if t)

        return super().__new__(cls, unique_non_empty)

    def __contains__(self, column_value):
        return super().__contains__(column_value.lower())


class GffSourceMatch(GffColumnMatch):
    def __repr__(self):
        return f'<Sources (col. 2) allowed: "{", ".join(self)}">'


class GffTypeMatch(GffColumnMatch):
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
    """For evaluating numerical values against a list and/or range of desired values

    FeatureSelector will determine a match by checking the frozenset
    superclass's __contains__() function. Frozensets have
    marginally faster lookup for wider ranges of values.
    """

    def __new__(cls, values):
        # Remove whitespace at start, end, "-", and ","
        ws_around_delims = r'(\s*)([,\-])(\s*)'
        cleaned = re.sub(ws_around_delims, r"\2", values.strip())
        cls.validate_definition(cleaned, values)

        # Supports intermixed lists and ranges
        rule, values = values.split(','), []
        for piece in rule:
            if '-' in piece:
                for lo, hi in re.findall(r"(\d+)-(\d+)", piece):
                    values.extend([*range(int(lo), int(hi) + 1)])
            else:
                values.append(int(piece))

        # Call frozenset superclass constructor
        return super().__new__(cls, values)

    @staticmethod
    def validate_definition(cleaned: str, values: str):
        assert re.match(r'^\d+(-\d+)?(,\d+(-\d+)?)*$', cleaned) is not None, \
            f'Invalid length selector: "{values}"'


class EditMatch(NumericalMatch):
    """Evaluates that the edit distance reported in the alignment
    dictionary's NM field is within the desired range of values."""

    def __new__(cls, values):
        return super().__new__(cls, values)

    def __contains__(self, alignment):
        return super().__contains__(alignment['NM'])


class AdarEditMatch(NumericalMatch):
    """Evaluates whether the alignment's mismatches are A -> G from ADAR enzyme editing.
    The count of A -> G mismatches must exactly equal the edit distance reported in the
    NM tag. Any other mismatches, including insertions and deletions, are disqualifying.
    """

    def __new__(cls, lengths):
        return super().__new__(cls, lengths)

    def __contains__(self, alignment):
        """Method for evaluating ADAR edit pattern.

        Args:
            alignment: a dictionary with information from the alignment's SAM record, with fields:
              - "Seq": The alignment's SEQ field, which is always 5' -> 3' regardless of strand
              - "NM": the edit distance between the aligned sequence and reference sequence.
              - "MD": the MD tag for determining reference bases
              - "Strand": True for (+) strand, False if (-) strand

        Returns: True if the A -> G mismatch count is equal to the NM value.
        """

        # Check that NM is in the desired range
        if not super().__contains__(alignment['NM']):
            return False

        nm, md = alignment['NM'], alignment['MD']
        read_seq = alignment['Seq']
        strand = alignment['Strand']

        # Tokenize the MD tag into match/mismatch runs
        md_ops = re.findall(r"\d+|\D+", md)

        # Account for SEQ aligned to reverse strand
        need_ref_nt = "A" if strand is True else "T"
        need_read_nt = "G" if strand is True else "C"

        read_pos, total_mm = 0, 0
        for run in md_ops:
            try:
                # op is a run of matches
                run = int(run)
                read_pos += run
            except ValueError:
                # op is a mismatch or run of deletions
                for reference_nt in run:
                    if reference_nt != need_ref_nt: return False
                    if read_seq[read_pos] != need_read_nt: return False

                    total_mm += 1
                    read_pos += 1

        # nm > total_mm if the read contains insertions
        return total_mm == nm


class TutEditMatch(NumericalMatch):
    """Evaluates whether the alignment's mismatches are due to 3' terminal uridylation.
    The uridylation mismatch count must exactly match the edit distance reported in the
    NM tag. Any other mismatches, including insertions and deletions, are disqualifying.
    """

    def __new__(cls, values):
        return super().__new__(cls, values)

    def __contains__(self, alignment: dict):
        """Method for evaluating TUT edit pattern.

        Args:
            alignment: a dictionary with information from the alignment's SAM record, with fields:
              - "Seq": The alignment's SEQ field, which is always 5' -> 3' regardless of strand
              - "NM": the edit distance between the aligned sequence and reference sequence.
              - "MD": the MD tag for determining reference bases
              - "Strand": True for (+) strand, False if (-) strand

        Returns: True if the 3' terminal uridylation mismatch count is equal to the NM value.
        """

        # Check that NM is in the desired range
        if not super().__contains__(alignment['NM']):
            return False

        nm, md = alignment['NM'], alignment['MD']
        strand = alignment['Strand']
        read_seq = alignment['Seq']

        # Iteration direction and read NT depend on strand
        # since SEQ, MD, and reference are always 5' -> 3'
        if strand is True:
            iter_3p_md = reversed(md)
            iter_read_md = reversed(read_seq)
            need_read_nt = "T"  # DNA
        else:
            iter_3p_md = iter(md)
            iter_read_md = read_seq
            need_read_nt = "A"  # DNA

        consecutive_3p_u = 0
        for md_3p_char, read_nt in zip(iter_3p_md, iter_read_md):
            while md_3p_char == "0":
                md_3p_char = next(iter_3p_md)        # Mismatch bases are delimited and flanked by 0

            if md_3p_char.isdigit(): break           # End of terminal run is indicated by a digit (match)
            if md_3p_char == "^": return False       # Deletions automatically disqualify the alignment
            if read_nt != need_read_nt: break        # End of terminal run is indicated by a non-uridylation nucleotide

            consecutive_3p_u += 1

        # nm == consecutive_3p_u if mismatches are only from 3' terminal uridylation
        return consecutive_3p_u == nm


# Used in IntervalSelector
class IllegalShiftError(Exception):
    def __init__(self, iv, shift, subtype):
        self.subtype = subtype
        self.shift = shift
        self.iv = iv

        self.args = (f"The interval {iv} cannot be shifted by {shift} "
                     f"(results in {self.subtype} interval)",)


class IntervalSelector:
    __slots__ = ("start", "end")

    """IntervalSelector classes use __slots__ rather than class attributes to
    reduce memory footprint. As a result, each instance requires only 128 bytes
    (compare to: 184 bytes empty class without slots, or 248 bytes empty dict).
    Instances of these classes may be numerous depending on the user's GFF file;
    one IntervalSelector is created for each unique match definition, for each
    interval associated with each retained feature."""

    def __init__(self, iv: HTSeq.GenomicInterval):
        """Descendants only need to know the start and end coordinates of the target feature"""
        self.start = iv.start
        self.end = iv.end

    def __hash__(self):
        """Descendants must be hashable in order to be stored in a GenomicArrayOfSets"""
        return self.start ^ self.end

    def __eq__(self, other):
        return self.__class__.__name__ == other.__class__.__name__ and \
               self.start == other.start and \
               self.end == other.end

    @classmethod
    def get_shifted_interval(cls, shift_defn: str, iv: HTSeq.GenomicInterval):
        """Shifts the interval's 5' and 3' ends according to the shift definition

        Positive values shift the interval in the 3' direction and negative values
        shift in the 5' direction. Both values must be specified but zero can be
        provided if no change is desired.

        Args:
            shift_defn: A string of two signed numbers, `M` and `N`, comma separated.
                `M` shifts the interval at the 5' end and `N` shifts the interval
                at the 3' end.
            iv: The interval to shift
        """

        cls.validate_shift_params(shift_defn)
        split = shift_defn.split(',', 1)
        shift = shift_5, shift_3 = int(split[0]), int(split[1])

        if iv.strand == '+':
            start, end = iv.start + shift_5, iv.end + shift_3
        elif iv.strand == '-':
            start, end = iv.start - shift_3, iv.end - shift_5
        else:  # iv.strand == '.':
            shift_x = max(shift)
            start, end = iv.start - shift_x, iv.end + shift_x
            shift = (shift_x,) * 2

        if start == end:
            raise IllegalShiftError(iv, shift, "null")
        if start > end:
            raise IllegalShiftError(iv, shift, "inverted")
        if start < 0:
            raise IllegalShiftError(iv, shift, "negative start")

        return HTSeq.GenomicInterval(iv.chrom, start, end, iv.strand)

    @staticmethod
    def validate_shift_params(defn):
        assert re.match(r'^-?\d+\s*,\s*-?\d+$', defn.strip()) is not None, \
            f'Invalid overlap shift parameters: "{defn}"'


class IntervalPartialMatch(IntervalSelector):
    """This is a no-op filter

    Since StepVector only returns features that overlap each alignment by
    at least one base, no further evaluation is required when this filter
    is used by FeatureSelector.choose()"""

    __slots__ = ()  # Base class holds attributes
    __contains__ = Wildcard.__contains__

    def __init__(self, iv: HTSeq.GenomicInterval):
        super().__init__(iv)

    def __repr__(self):
        return f"<Any overlap [{self.start}, {self.end})>"


class IntervalNestedMatch(IntervalSelector):
    __slots__ = ()  # Base class holds attributes

    def __init__(self, iv: HTSeq.GenomicInterval):
        super().__init__(iv)

    def __contains__(self, alignment):
        return self.start <= alignment['Start'] and alignment['End'] <= self.end

    def __repr__(self):
        return f"<Between [{self.start}, {self.end})>"


class IntervalExactMatch(IntervalSelector):
    __slots__ = ()  # Base class holds attributes

    def __init__(self, iv: HTSeq.GenomicInterval):
        super().__init__(iv)

    def __contains__(self, alignment):
        return self.start == alignment['Start'] and alignment['End'] == self.end

    def __repr__(self):
        return f"<Exactly [{self.start}, {self.end})>"


class IntervalAnchorMatch(IntervalSelector):
    """Evaluates whether either end of the alignment's interval is anchored to
    the feature with the non-anchored end nested within the feature's interval"""

    __slots__ = ()  # Base class holds attributes

    def __init__(self, iv: HTSeq.GenomicInterval):
        super().__init__(iv)

    def __contains__(self, alignment: dict):
        """The following diagram demonstrates unstranded anchored matching semantics.

                    Match  |   <------->|
                    Match  |<------->   |
                    Match  |<---------->|
                 No match  |  <---------|->
               <-----------|<==feat_A==>|----------->

        Args:
            alignment: An alignment dictionary containing `Start` and `End` keys
        """

        return (alignment['Start'] == self.start and alignment['End'] <= self.end) \
            or (alignment['End'] == self.end and alignment['Start'] >= self.start)

    def __repr__(self):
        return f"<Unstranded anchor from {self.start} to {self.end}>"


class Interval5pMatch(IntervalSelector):
    """Evaluates whether an alignment's 5' end is anchored to the corresponding terminus
    of the feature, and the alignment's 3' end is nested within the feature's interval."""

    __slots__ = ()  # Base class holds attributes

    def __init__(self, iv: HTSeq.GenomicInterval):
        super().__init__(iv)

    def __contains__(self, alignment: dict):
        """The following diagram demonstrates 5' anchored matching semantics.
        Each "(no) match" label applies to BOTH feat_A AND feat_B.

                 No match  |   -------->|
                 No Match  |------------|-->
                    Match  |----------->|
                    Match  |-------->   |
        (+) 5' ------------|==feat_A===>|-----------> 3'

        (-) 3' <-----------|<===feat_B==|------------ 5'
                           |   <--------|  Match
                           |<-----------|  Match
                        <--|------------|  No Match
                           |<--------   |  No match
        Args:
            alignment: An alignment dictionary containing `Strand`, `Start` and `End` keys
        """

        if alignment['Strand'] is True:
            return alignment['Start'] == self.start and alignment['End'] <= self.end
        else:
            return alignment['End'] == self.end and alignment['Start'] >= self.start

    def __repr__(self):
        return f"<5' anchored to {self.start} nested in {self.end} on (+) / " \
               f"anchored to {self.end} nested in {self.start} on (-)>"


class Interval3pMatch(IntervalSelector):
    """Evaluates whether an alignment's 3' end is anchored to the corresponding terminus
    of the feature, and the alignment's 5' end is nested within the feature's interval."""

    __slots__ = ()  # Base class holds attributes

    def __init__(self, iv: HTSeq.GenomicInterval):
        super().__init__(iv)

    def __contains__(self, alignment):
        """The following diagram demonstrates 3' anchored matching semantics.
        Each "(no) match" label applies to BOTH feat_A AND feat_B.

               No match    |-------->   |
               No Match  --|----------->|
                  Match    |----------->|
                  Match    |   -------->|
        (+) 5' ------------|==feat_A===>|-----------> 3'

        (-) 3' <-----------|<===feat_B==|------------ 5'
                           |<--------   |    Match
                           |<-----------|    Match
                           |<-----------|--  No Match
                           |   <--------|    No match
        Args:
            alignment: An alignment dictionary containing `Strand`, `Start` and `End` keys
        """

        if alignment['Strand'] is True:
            return alignment['End'] == self.end and alignment['Start'] >= self.start
        else:
            return alignment['Start'] == self.start and alignment['End'] <= self.end

    def __repr__(self):
        return f"<3' anchored to {self.end} nested in {self.start} on (+) / " \
               f"anchored to {self.start} nested in {self.end} on (-)>"
