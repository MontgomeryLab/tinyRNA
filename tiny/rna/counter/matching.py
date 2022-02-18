"""Selector classes which determine a match via __contains__()"""

import re
import HTSeq


class Wildcard:
    @staticmethod
    def __contains__(_): return True
    def __repr__(_): return "all"


class StrandMatch:
    """Evaluates BOTH the alignment's strand and the feature's strand for a match

    If sense: alignment's strand == feature strand for a match
    If antisense: alignment's strand != feature strand for a match
    """
    def __init__(self, strand):
        self.strand = strand.lower()
        self.select = (self.strand == 'sense')

    def __contains__(self, x):
        return self.select == (x[0] == x[1])

    def __repr__(self): return str(self.strand)


class NtMatch(tuple):
    """For evaluating a single nucleotide against a list of desired bases

    FeatureSelector will determine a match by checking the tuple
    superclass's __contains__() function via the `in` keyword.
    Tuples are more lightweight for a small range of possible values.
    """

    def __new__(cls, nts):
        if type(nts) is not tuple:
            # Supports single or comma separated values
            nts = map(lambda x: x.strip().upper(), nts.split(','))

        # Call tuple superclass constructor
        return super().__new__(cls, nts)


class NumericalMatch(frozenset):
    """For evaluating sequence length against a list and/or range of desired values

    FeatureSelector will determine a match by checking the frozenset
    superclass's __contains__() function. Frozensets have
    marginally faster lookup for wider ranges of values.
    """

    def __new__(cls, lengths):
        if type(lengths) is not list:
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


class IntervalPartialMatch(Wildcard):
    """This is a no-op filter

    Since StepVector only returns features that overlap each alignment by
    at least one base, no further evaluation is required when this filter
    is used by FeatureSelector.choose()"""

    __slots__ = ("start", "end")

    def __init__(self, iv: HTSeq.GenomicInterval):
        self.start = iv.start
        self.end = iv.end

    def __repr__(self):
        return f"Any overlap {self.start}:{self.end}"


class IntervalFullMatch:
    __slots__ = ("start", "end")

    def __init__(self, iv: HTSeq.GenomicInterval):
        self.start = iv.start
        self.end = iv.end

    def __contains__(self, alignment):
        return self.start <= alignment['start'] and alignment['end'] <= self.end

    def __repr__(self):
        return f"Between {self.start}:{self.end}"


class IntervalExactMatch:
    __slots__ = ("start", "end")

    def __init__(self, iv: HTSeq.GenomicInterval):
        self.start = iv.start
        self.end = iv.end

    def __contains__(self, alignment):
        return self.start == alignment['start'] and alignment['end'] == self.end

    def __repr__(self):
        return f"Exactly {self.start}:{self.end}"


class Interval5pMatch:
    __slots__ = ("start", "strand")

    def __init__(self, iv: HTSeq.GenomicInterval):
        self.strand = iv.strand
        self.start = iv.start

    def __contains__(self, alignment):
        if self.strand == alignment['strand']:
            return self.start == alignment['start']
        else:
            return self.start == alignment['end']

    def __repr__(self):
        return f"5' anchored to {self.start} ({self.strand})"


class Interval3pMatch:
    __slots__ = ("end", "strand")

    def __init__(self, iv: HTSeq.GenomicInterval):
        self.strand = iv.strand
        self.end = iv.end

    def __contains__(self, alignment):
        if self.strand == alignment['strand']:
            return self.end == alignment['end']
        else:
            return self.end == alignment['start']

    def __repr__(self):
        return f"3' anchored to {self.end} ({self.strand})"