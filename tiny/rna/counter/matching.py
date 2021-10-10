"""Selector classes which determine a match via __contains__()"""

import re


class Wildcard:
    @staticmethod
    def __contains__(x): return True
    def __repr__(self): return "all"


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
    superclass's __contains__() function. Frozensets have marginally faster lookup for wider ranges of values.
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