"""
Collapse sequences from a fastq file to a fasta file. Headers in the output fasta file
will contain the number of times each sequence occurred in the input fastq file, and
an ID which indicates the relative order in which each sequence was first encountered.
Gzipped files are automatically supported for fastq inputs, and compressed fasta outputs
are available by request.
"""

import argparse
import builtins
import gzip
import os
import re

from collections import Counter
from functools import partial
from typing import Tuple, Iterable

try:
    from _collections import _count_elements  # Load Counter's C helper function if it is available
except ImportError:
    from collections import _count_elements   # Slower mapping[elem] = mapping.get(elem,default_val)+1

# The GZIP read/write interface used by seq_counter() and seq2fasta()
gz_f = partial(gzip.GzipFile, compresslevel=6, fileobj=None, mtime=0)


def get_args() -> 'argparse.NameSpace':
    """Get command line arguments"""

    parser = argparse.ArgumentParser(description=__doc__)
    required_group = parser.add_argument_group("required arguments")

    # Required arguments
    required_group.add_argument(
        '-i', '--input-file', metavar='FASTQFILE', required=True, help=
        'The input fastq(.gz) file to collapse'
    )

    required_group.add_argument(
        '-o', '--out-prefix', metavar='OUTPREFIX', required=True, help=
        'The prefix for output files {prefix}_collapsed.fa and, if '
        'counts fall below threshold, {prefix}_collapsed_lowcounts.fa'
    )

    def positive_number(t):
        if int(t) >= 0:
            return int(t)
        else:
            raise argparse.ArgumentTypeError("Numerical arguments must be >= 0")

    # Optional arguments
    parser.add_argument(
        '-t', '--threshold', default=0, required=False, type=positive_number,
        help='Sequences <= THRESHOLD will be omitted from {prefix}_collapsed.fa '
        'and will instead be placed in {prefix}_collapsed_lowcounts.fa'
    )

    parser.add_argument(
        '-c', '--compress', required=False, action='store_true',
        help='Use gzip compression when writing fasta outputs'
    )

    parser.add_argument(
        '--3p-trim', metavar='N BASES', required=False, type=positive_number,
        help="Treat N BASES at the 3' end of each sequence as a UMI and trim them"
    )

    parser.add_argument(
        '--5p-trim', metavar='N BASES', required=False, type=positive_number,
        help="Treat N BASES at the 5' end of each sequence as a UMI and trim them"
    )

    parser.add_argument(
        '-d', '--discard_umi_duplicates', action='store_true',
        help="Discard duplicate UMIs. The --3p-trim and/or --5p-trim arguments "
             "define which part of each sequence to treat as the UMI sequence."
    )

    args = parser.parse_args()
    umi = {
        'dedup': args.discard_umi_duplicates,
        'trim5': getattr(args, '5p_trim', 0),
        'trim3': getattr(args, '3p_trim', 0)
    }

    if umi['dedup'] and not (umi['trim3'] or umi['trim5']):
        raise argparse.ArgumentError(None, "UMI deduplication requires the length and location of the UMI(s) to be "
                                           "specified with --3p-trim and/or --5p-trim.")

    return args, umi

def seq_counter(fastq_file: str, file_reader: callable = builtins.open, umi: dict = None) -> 'Counter':
    """Counts the number of times each sequence appears

    Args:
        fastq_file: A trimmed, quality filtered, optionally gzip compressed fastq file.
        file_reader: The file context manager to use. Must support .readline() and 'rb'
        umi: A dictionary specifying 5'/3' trim lengths and deduplication preferences

    Returns: An ordered dictionary of unique sequences with associated counts.
    """

    with file_reader(fastq_file, 'rb') as f:
        def sequences():    # Generator function for every 4th line (fastq sequence line) of file
            while f.readline():  # Sequence identifier
                # Sequence (Binary -> ASCII extract every 4th from 1st line, newline removed)
                yield f.readline()[:-1].decode("utf-8")
                f.readline()     # "+"
                f.readline()     # Quality Score

        # Switch file_reader interface if reading gzipped fastq files
        if f.read(2) == b'\x1F\x8B': return seq_counter(fastq_file, gz_f, umi)

        # Count occurrences of unique sequences while maintaining insertion order
        counts = Counter(sequences() if umi is None else with_umi(sequences(), umi))

    counts.pop("", None)  # Remove blank line counts from the dictionary
    return counts


def with_umi(seq_iter: Iterable, umi_spec: dict) -> Iterable:
    """Restructure each FASTQ sequence for counting in accordance with the umi_spec

    If deduplication is requested, the resulting sequence

    Args:
        seq_iter: iterator of sequences read from FASTQ file
        umi_spec: a dictionary specifying 5'/3' trim lengths and a boolean flag
            to indicate whether deduplication should be performed

    Returns: nested generator yielding the appropriate UMI tuple
    """

    trim5 = umi_spec['trim5']
    trim3 = umi_spec['trim3']
    dedup = umi_spec['dedup']

    # I know this looks terrible, but hear me out:
    # The decision tree is evaluated just once then
    # settles into/resumes the appropriate loop.
    if dedup:
        if trim5 and trim3:
            for seq in seq_iter:
                yield seq[:trim5], seq[trim5:-trim3], seq[-trim3:]
        elif trim5:
            for seq in seq_iter:
                yield seq[:trim5], seq[trim5:]
        elif trim3:
            for seq in seq_iter:
                yield seq[:-trim3], seq[-trim3:]
    else:
        if trim5 and trim3:
            for seq in seq_iter:
                yield seq[trim5:-trim3]
        elif trim5:
            for seq in seq_iter:
                yield seq[trim5:]
        elif trim3:
            for seq in seq_iter:
                yield seq[:-trim3]


def with_umi_regex(seq_iter: Iterable, umi_spec: dict) -> Iterable:
    """Marginally more readable (and marginally slower) version for testing"""

    trim5 = umi_spec['trim5']
    trim3 = umi_spec['trim3']
    dedup = umi_spec['dedup']

    def umigrp(length):
        if length == 0: return ""
        pattern = r"(.{" + length + "})"
        if not dedup: pattern = pattern[1:-1]
        return pattern

    umipat = umigrp(trim5) + r"(.*)" + umigrp(trim3)
    umipat = re.compile(umipat)

    for seq in seq_iter:
        yield umipat.findall(seq)[0]


def seq2fasta(seqs: dict, out_prefix: str, thresh: int = 0, gz: bool = False) -> None:
    """Converts a sequence count dictionary to a fasta file, with count filtering

    If a threshold is specified, sequences with count > thresh will be written to
    {out_prefix}_collapsed.fa, and sequences with count <= thresh will be written
    to {out_prefix}_collapsed_lowcounts.fa. If the specified threshold results in
    an empty collection for either output file, a blank file will still be created
    under the corresponding name.

        The fasta header is formatted as:
            >ID_count=COUNT

    Headers indicate the ID of the sequence and the sequence count. The first
    unique sequence is assigned ID 0, and the second unique sequence, (ID 1),
    may have n repetitions of sequence ID 0 before it in the fastq file, but its
    ID will be 1 not n+1.

    Args:
        seqs: A dictionary containing sequences and associated counts
        out_prefix: A prefix name for the output fasta files
        thresh: Sequences with count LE thresh will placed in a separate file
        gz: If true, fasta outputs will be gzip compressed

    Returns: None
    """

    assert out_prefix is not None, "Collapser critical error: an output file prefix must be specified."
    assert thresh >= 0, "An invalid threshold was specified."

    writer, encoder, mode = fasta_interface(gz)
    out_file, low_count_file = look_before_you_leap(out_prefix, gz)

    def to_fasta_record(x):
        # x[0]=ID, x[1][1]=sequence count, x[1][0]=sequence
        return ">%d_count=%d\n%s" % (x[0], x[1][1], x[1][0])

    above_thresh = filter(lambda x: x[1][1] > thresh, enumerate(seqs.items()))
    below_thresh = filter(lambda x: x[1][1] <= thresh, enumerate(seqs.items()))

    with writer(out_file, mode) as fasta:

        if thresh == 0:  # No filtering required
            fasta.write(encoder('\n'.join(map(to_fasta_record, enumerate(seqs.items())))))
        else:
            with writer(low_count_file, mode) as lowfa:
                fasta.write(encoder('\n'.join(map(to_fasta_record, above_thresh))))
                lowfa.write(encoder('\n'.join(map(to_fasta_record, below_thresh))))


def look_before_you_leap(out_prefix: str, gz: bool) -> (str, str):
    """Check that we'll be able to write results before we spend time on the work"""

    ext = '.fa.gz' if gz else '.fa'
    candidates = tuple(f"{out_prefix}{file}{ext}" for file in ["_collapsed", "_collapsed_lowcounts"])
    for file in candidates:
        if os.path.isfile(file):
            raise FileExistsError(f"Collapser critical error: {file} already exists.")

    return candidates


def fasta_interface(gz: bool) -> Tuple[callable, callable, str]:
    """Switches to writing via gzip.GzipFile() if fasta compression is specified"""

    if gz:
        # Writing gzip requires byte array input
        def encoder(x): return x.encode('utf-8')
        writer, mode = gz_f, 'wb'
    else:
        # No gzip, no conversion
        def encoder(x): return x
        writer, mode = open, 'w'

    return writer, encoder, mode


def main():
    # Get command line arguments
    args, umi_spec = get_args()
    # Ensure that the provided prefix will not result in overwritten output files
    look_before_you_leap(args.out_prefix, args.compress)
    # Count unique sequences in input fastq file
    seqs = seq_counter(args.input_file)
    # Write counted sequences to output file(s)
    seq2fasta(seqs, args.out_prefix, args.threshold, args.compress)


if __name__ == '__main__':
    main()
