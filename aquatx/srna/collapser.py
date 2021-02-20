"""
Collapse sequences from a fastq file to a fasta file. Headers in the output fasta file
will contain the number of times each sequence occurred in the input fastq file, and
an ID which indicates the relative order in which each sequence was first encountered.
"""

import argparse
import os

from collections import OrderedDict

try:
    # Load Counter's C helper function if it is available. Uses 30% less memory than Counter()
    from _collections import _count_elements
except ImportError:
    # Uses slower mapping[elem] = mapping.get(elem,default_val)+1 in loop
    from collections import _count_elements


def get_args() -> 'argparse.NameSpace':
    """Get command line arguments"""

    def positive_threshold(t):
        if int(t) >= 0:
            return int(t)
        else:
            raise argparse.ArgumentTypeError("Threshold must be >= 0")

    parser = argparse.ArgumentParser(description=__doc__)
    required_group = parser.add_argument_group("required arguments")

    # Required arguments
    required_group.add_argument(
        '-i', '--input-file', metavar='FASTQFILE', required=True, help=
        'The input fastq file to collapse'
    )

    required_group.add_argument(
        '-o', '--out-prefix', metavar='OUTPREFIX', required=True, help=
        'The prefix for output files {prefix}_collapsed.fa and, if '
        'counts fall below threshold, {prefix}_collapsed_lowcounts.fa'
    )

    # Optional arguments
    parser.add_argument(
        '-t', '--threshold', default=0, required=False, type=positive_threshold,
        help='Sequences <= threshold will be omitted from {prefix}_collapsed.fa '
        'and will instead be placed in {prefix}_collapsed_lowcounts.fa'
    )

    return parser.parse_args()


def seq_counter(fastq_file: str) -> 'OrderedDict':
    """Counts the number of times each sequence appears

    Args:
        fastq_file: A trimmed, quality filtered fastq file
                    containing sequences to count.

    Returns: An ordered dictionary of unique sequences with associated counts.
    """

    with open(fastq_file, 'rb') as f:
        def line_generator():    # Generator function for every 4th line (fastq sequence line) of file
            while f.readline():  # Sequence identifier
                # Sequence (Binary -> ASCII extract every 4th from 1st line, newline removed)
                yield f.readline()[:-1].decode("utf-8")
                f.readline()     # "+"
                f.readline()     # Quality Score

        # Count occurrences of unique sequences while maintaining insertion order
        seqs = OrderedDict()
        _count_elements(seqs, line_generator())

    seqs.pop("", None)  # Remove blank line counts from the dictionary
    return seqs


def seq2fasta(seqs: dict, out_prefix: str, thresh: int = 0) -> None:
    """Converts a sequence count dictionary to a fasta file, with count filtering

    If a threshold is specified, sequences with count > thresh will be written to
    {out_prefix}_collapsed.fa, and sequences with count <= thresh will be written
    to {out_prefix}_collapsed_lowcounts.fa. If the specified threshold results in
    an empty collection for either output file, a blank file will still be created
    under the corresponding name.

        The fasta header is formatted as:
            >ID_count=COUNT

    Headers indicate the ID of the sequence and the sequence count. The first
    unique sequence is assigned ID 0, and the second unique sequence (ID 1)
    may have n repetitions of sequence ID 0 before it in the fastq file, but its
    ID will be 1, not n+1.

    Args:
        seqs: A dictionary containing sequences and associated counts
        out_prefix: A prefix name for the output fasta files
        thresh: Sequences with count <= thresh will placed in a separate file

    Returns: None
    """

    assert out_prefix is not None, "Collapser critical error: an output file prefix must be specified."
    assert thresh >= 0, "An invalid threshold was specified."

    def to_fasta_record(x):
        # x[0]=ID, x[1][1]=sequence count, x[1][0]=sequence
        return ">%d_count=%d\n%s\n" % (x[0], x[1][1], x[1][0])

    out_file, low_count_file = look_before_you_leap(out_prefix)
    above_thresh = filter(lambda x: x[1][1] > thresh, enumerate(seqs.items()))
    below_thresh = filter(lambda x: x[1][1] <= thresh, enumerate(seqs.items()))

    with open(out_file, 'w') as fasta:
        if thresh == 0:  # No filtering required
            fasta.write(map(to_fasta_record, enumerate(seqs.items())))
        else:
            with open(low_count_file, 'w') as lowfa:
                fasta.write(map(to_fasta_record, above_thresh))
                lowfa.write(map(to_fasta_record, below_thresh))


def look_before_you_leap(out_prefix: str) -> (str, str):
    """Check that we'll be able to write results before we spend time on the work"""

    out_file, low_count_file = f"{out_prefix}_collapsed.fa", f"{out_prefix}_collapsed_lowcounts.fa"
    for file in [out_file, low_count_file]:
        if os.path.isfile(file):
            raise FileExistsError(f"Collapser critical error: {file} already exists.")

    return out_file, low_count_file


def main():
    # Get command line arguments
    args = get_args()
    # Ensure that the provided prefix will not result in overwritten output files
    look_before_you_leap(args.out_prefix)
    # Count unique sequences in input fastq file
    seqs = seq_counter(args.input_file)
    # Write counted sequences to output file(s)
    seq2fasta(seqs, args.out_prefix, args.threshold)


if __name__ == '__main__':
    main()
