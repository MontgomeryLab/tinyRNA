""" Collapse sequences from a fastq file to a fasta file.
Headers of the final fasta file will contain the count
"""
import argparse
import os

from typing import Optional
from collections import OrderedDict

try:
    # Load Counter's C helper function if it is available. Uses 30% less memory than Counter()
    from _collections import _count_elements
except ImportError:
    # Uses slower mapping[elem] = mapping.get(elem,default_val)+1 in loop
    from collections import _count_elements


def get_args():
    """
    Get input arguments for collapser
    functions.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input-file', metavar='FASTQFILE', required=True,
                        help='input fastq file to collapse')
    parser.add_argument('-o', '--out-file', metavar='OUTPUTFILE', required=True,
                        help='output file name to use')
    parser.add_argument('-t', '--threshold', type=int, default=0,
                        help='number of sequences needed to keep in'
                        'final fasta file')

    args = parser.parse_args()

    return args


def seq_counter(fastq_file: str) -> dict:
    """Counts number of times each sequence appears.
    Skip quality scores since it is quality filtered
    data.

    Inputs:
        fastq_file: A trimmed, quality filtered fastq
                    file

    Outputs:
        seqs: return a dict of {seq: count}
    """

    with open(fastq_file, 'rb') as f:
        def line_generator():    # Generator function for every 4th line (fastq sequence line) of file
            while f.readline():  # Sequence identifier
                # Sequence (Binary -> ASCII extract every 4th from 1st line, newline removed)
                yield f.readline()[:-1].decode("utf-8")
                f.readline()     # "+"
                f.readline()     # Quality Score

        # Count occurrences of unique elements, and populate (merge) with seqs
        seqs = OrderedDict()
        _count_elements(seqs, line_generator())

    seqs.pop("", None)  # Remove blank line counts from the dictionary
    return seqs


def seq2fasta(seqs: dict, out_prefix: str, thresh: int = 0) -> None:
    """
    Turns a sequence count dict into a fasta file.

    Inputs:
        seqs: dictionary containing sequences and counts
        out_file: string to name the output file

    Outputs:
        Writes a fasta file with the format:
            >seq_N_xCOUNTS
            SEQUENCE
    """

    # TODO: move these same checks to main(). Current state is sufficient for script use, but command line use will still spend time before notifying about the problem
    assert out_prefix is not None, "Collapser critical error: an output file must be specified."
    assert thresh >= 0, "An invalid threshold was specified."

    # Check that we'll be able to write results before we spend time on the work
    out_file = out_prefix + "_collapsed.fa"
    low_count_file = out_prefix + "_collapsed_lowcounts.fa"
    for file in [out_file, low_count_file]:
        if os.path.isfile(file):
            print(f"Collapser critical error: {file} already exists.")
            return

    # >seq_INDEX_xCOUNT
    # SEQUENCE
    def to_fasta_record(x):
        # x[0]=index, x[1][1]=sequence count, x[1][0]=sequence
        return '\n'.join([
            f">seq_{x[0]}_x{x[1][1]}",
            f"{x[1][0]}"
        ])

    def get_above_thresh(x):
        return x[1][1] > thresh

    def get_below_thresh(x):
        return x[1][1] <= thresh

    above_thresh = filter(get_above_thresh, enumerate(seqs.items()))
    below_thresh = filter(get_below_thresh, enumerate(seqs.items()))

    with open(out_file, 'w') as fasta:
        if thresh == 0:  # No filtering required
            fasta.write('\n'.join(map(to_fasta_record, enumerate(seqs.items()))))
        else:
            with open(low_count_file, 'w') as lcf:
                fasta.write('\n'.join(map(to_fasta_record, above_thresh)))
                lcf.write('\n'.join(map(to_fasta_record, below_thresh)))


def main():
    """
    main routine
    """
    args = get_args()
    seqs = seq_counter(args.input_file)
    seq2fasta(seqs, args.out_file, args.threshold)


if __name__ == '__main__':
    main()
