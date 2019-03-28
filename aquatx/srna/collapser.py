""" Collapse sequences from a fastq file to a fasta file.
Headers of the final fasta file will contain the count
"""

from collections import Counter
import argparse

def get_args():
    """
    Get input arguments for collapser
    functions.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-file', metavar='FASTQFILE', required=True,
                        help='input fastq file to collapse')
    parser.add_argument('-o', '--out-file', metavar='OUTPUTFILE',
                        help='output file name to use')
    parser.add_argument('-t', '--threshold', type=int, default=0,
                        help='number of sequences needed to keep in'
                        'final fasta file')
    parser.add_argument('-k', '--keep-low-counts', metavar='FILENAME', 
                        help='keep sequences not meeting threshold in'
                        'a separate file')

    args = parser.parse_args()

    return args

def seq_counter(fastq_file):
    """
    Counts number of times each sequence appears.
    Skip quality scores since it is quality filtered
    data.

    Inputs:
        fastq_file: A trimmed, quality filtered fastq
                    file

    Outputs:
        seqs: return a dict of {seq: count}
    """
    seqs = Counter()

    with open(fastq_file) as f:
        while True:
            # Read 4 lines per loop for a fastq format
            _ = f.readline() # unused
            seq = f.readline()
            _ = f.readline() # unused
            qual = f.readline() # qual should be last line for a non-corrupt fq

            # Indicates EOF
            if not qual:
                break
            
            seqs[seq] += 1

    return seqs

def seq2fasta(seqs, out_file, thresh=0, low_count_file=None):
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
    # Count keeps track of sequence number for headers
    count = 0
    lowseq = Counter()
    lowseqname = dict()
    # Append to a fasta file - need to add a check to see if file exists
    with open(out_file, "w") as f:
        for key, value in seqs.items():
            if value > thresh:
                f.writelines('>seq_{0}_x{1}\n{2}'.format(count, value, key))
            else:
                lowseq[key] += 1
                lowseqname[key] = count
            count += 1
    if low_count_file is not None:
        with open(low_count_file, "w") as f:
            for key, value in lowseq.items():
                f.writelines('>seq_{0}_x{1}\n{2}'.format(lowseqname[key], value, key))

def main():
    """
    main routine
    """
    args = get_args()
    seqs = seq_counter(args.input_file)

    if args.keep_low_counts:
        seq2fasta(seqs, args.out_file, args.threshold, args.keep_low_counts)
    else:
        seq2fasta(seqs, args.out_file, args.threshold)


if __name__ == '__main__':
    main()
