""" 
Pre-process GTF files prior to using the small RNA tool.
Takes a GTF file and collapses to unique features based
on position. 
"""

import argparse
import os.path
import pandas as pd

def get_args():
    """
    Get input arguments

    Requires the user provide a GFF file
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-file', metavar='GTFFILE',
                        help='GTF file to pre-process')

    args = parser.parse_args()

    return args

def check_chr_labels(features, alignments):
    """
    Compare the chromosome labels between the input
    reference annotation and the alignment file
    to make sure the chromosome labeling is consistent
    """

    chr1 = features["chr"].unique()
    chr1.sort()

    chr2 = alignments["chr"].unique()
    chr2.sort()

    is_equal = np.array_equal(chr1, chr2)

    return is_equal

def swap_chroms(features_file):
    """ 
    Make the chromosomes in the feature file
    match the chromosomes in the alignment file
    """
    header = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attr"] 
    features = pd.read_csv(features_file, sep="\t", header=0, names=header)
    features_uniq = features.drop_duplicates()
    chrdict = {'1': 'I', '2': 'II', '3': 'III', '4': 'IV', '5': 'V', '6': 'X', '7': 'MtDNA'}
    features['chr'] = features['chr'].apply(lambda x: "CHROMOSOME_" + chrdict[str(x)])
    features.to_csv('fixed_'+ os.path.basename(features_file), sep='\t', header=header, index=False)

def main():
    """ 
    main routine 
    """
    args = get_args()
    swap_chroms(args.input_file)

if __name__ == '__main__':
    main()


