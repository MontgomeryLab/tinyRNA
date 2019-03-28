#! /usr/bin/env python
"""
Small RNA counter

This script counts small RNA sequencing data using the HTSeq API. It assumes the
format of feature files are GFF3 and can use a SAM/BAM alignment file. It allows
for multiple feature file inputs, associated mask files to avoid double counting
certain features (ie miRNA within a coding region), and whether or not to count
both sense and antisense reads. The output is appropriate for use in other DEG
programs such as DESeq2. Summary statistics are also produced.
"""
from collections import Counter
import argparse
import numpy as np
import pandas as pd
import HTSeq

def get_args():
    """
    Get input arguments from the user/command line.

    Requires the user provide a SAM file, a GFF file, and an output
    prefix to save output tables and plots using.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-file', metavar='SAMFILE', required=True,
                        help='input sam file to count features for')
    parser.add_argument('-r', '--ref-annotations', metavar='GTFFILE', nargs='+', required=True,
                        help='reference gff3 files with annotations to count.')
    parser.add_argument('-m', '--mask-file', metavar='MASKFILE', nargs='+', default=None,
                        help='reference gff3 files with annotations to mask from counting.')
    parser.add_argument('-o', '--out-prefix', metavar='OUTPUTPREFIX',
                        help='output prefix to use for file names')
    parser.add_argument('-a', '--antisense', nargs='+', default=None,
                        help='also count reads that align to the antisense'
                             'strand and store in a separate file.')
    parser.add_argument('-t', '--intermed-file', action='store_true',
                        help='Save the intermediate file containing all alignments and'
                             'associated features.')

    args = parser.parse_args()

    return args

def create_ref_array(ref_file, class_counts, feat_counts, mask_file=None, stranded=True):
    """
    Creates the array of features to count from a reference gff3 file. Masks reads from array
    if desired.

    Inputs:
      ref_file: The reference gff3 file with features to counts.
      class_counts: The dictionary for counting classes to assign a value of 0 to
      feat_counts: The dictionary for counting features to assign a value of 0 to
      mask_file: The associated file with features to mask from counting. Default: None
      stranded - Boolean indicating if only sense of a feature is counted. Default: True

    Outputs:
      ref_array - the HTSeq Genomic array of sets containing features and mask features.
    """
    # Read in the gff files
    feat_gff = HTSeq.GFF_Reader(ref_file)
    # Initialize feature array
    feat_array = HTSeq.GenomicArrayOfSets("auto", stranded=stranded)

    # Add all features in the feature file to the array along with class information
    for feat in feat_gff:
        feat_array[feat.iv] += "class_" + feat.type + "_feature_" + feat.attr["ID"]
        # Set value in Counter dicts to 0 so the final output contains all features, even if
        # a library contains no reads for that feature. Required for future normalization.
        class_counts[feat.type] = 0
        feat_counts[feat.attr["ID"]] = 0
    # Add mask features so intervals that overlap have > 1 feature and aren't counted
    if mask_file is not None:
        mask_gff = HTSeq.GFF_Reader(mask_file)
    # Add all masked features that overlap with existing features in array
        for mask in mask_gff:
            # mark features as mask to distinguish them later
            # might make sense to to step through feature array & only add if the mask overlaps
            # with a features
            feat_array[mask.iv] += "class_" + mask.type + "_mask_" + mask.attr["ID"]

    return feat_array, class_counts, feat_counts

def create_ref_dict(ref_files, stranded=None, mask_files=None):
    """
    Creates a dictionary of reference genomic arrays for multiple inputs to later use for
    assigning counts to features.

    Inputs:
        ref_files: List of reference files to count features of.
        mask_file: List of reference files to mask features from associated ref_file. Must be
                  in the same position as ref_file to mask from. Use None for no mask file.
                  Default is None when no mask files are used.
        stranded: List of booleans indicating whether these features should be counted stranded
                  or not. Default is only count sense strands
    Output:
        ref_array_dict: a dictionary containing all feature arrays to be counted.
    """
    ref_array_dict = {}
    class_counts = Counter()
    feat_counts = Counter()

    # Set up mask files list
    if mask_files is not None:
        mask_files = [None if m in 'None' else m for m in mask_files]
    else:
        mask_files = [None for m in ref_files]

    if stranded is None:
        stranded = [True for s in ref_files]
    else:
        stranded = [True if m == 'true' else False for m in stranded]

    try:
        ref_mask_files = zip(ref_files, mask_files, stranded)
    except (ValueError, TypeError) as er:
        print("Length of reference files, mask files, and strand input lists are uneven.")
        raise er

    # populate dict with reference arrays
    for rf, mf, st in ref_mask_files:
        ref_array_dict[rf], class_counts, feat_counts = create_ref_array(rf, class_counts,
                                                                         feat_counts, mf, st)

    return ref_array_dict, class_counts, feat_counts

def assign_features(aln, ref_array_dict):
    """
    Finds a class and a feature that overlaps with the alignment of interest

    Inputs:
        aln: the alignment
        ref_array_dict: the dictionary of feature arrays to check

    Output:
        aln_feats: List of unique features that the alignment corresponds to
        aln_classes: List of unique classes that the alignment corresponds to
    """
    aln_feats = list()
    aln_classes = list()

    # Check all reference arrays for overlapping features
    for ref_file, ref_array in ref_array_dict.items():
        gene_ids = set()
        for iv, val in ref_array[aln.iv].steps():
            gene_ids |= val

            # Assign only if it's one feature per interval
            if len(gene_ids) == 1:
                gene_id = list(gene_ids)[0]
                if gene_id.split('_')[2] == 'mask':
                    continue
                else:
                    aln_feats.append(gene_id.split('_')[3])
                    aln_classes.append(gene_id.split('_')[1])

    # Assign category if no features are found
    if not aln_feats:
        aln_feats.append('_no_feature')
        aln_classes.append('_no_class')

    # Drop duplicate features or classes
    aln_feats = np.unique(np.array(aln_feats))
    aln_classes = np.unique(np.array(aln_classes))

    return aln_feats, aln_classes

def tally_feature_counts(sam_alignment, ref_array_dict, class_counts, feat_counts,
                         stats_out, write=False, outfile=None):
    """
    Tally the counts appropriately for different features and classes of small RNAs.

    Inputs:
        sam_alignment: The sam/bam alignment file
        ref_array_dict: the dictionary containing reference genomic arrays
        stats_out: file to write summary stats to
        write: boolean indicating whether the full feature information should be written
               Default is False.
        outfile: the file handle to write to. Default is none, write must be True to write.

    Outputs:
        class_counts: A dataframe containing counts per possible class
        feature_counts: A dataframe containing counts per feature
        nt_len_mat: A dataframe containing counts per 5' nt x length
    """
    nt_len_mat = {'A': Counter(),
                  'C': Counter(),
                  'T': Counter(),
                  'G': Counter()}
    stats_counts = Counter()

    for aln_bundle in HTSeq.bundle_multiple_alignments(sam_alignment):
        # Calculate counts for multimapping
        dup_counts = int(aln_bundle[0].read.name.split('_x')[1])
        cor_counts = dup_counts/len(aln_bundle)
        stats_counts['_unique_sequences_aligned'] += 1
        stats_counts['_aligned_reads'] += dup_counts
        if len(aln_bundle) > 1:
            stats_counts['_aligned_reads_multi_mapping'] += dup_counts
        else:
            stats_counts['_aligned_reads_unique_mapping'] += dup_counts

        # fill in 5p nt/length matrix
        nt_len_mat[str(aln_bundle[0].read)[0]][len(aln_bundle[0].read)] += dup_counts

        # bundle counts
        bundle_feats = Counter()
        bundle_class = Counter()

        for aln in aln_bundle:
            aln_feats, aln_classes = assign_features(aln, ref_array_dict)
            if write and outfile is not None:
                aln_str = '\t'.join([str(aln.read), str(cor_counts), aln.iv.strand,
                                     str(aln.iv.start), str(aln.iv.end), ';'.join(aln_classes),
                                     ';'.join(aln_feats)])
                outfile.write(aln_str + '\n')

            if len(aln_classes) > 1:
                bundle_class["ambiguous"] += cor_counts
            elif len(aln_feats) > 1:
                for feat in aln_feats:
                    bundle_feats[feat] += cor_counts/len(aln_feats)
            else:
                if aln_classes.item() == '_no_feature':
                    stats_counts['_no_feature'] += cor_counts
                else:
                    bundle_class[aln_classes.item()] += cor_counts
                    bundle_feats[aln_feats.item()] += cor_counts

        if len(bundle_class) > 1:
            class_counts["ambiguous"] += sum(bundle_class.values())
            stats_counts['_ambiguous_alignments_classes'] += 1
            stats_counts['_ambiguous_reads_classes'] += dup_counts
        elif len(bundle_feats) > 1:
            key = next(iter(bundle_class))
            stats_counts['_ambiguous_alignments_features'] += 1
            stats_counts['_ambiguous_reads_features'] += dup_counts
            class_counts[key] += bundle_class[key]
            for key, value in bundle_feats.items():
                feat_counts[key] += value
        else:
            try:
                key = next(iter(bundle_class))
                class_counts[key] += bundle_class[key]
                key = next(iter(bundle_feats))
                feat_counts[key] += bundle_feats[key]
                stats_counts['_alignments_unique_features'] += 1
                stats_counts['_reads_unique_features'] += 1
            except StopIteration:
                pass

    with open(stats_out, 'w') as out:
        out.write('Summary Statistics\n')
        for key, value in stats_counts.items():
            out.write('\t'.join([key, str(value) + '\n']))

    return class_counts, feat_counts, nt_len_mat

def main():
    """
    Main routine for small RNA counter script
    """
    # Step 1: Get command line arguments.
    args = get_args()

    # Step 2: Read in SAM or BAM file
    sam_alignment = HTSeq.SAM_Reader(args.input_file)

    # Step 3: Create feature arrays from GFF files
    ref_array_dict, class_counts, feat_counts = create_ref_dict(args.ref_annotations,
                                                                args.antisense,
                                                                args.mask_file)
    print("Processed feature arrays...")

    # Step 4: Assign alignment counts to features
    stats_out = args.out_prefix + '_stats.txt'

    # Save an intermediate file with all assigned features
    if args.intermed_file:
        aln_int_file = args.out_prefix + '_out_aln_table.txt'
        aln_header = '\t'.join(['seq', 'counts', 'strand', 'start', 'end', 'classes', 'features'])
        with open(aln_int_file, 'w') as outfile:
            outfile.write(aln_header + '\n')
            class_counts, feat_counts, nt_len_mat = tally_feature_counts(sam_alignment,
                                                                         ref_array_dict,
                                                                         class_counts,
                                                                         feat_counts,
                                                                         stats_out,
                                                                         write=True,
                                                                         outfile=outfile)
    else:
        # assign features
        class_counts, feat_counts, nt_len_mat = tally_feature_counts(sam_alignment,
                                                                     ref_array_dict,
                                                                     class_counts,
                                                                     feat_counts,
                                                                     stats_out)

    print("Completed feature assignment...")
    class_counts_df = pd.DataFrame.from_dict(class_counts, orient='index').reset_index()
    feat_counts_df = pd.DataFrame.from_dict(feat_counts, orient='index').reset_index()

    print("Writing final count files...")
    class_counts_df.to_csv(args.out_prefix + '_out_class_counts.csv', index=False, header=False)
    feat_counts_df.to_csv(args.out_prefix + '_out_feature_counts.txt', sep='\t', index=False, header=False)
    pd.DataFrame(nt_len_mat).to_csv(args.out_prefix + '_out_nt_len_dist.csv')

if __name__ == '__main__':
    main()
