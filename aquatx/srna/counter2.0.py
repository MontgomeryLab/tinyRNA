import argparse
import itertools
import re
import sys

import HTSeq

_re_attr_main = re.compile("\s*([^\s=]+)[\s=]+(.*)")
_re_attr_empty = re.compile("^\s*$")

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
    parser.add_argument('-o', '--out-prefix', metavar='OUTPUTPREFIX',
                        help='output prefix to use for file names')
    parser.add_argument('-t', '--intermed-file', action='store_true',
                        help='Save the intermediate file containing all alignments and'
                             'associated features.')
    parser.add_argument('-s', '--selectors', nargs='+', required=True,
                        help='selection strings, separated by spaces.')

    args = parser.parse_args()

    return args

def parse_selection_string(selector: str) -> list:
    # "CSR;antisense;/path/to/file;1;C,G,U;all"
    return selector.split(';')

def parse_GFF_attribute_string(attrStr, extra_return_first_value = False):
    """Parses a GFF attribute string and returns it as a dictionary.

    This is a slight modification of the same method found in HTSeq.features.
    It has been adapted to allow features to have multiple classes, which are
    stored as a list rather than a comma separated string. This should save
    some CPU time down the road.

    If 'extra_return_first_value' is set, a pair is returned: the dictionary
    and the value of the first attribute. This might be useful if this is the
    ID.
    """
    if attrStr.endswith("\n"):
        attrStr = attrStr[:-1]
    d = {}
    first_val = "_unnamed_"
    for (i, attr) in zip(
            itertools.count(),
            HTSeq._HTSeq.quotesafe_split(attrStr.encode())):
        attr = attr.decode()
        if _re_attr_empty.match(attr):
            continue
        if attr.count('"') not in (0, 2):
            raise ValueError(
                "The attribute string seems to contain mismatched quotes.")
        mo = _re_attr_main.match(attr)
        if not mo:
            raise ValueError("Failure parsing GFF attribute line")
        val = mo.group(2)
        if val.startswith('"') and val.endswith('"'):
            val = val[1:-1]
    # Modification begin
        if mo.group(1) == 'Class':
            val = val.split(',')
            d[sys.intern(mo.group(1))] = val
        else:
            d[sys.intern(mo.group(1))] = sys.intern(val)
    # Modification end
        if extra_return_first_value and i == 0:
            first_val = val
    if extra_return_first_value:
        return (d, first_val)
    else:
        return d

def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2

def main():
    # Step 1: Get command line arguments.
    #args = get_args()

    # Need to patch the GFF attribute parser for our particular GFF structure (list fields)
    setattr(HTSeq.features, 'parse_GFF_attribute_string', parse_GFF_attribute_string)
    gff = HTSeq.GFF_Reader("../../tests/testdata/cel_ws279/c_elegans.PRJNA13758.WS279.chr1.gff3")

    attrs = ['Name', 'interpolated_map_position', 'sequence_name', 'biotype', 'so_term_name', 'curie', 'Alias', 'Class']
    feature_scan = HTSeq.make_feature_genomicarrayofsets(gff,'ID', additional_attributes=attrs)
    features = feature_scan['features']
    attributes = feature_scan['attributes']
    feature_attr = sorted(attributes.keys())

    # CIGAR match characters (including alignment match, sequence match, and
    # sequence mismatch
    com = ('M', '=', 'X')
    counts = {feat: 0 for feat in feature_attr}
    empty = ambiguous = notaligned = nonunique = 0

    read_seq = HTSeq.BAM_Reader('/home/t8/Projects/PycharmProjects/aquatx-srna/tests/run_directory/Lib303_test_aligned_seqs.sam')

    # The previous means of reading in our SAM file, placed here only for testing purposes
    for aln_bundle in HTSeq.bundle_multiple_alignments(read_seq):
        print(len(aln_bundle))

    read: HTSeq.SAM_Alignment
    for read in read_seq:
        if not read.aligned:
            notaligned += 1
            continue
        try:
            if read.optional_field('NH') > 1:
                print(f"Not unique: {read.read}")
                nonunique += 1
        except KeyError:
            pass
        iv_seq = (co.ref_iv for co in read.cigar if co.type in com and co.size > 0)

        # Multi-mapping
        fs = set()
        for iv in iv_seq:
            if iv.chrom not in features.chrom_vectors:
                empty += 1
            for iv2, fs2 in features[iv].steps():
                fs = fs.union(fs2)

        if fs is None or len(fs) == 0:  # No feature
            empty += 1
        elif len(fs) > 1:
            ambiguous += 1

        # Finally, count
        for fsi in list(fs):
            counts[fsi] += 1

    # print(counts)
    print(f"empty {empty} ambiguous {ambiguous} notaligned {notaligned} nonunique {nonunique}")

if __name__ == '__main__':
    main()