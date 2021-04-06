import HTSeq
import sys
import re

# For parse_GFF_attribute_string()
_re_attr_main = re.compile(r"\s*([^\s=]+)[\s=]+(.*)")
_re_attr_empty = re.compile(r"^\s*$")


class Alignment:
    class Sequence:
        def __init__(self, name, seq):
            self.name = name
            self.seq = seq
            self.len = len(seq)

        def __repr__(self):
            return f"{self.name}, {self.seq}, {self.len}"

        def __len__(self):
            return self.len

    def __init__(self, iv, name, seq):
        self.iv = iv
        self.read = self.Sequence(name, seq)

    def __repr__(self):
        return f"{self.read}, {self.iv}"


def read_SAM(file):
    with open(file, 'rb') as f:
        line = f.readline()

        # Skip headers
        while line[0] == ord('@'):
            line = f.readline()

        while line:
            cols = line.split(b'\t')

            strand = "+" if (int(cols[1]) & 16) >> 4 == 0 else "-"
            chrom = cols[2].decode('utf-8')
            name = cols[0].decode('utf-8')
            start = int(cols[3]) - 1
            seq = cols[9]

            iv = HTSeq.GenomicInterval(chrom, start, start + len(seq), strand)
            yield Alignment(iv, name, seq)

            line = f.readline()


def parse_GFF_attribute_string(attrStr, extra_return_first_value=False):
    """Parses a GFF attribute string and returns it as a dictionary.

    This is a slight modification of the same method found in HTSeq.features.
    It has been adapted to allow features to have multiple classes, which are
    split and stored as a tuple rather than a comma separated string. This
    should save some CPU time down the road.

    If 'extra_return_first_value' is set, a pair is returned: the dictionary
    and the value of the first attribute. This might be useful if this is the
    ID.
    """
    attribute_dict = {}
    first_val = "_unnamed_"
    for i, attr in enumerate(HTSeq._HTSeq.quotesafe_split(attrStr.rstrip().encode())):
        attr = attr.decode()
        if _re_attr_empty.match(attr):
            continue
        if attr.count('"') not in (0, 2):
            raise ValueError(
                "The attribute string seems to contain mismatched  quotes.")
        mo = _re_attr_main.match(attr)
        if not mo:
            raise ValueError("Failure parsing GFF attribute line")
        idt = mo.group(1)
        val = mo.group(2)
        if val.startswith('"') and val.endswith('"'):
            val = val[1:-1]
        # Modification: allow for comma separated attribute values
        attribute_dict[sys.intern(idt)] = (sys.intern(val),) \
            if ',' not in val \
            else tuple(c.strip() for c in val.split(','))
        if extra_return_first_value and i == 0:
            first_val = val
    if extra_return_first_value:
        return attribute_dict, first_val
    else:
        return attribute_dict