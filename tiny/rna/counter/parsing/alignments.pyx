# cython: language_level = 3
# cython: profile=False

cimport cython
from libc.stdint cimport uint32_t
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from pysam.libcalignedsegment cimport pysam_bam_get_qname, pysam_bam_get_cigar, pysam_get_n_cigar
from pysam.libchtslib cimport (
    BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CREF_SKIP, BAM_CPAD,  # CIGAR operator constants
    bam_aux_get, bam_aux2i, bam_aux2Z,  # tags
    BAM_FUNMAP, BAM_FREVERSE,  # flags
    bam_cigar_op,  # Decode CIGAR operator
)

from typing import Callable

cdef class AlignmentIter:

    cdef AlignmentFile reader
    cdef tuple references
    cdef object dc_callback
    cdef list dc_queue
    cdef bint detailed_mismatch
    cdef bint decollapse
    cdef int record_num


    def __init__(
            self,
            pysam_reader: AlignmentFile,
            expected_tags: tuple,
            dc_callback: Callable,
            dc_queue: list,
            detailed_mismatch: bint = False
    ):
        self.reader = pysam_reader
        self.references = pysam_reader.header.references
        self.decollapse = dc_callback is not None
        self.dc_callback = dc_callback
        self.dc_queue = dc_queue
        self.detailed_mismatch = detailed_mismatch
        self.record_num = 0

        if "NM" not in expected_tags:
            raise NotImplementedError("Alignments must include the NM tag.")
        if "MD" not in expected_tags and detailed_mismatch is True:
            raise ValueError(f"Alignments must include the MD tag for detailed mismatches.")

    def __iter__(self):
        return self

    def __next__(self):
        return self._next_alignment()

    cdef dict _next_alignment(self):
        cdef:
            AlignedSegment aln = next(self.reader)                  # Equivalent Python API calls:
            int flag = aln._delegate.core.flag                          # aln.flag
            int start, length
            str seq, nt5, name
            bint strand

        self.record_num += 1

        # Skip unmapped alignments
        while flag & BAM_FUNMAP:                                        # aln.is_unmapped
            aln = next(self.reader)
            flag = aln._delegate.core.flag
            self.record_num += 1

        # Check for incompatible CIGAR ops
        _validate_alignment(aln)

        # Extract relevant info from Pysam's AlignedSegment
        seq = aln.query_alignment_sequence
        start = aln._delegate.core.pos + aln.query_alignment_start      # aln.reference_start
        strand = (flag & BAM_FREVERSE) == 0                             # aln.is_forward
        name = pysam_bam_get_qname(aln._delegate).decode()              # aln.query_name
        mismatches = _get_mismatches_from_nm(aln)
        nt5 = _get_nt5(seq, strand)
        length = len(seq)

        alignment = {
            "Name": name,
            "Length": length,
            "Seq": seq,
            "nt5end": nt5,
            "Chrom": self.references[aln._delegate.core.tid],           # aln.reference_name
            "Start": start,
            "End": start + length,
            "Strand": strand,
            "NM": mismatches
        }

        if self.detailed_mismatch:
            alignment["MD"] = _get_md_string(aln)

        if self.decollapse:
            self._append_to_dc_queue(aln)

        return alignment

    cdef _append_to_dc_queue(self, AlignedSegment aln):
        self.dc_queue.append(aln)
        if len(self.dc_queue) > 100_000:
            self.dc_callback()


cdef int _get_mismatches_from_nm(AlignedSegment aln):                   # aln.get_tag("NM")
    nm = bam_aux_get(aln._delegate, b"NM")
    if nm == NULL:
        raise ValueError("Alignment is missing a required MD tag.")

    return bam_aux2i(nm)

cdef object _get_mismatches_from_cigar(AlignedSegment aln):
    # This method doesn't disambiguate the M operator. Revision required.
    raise NotImplementedError("Mismatches can't be calculated from CIGAR in this version.")

cdef str _get_md_string(AlignedSegment aln):                            # aln.get_tag("MD")
    md = bam_aux_get(aln._delegate, b"MD")
    if md == NULL:
        raise ValueError("Alignment is missing a required MD tag.")

    return bam_aux2Z(md).decode()

cdef str _get_nt5(str seq, bint strand):
    cdef str nt

    if strand:
        return seq[0]
    else:
        nt = seq[-1]
        if nt == "A": return "T"
        elif nt == "T": return "A"
        elif nt == "G": return "C"
        elif nt == "C": return "G"
        else: return nt
        
cpdef _validate_alignment(AlignedSegment aln):
        cdef uint32_t * cigar_p = pysam_bam_get_cigar(aln._delegate)
        cdef uint32_t n_cigar = pysam_get_n_cigar(aln._delegate)
        cdef uint32_t i, op

        for i in range(n_cigar):
            op = bam_cigar_op(cigar_p[i])

            if op == BAM_CSOFT_CLIP:  # S
                raise ValueError(f"Soft clipped bases are not supported at this time")
            elif op == BAM_CHARD_CLIP:  # H
                raise ValueError(f"Hard clipped bases are not supported at this time")
            elif op == BAM_CREF_SKIP:  # N
                raise ValueError(f"Skipped bases are not supported at this time")
            elif op == BAM_CPAD:  # P
                raise ValueError(f"Padded bases are not supported at this time")
