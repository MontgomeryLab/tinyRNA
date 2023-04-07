# cython: language_level = 3
# cython: profile=False

cimport cython
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from pysam.libcalignedsegment cimport pysam_bam_get_qname
from pysam.libchtslib cimport bam_aux_get, bam_aux2i, BAM_FUNMAP, BAM_FREVERSE

from typing import Callable

cdef tuple cigar_mismatch = (1, 2, 8)

cdef class AlignmentIter:

    cdef AlignmentFile reader
    cdef tuple references
    cdef object dc_callback
    cdef list dc_queue
    cdef bint decollapse
    cdef bint has_nm

    def __init__(self, pysam_reader: AlignmentFile, has_nm: bool, dc_callback: Callable, dc_queue: list):
        self.reader = pysam_reader
        self.references = pysam_reader.header.references
        self.decollapse = dc_callback is not None
        self.dc_callback = dc_callback
        self.dc_queue = dc_queue
        self.has_nm = has_nm

    def __iter__(self):
        return self

    def __next__(self):
        return self._next_alignment()

    cdef dict _next_alignment(self):
        cdef:
            AlignedSegment aln = next(self.reader)                  # Equivalent Python API calls:
            int flag = aln._delegate.core.flag                          # aln.flag
            str seq = aln.query_sequence
            int start = aln._delegate.core.pos                          # aln.reference_start
            int length = aln._delegate.core.l_qseq                      # aln.query_length
            bint strand = (flag & BAM_FREVERSE) == 0                    # aln.is_forward
            str nt5 = AlignmentIter._get_nt5(seq, strand)
            str name = pysam_bam_get_qname(aln._delegate).decode()      # aln.query_name

        if (flag & BAM_FUNMAP) != 0:                                    # aln.is_unmapped
            return self._next_alignment()

        if self.decollapse:
            self._append_to_dc_queue(aln)

        if self.has_nm:
            mismatches = AlignmentIter._get_mismatches_from_nm(aln)     # aln.get_tag("NM")
        else:
            mismatches = AlignmentIter._get_mismatches_from_cigar(aln)

        return {
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

    @staticmethod
    cdef _get_mismatches_from_nm(AlignedSegment aln):
        nm = bam_aux_get(aln._delegate, b"NM")
        if nm == NULL: return 0  # Assume missing tag means NM:i:0
        return bam_aux2i(nm)

    @staticmethod
    cdef _get_mismatches_from_cigar(AlignedSegment aln):
        # Calculate mismatches using the CIGAR string's I, D, and X operations
        return sum(op_len for op, op_len in aln.cigartuples if op in cigar_mismatch)

    @staticmethod
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

    cdef _append_to_dc_queue(self, AlignedSegment aln):
        self.dc_queue.append(aln)
        if len(self.dc_queue) > 100_000:
            self.dc_callback()