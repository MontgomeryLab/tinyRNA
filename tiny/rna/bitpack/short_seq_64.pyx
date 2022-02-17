# cython: language_level = 3, language=c++, profile=False, linetrace=False
import sys

cimport cython
from cython.operator cimport dereference as deref

"""
ShortSeq64: packs sequences up to 32 bases in length
into fixed-size objects using 2-bit encoding.

For a 32 base sequence:
    Memory footprint: 32 bytes
    PyBytes equivalent: 65 bytes (51% reduction)
    PyUnicode equivalent: 81 bytes (60% reduction)
"""

cdef class ShortSeq64:
    @property
    def _length(self) -> uint8_t:
        return self._length

    @_length.setter
    def _length(self, uint8_t val):
        self._length = val

    @property
    def _packed_seq(self) -> uint64_t:
        return self._packed_seq

    @_packed_seq.setter
    def _packed_seq(self, val):
        self._packed_seq = val

    def __hash__(self) -> uint64_t:
        return self._packed_seq

    def __eq__(self, ShortSeq64 other):
        return self._length == other._length and self._packed_seq == other._packed_seq

    def __str__(self):
        return _unmarshall_bytes_64(self._packed_seq, self._length)

    def __len__(self):
        return self._length

# This avoids the (relatively) massive overhead of calling __cinit__()
@cython.always_allow_keywords(False)
cdef inline ShortSeq64 make_ShortSeq64(bytes sequence):
    cdef ShortSeq64 obj = ShortSeq64.__new__(ShortSeq64)
    cdef uint8_t length = Py_SIZE(sequence)
    obj._packed_seq = _marshall_bytes_64(sequence, length)
    obj._length = length
    return obj

@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline uint64_t _marshall_bytes_64(bytes seq_bytes, uint8_t length) nogil:
    cdef uint8_t* sequence = <uint8_t*>(deref(<PyBytesObject*>seq_bytes).ob_sval)
    cdef uint64_t hashed = 0L
    cdef uint8_t i

    for i in reversed(range(length)):
        hashed = (hashed << 2) | table_91[sequence[i]]

    return hashed

@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline unicode _unmarshall_bytes_64(uint64_t enc_seq, size_t length):
    cdef uint8_t i

    for i in reversed(range(length)):
        out_ascii_buffer_32[i] = charmap[enc_seq & mask]
        enc_seq >>= 2

    return PyUnicode_DecodeASCII(out_ascii_buffer_32, length, NULL)