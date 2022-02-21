# cython: language_level = 3, language=c++, profile=False, linetrace=False

cimport cython

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
    def _packed(self) -> uint64_t:
        return self._packed

    @_packed.setter
    def _packed(self, val):
        self._packed = val

    def __hash__(self) -> uint64_t:
        return self._packed

    def __eq__(self, other):
        if type(other) is ShortSeq64:
            return self._length == (<ShortSeq64>other)._length and \
                   self._packed == (<ShortSeq64>other)._packed
        elif type(other) is ShortSeq128:
            return self._length == (<ShortSeq128>other)._length and \
                   self._packed == (<ShortSeq128>other)._packed
        else:
            return False

    def __str__(self):
        return _unmarshall_bytes_64(self._packed, self._length)

    def __len__(self):
        return self._length

@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline uint64_t _marshall_bytes_64(uint8_t* sequence, uint8_t length) nogil:
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

    for i in range(length):
        out_ascii_buffer_32[i] = charmap[enc_seq & mask]
        enc_seq >>= 2

    return PyUnicode_DecodeASCII(out_ascii_buffer_32, length, NULL)