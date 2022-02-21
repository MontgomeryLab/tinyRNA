# cython: language_level=3, language=c++, profile=False, linetrace=False

cimport cython

"""
ShortSeq128: packs sequences up to 64 bases in length
into fixed-size objects using 2-bit encoding.

Memory footprint: 48 bytes
PyBytes equivalent: 97 bytes (51% reduction)
PyUnicode equivalent: 113 bytes (58% reduction)
"""
# Todo
"""
Consider encoding length:
    Max sequence length: 60 bases
    Memory footprint: 32 bytes
    PyBytes equivalent: 93 bytes (66% reduction)
    PyUnicode equivalent: 109 bytes (71% reduction)
"""

cdef class ShortSeq128:
    def __hash__(self):
        return <uint64_t>self._packed

    def __len__(self):
        return self._length

    def __eq__(self, other):
        if type(other) is ShortSeq64:
            return self._length == (<ShortSeq64> other)._length and \
                   self._packed == (<ShortSeq64> other)._packed
        elif type(other) is ShortSeq128:
            return self._length == (<ShortSeq128> other)._length and \
                   self._packed == (<ShortSeq128> other)._packed
        else:
            return False

    def __str__(self):
        return _unmarshall_bytes_128(self._packed, self._length)

@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline uint128_t _marshall_bytes_128(uint8_t* sequence, uint8_t length, bint with_length = False) nogil:
    cdef uint128_t hashed = 0LL
    cdef uint8_t i

    for i in reversed(range(length)):
        hashed = (hashed << 2) | table_91[sequence[i]]

    if with_length:
        hashed = (hashed << 8) | length

    return hashed

@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline unicode _unmarshall_bytes_128(uint128_t enc_seq, uint8_t length = 0, bint with_length = False):
    cdef uint8_t i

    if with_length:
        length = enc_seq & 0xFF
        enc_seq >>= 8

    for i in range(length):
        out_ascii_buffer_64[i] = charmap[enc_seq & mask]
        enc_seq >>= 2

    return PyUnicode_DecodeASCII(out_ascii_buffer_64, length, NULL)