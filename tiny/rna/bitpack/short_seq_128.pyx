# cython: language_level=3, language=c++, profile=False, linetrace=False

cimport cython
from cython.operator cimport dereference as deref

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
        return <uint64_t>self._packed_seq

    def __len__(self):
        return self._length

    def __eq__(self, ShortSeq128 other):
        return self._length == other._length and self._packed_seq == other._packed_seq

    def __str__(self):
        return _unmarshall_128(self._packed_seq, self._length)

# This avoids the (relatively) massive overhead of calling __cinit__()
@cython.always_allow_keywords(False)
cdef inline ShortSeq128 make_ShortSeq128(bytes sequence):
    cdef ShortSeq128 obj = ShortSeq128.__new__(ShortSeq128)
    cdef uint8_t length = Py_SIZE(sequence)
    obj._packed_seq = _marshall_bytes_128(sequence, length, with_length=False)
    obj._length = length
    return obj

@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline uint128_t _marshall_bytes_128(bytes seq_bytes, uint8_t length, bint with_length = False):
    cdef uint8_t* sequence = <uint8_t*>(deref(<PyBytesObject*>seq_bytes).ob_sval)
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
cdef inline unicode _unmarshall_128(uint128_t enc_seq, uint8_t length = 0, bint with_length = False):
    cdef uint8_t i

    if with_length:
        length = enc_seq & 0xFF
        enc_seq >>= 8

    for i in range(length):
        out_ascii_buffer_64[i] = charmap[enc_seq & mask]
        enc_seq >>= 2

    return PyUnicode_DecodeASCII(out_ascii_buffer_64, length, NULL)