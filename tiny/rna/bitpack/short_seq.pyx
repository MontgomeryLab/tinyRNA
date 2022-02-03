# cython: language_level = 3
# cython: profile=False
# language: c++

import cython
from cython.operator cimport dereference as deref, preincrement as inc

cdef class ShortSeq:                 # 16 bytes (PyObject_HEAD)
    cdef uPY_LONG_LONG _PyLong_hash  # 8 bytes
    cdef uint8_t _length             # 1 byte
    # cdef uint8_t pad[7]            # Remain: 7 bytes pad/free

    def __cinit__(self, bytes sequence):
        cdef Py_ssize_t length
        cdef uPY_LONG_LONG seqhash

        length = Py_SIZE(sequence)
        seqhash = _marshall_bytes(sequence, length)

        self._PyLong_hash = seqhash
        self._length = length

    @property
    def _length(self) -> uint8_t:
        return self._length

    @_length.setter
    def _length(self, uint8_t val):
        self._length = val

    @property
    def _PyLong_hash(self) -> uPY_LONG_LONG:
        return self._PyLong_hash

    @_PyLong_hash.setter
    def _PyLong_hash(self, uPY_LONG_LONG val):
        self._PyLong_hash = val

    def __hash__(self) -> uPY_LONG_LONG:
        return self._PyLong_hash

    def __eq__(self, ShortSeq other):
        return self._length == other._length and self._PyLong_hash == other._PyLong_hash

    def __str__(self):
        return _unmarshall_bytes(self._PyLong_hash, self._length)

    def __len__(self):
        return self._length


cdef uint8_t table_91[91]
table_91[:] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0]


@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline uPY_LONG_LONG _marshall_bytes(bytes seq_bytes, uint8_t length):
    cdef uint8_t* sequence = <uint8_t*>(deref(<PyBytesObject*>seq_bytes).ob_sval)
    cdef uint64_t hashed = 0L
    cdef uint8_t i

    for i in reversed(range(length)):
        hashed = (hashed << 2) | table_91[sequence[i]]

    return hashed


cdef uint8_t mask = 0x3
cdef uint8_t encoded = 0
cdef char out_buff[32]
cdef char charmap[4]
charmap = [b'A', b'C', b'T', b'G']

@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef unicode _unmarshall_bytes(uPY_LONG_LONG enc_seq, size_t length):
    cdef uint8_t i

    for i in reversed(range(length)):
        out_buff[i] = charmap[enc_seq & mask]
        enc_seq >>= 2

    return PyUnicode_DecodeASCII(out_buff, length, NULL)