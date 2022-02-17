from libc.stdint cimport uint8_t, uint64_t

from cpython.object cimport Py_SIZE, PyObject
from cpython.ref cimport Py_XDECREF, Py_XINCREF
from cpython.unicode cimport PyUnicode_DecodeASCII

# Constants
cdef uint8_t mask = 0x3
cdef char[4] charmap

"""
This is used for encoding ASCII nucleotide values to 2-bit form.
Indices correspond to the ASCII offset, and table values represent the
encoded value. While table lookup seems like it would be slower than pure 
bitwise ops, this consistently performed better (probably because the useful
portion of the table ends up in a cache line). Supports A,T,U,G,C.
"""
cdef uint8_t[91] table_91

"""
Useful for quick access to ob_sval without interacting with the CPython API
"""
ctypedef struct PyBytesObject:
    PyObject ob_base
    Py_ssize_t ob_size
    Py_hash_t ob_shash
    char ob_sval[1]

"""
For vectorized operations. 
"""
cdef extern from "x86intrin.h" nogil:
    uint64_t _pext_u64(uint64_t __X, uint64_t __Y)