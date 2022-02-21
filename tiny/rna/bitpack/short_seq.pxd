import cython

from cpython.long cimport PyLong_FromSize_t
from cpython.list cimport PyList_GET_ITEM
from cpython.exc cimport PyErr_Occurred

from libcpp.typeinfo cimport type_info
from libcpp.vector cimport vector

from .short_seq_util cimport *
from .short_seq_128 cimport *
from .short_seq_64 cimport *
from .fast_read cimport *

"""
Private dictionary fast-path methods not currently offered by the Cython wrapper
"""

cdef extern from "Python.h":
    dict _PyDict_NewPresized(int minused)
    PyObject* _PyDict_GetItem_KnownHash(object mp, object key, Py_hash_t hash)
    PyObject* _PyDict_Pop_KnownHash(object mp, object key, Py_hash_t hash, object deflt)
    int _PyDict_SetItem_KnownHash(object mp, object key, object item, Py_hash_t hash)
    int _PyDict_DelItem_KnownHash(object mp, object key, Py_hash_t hash)
    bint _PyDict_Contains_KnownHash(object mp, object key, Py_hash_t hash)


"""
Forward declarations for importing
"""

cpdef ShortSeqCounter fast_read_and_count(str filename)


cdef class ShortSeqCounter(dict):
    cdef count_short_seqs(self, vector[PyObject *])
    cdef count_items_cc(self, vector[char*] it)
    cpdef count_items_py(self, list it)
    cdef _count_sequence(self, object seq)


cdef class ShortSeqFactory:
    @staticmethod
    cdef inline object from_bytes(bytes seqbytes)

    @staticmethod
    cdef inline object from_chars(char* sequence)

    @staticmethod
    cdef inline object _factory(char* sequence, uint8_t length)
