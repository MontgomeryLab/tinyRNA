from cpython.long cimport PyLong_FromSize_t
from cpython.list cimport PyList_GET_ITEM
from cpython.exc cimport PyErr_Occurred

from .short_seq_util cimport *
from .short_seq_128 cimport *
from .short_seq_64 cimport *

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