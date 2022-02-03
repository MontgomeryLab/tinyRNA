from libc.stdint cimport uint8_t, uint64_t

from cpython.long cimport uPY_LONG_LONG, PyLong_FromSize_t
from cpython.dict cimport PyDict_New
from cpython.exc cimport PyErr_Occurred
from cpython.unicode cimport PyUnicode_DecodeASCII
from cpython.object cimport Py_SIZE, PyObject
from cpython.number cimport PyNumber_Add
from cpython.ref cimport Py_XDECREF, Py_XINCREF

cdef size_t MAX_CHARS = sizeof(unsigned long long) * 4

ctypedef struct PyBytesObject:
    PyObject ob_base
    Py_ssize_t ob_size
    Py_hash_t ob_shash
    char ob_sval[1]