from cython.operator cimport dereference as deref
from libcpp.vector cimport vector
from libc.stdio cimport *
from libc.string cimport strdup, memchr

from cpython.object cimport PyObject
from cpython.ref cimport Py_XINCREF

from .short_seq cimport ShortSeqFactory

cdef extern from "<fcntl.h>" nogil:
    # There is a performance advantage to notifying the kernel of our intent to use
    # sequential access pattern on Linux (posix_fadvise), but doing the same on
    # macOS (fnctl) actually seems to hurt performance...

    int F_RDADVISE
    int F_RDAHEAD

cdef void read_fastq(char* fname, vector[PyObject *] &out)
cdef void read_fastq_raw(char* fname, vector[char *] &out) nogil