# cython: language_level = 3, language=c++, profile=False, linetrace=False
import sys
import time

import cython  # For function decorators
from cython.operator cimport dereference as deref

# Avoids initial increment from zero (reuse existing one-object instead)
cdef object one = PyLong_FromSize_t(1)


cpdef ShortSeqCounter fast_read_and_count(str filename):
    cdef ShortSeqCounter counts = ShortSeqCounter()
    cdef vector[PyObject *] seqs
    # cdef vector[char *] seqs

    t1 = time.time()
    read_fastq(filename.encode('utf-8'), seqs)
    # read_fastq_raw(filename.encode('utf-8'), seqs)
    t2 = time.time()
    counts.count_short_seqs(seqs)
    # counts.count_items_cc(seqs)
    t3 = time.time()

    print(f"{t2-t1}s to read {seqs.size()} total seqs, and {t3 - t2}s to count {len(counts)} unique sequences")
    return counts


cdef class ShortSeqCounter(dict):
    def __init__(self, source=None):
        super().__init__()

        if type(source) is list:
            self.count_items_py(source)

    cdef count_short_seqs(self, vector[PyObject *] &short_seqs):
        for seq in short_seqs:
            self._count_sequence(<object>seq)

    @cython.boundscheck(False)
    cdef count_items_cc(self, vector[char *] &raw_lines):
        for seqchars in raw_lines:
            seq = ShortSeqFactory.from_chars(seqchars)
            self._count_sequence(seq)

    @cython.boundscheck(False)
    cpdef count_items_py(self, list it):
        cdef bytes seqbytes

        for i in range(len(it)):
            seqbytes = <bytes> PyList_GET_ITEM(it, i)
            seq = ShortSeqFactory.from_bytes(seqbytes)
            self._count_sequence(seq)

    cdef inline _count_sequence(self, object seq):
        cdef PyObject *oldval

        seqhash = deref(<ShortSeq*>seq)._packed
        oldval = _PyDict_GetItem_KnownHash(self, seq, seqhash)

        if oldval == NULL:
            if PyErr_Occurred():
                raise Exception("Something went wrong while retrieving sequence count.")
            if _PyDict_SetItem_KnownHash(self, seq, one, seqhash) < 0:
                raise Exception("Something went wrong while setting a new sequence count.")
        else:
            if _PyDict_SetItem_KnownHash(self, seq, <object> oldval + 1, seqhash) < 0:
                raise Exception("Something went wrong while setting an incremented sequence count.")


cdef class ShortSeqFactory:
    @staticmethod
    @cython.always_allow_keywords(False)
    cdef inline object from_bytes(bytes seqbytes):
        cdef char* sequence = deref(<PyBytesObject *> seqbytes).ob_sval
        cdef uint8_t length = Py_SIZE(seqbytes)
        return ShortSeqFactory._factory(sequence, length)

    @staticmethod
    @cython.always_allow_keywords(False)
    cdef inline object from_chars(char* sequence):
        cdef uint8_t length = strlen(sequence) - 1  # omit newline char
        return ShortSeqFactory._factory(sequence, length)

    @staticmethod
    cdef inline object _factory(char* sequence, uint8_t length):
        if 32 < length < 64:
            out128 = ShortSeq128.__new__(ShortSeq128)
            (<ShortSeq128>out128)._packed = _marshall_bytes_128(<uint8_t*>sequence, length)
            (<ShortSeq128>out128)._length = length
            return out128
        elif length <= 32:
            out64 = ShortSeq64.__new__(ShortSeq64)
            (<ShortSeq64>out64)._packed = _marshall_bytes_64(<uint8_t*>sequence, length)
            (<ShortSeq64>out64)._length = length
            return out64
        else:
            raise Exception("Sequences longer than 64 bases are not supported.")