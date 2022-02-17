# cython: language_level = 3, language=c++, profile=False, linetrace=False

import cython  # For function decorators

# Avoids initial increment from zero (reuse existing one-object instead)
cdef object one = PyLong_FromSize_t(1)

cdef class ShortSeqCounter(dict):
    def __init__(self, list iterable):
        super().__init__()
        self.count_items(iterable)

    @cython.boundscheck(False)
    cdef count_items(self, list iterable):
        cdef PyObject* oldval
        cdef bytes seqbytes
        cdef Py_ssize_t i
        cdef uint8_t length

        for i in range(len(iterable)):
            seqbytes = <bytes>PyList_GET_ITEM(iterable, i)
            length = Py_SIZE(seqbytes)

            if 32 < length < 64:
                seq = make_ShortSeq128(seqbytes)
            elif length <= 32:
                seq = make_ShortSeq64(seqbytes)
            else:
                raise Exception()

            seqhash = hash(seq)
            oldval = _PyDict_GetItem_KnownHash(self, seq, seqhash)

            if oldval == NULL:
                if PyErr_Occurred():
                    raise Exception("Something went wrong while retrieving sequence count.")
                if _PyDict_SetItem_KnownHash(self, seq, one, seqhash) < 0:
                    raise Exception("Something went wrong while setting a new sequence count.")
            else:
                if _PyDict_SetItem_KnownHash(self, seq, <object>oldval + 1, seqhash) < 0:
                    raise Exception("Something went wrong while setting an incremented sequence count.")
