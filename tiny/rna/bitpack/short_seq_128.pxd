from .short_seq_util cimport *

# Cython needs help recognizing 128bit integers
cdef extern from *:
    ctypedef unsigned long long uint128_t "__uint128_t"

# Reusable buffer for unmarshalling
cdef char out_ascii_buffer_64[64]

cdef class ShortSeq128:                 # 16 bytes (PyObject_HEAD)
    cdef uint128_t _packed_seq          # 16 bytes
    cdef uint8_t _length                # 1 byte
                                        # Remain: 15 bytes pad/free

cdef ShortSeq128 make_ShortSeq128(bytes sequence)