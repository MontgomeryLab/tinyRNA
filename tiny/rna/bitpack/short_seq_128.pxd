from .short_seq_util cimport *

# Cython needs help recognizing 128bit integers
cdef extern from *:
    ctypedef unsigned long long uint128_t "__uint128_t"

# Reusable buffer for unmarshalling
cdef char out_ascii_buffer_64[64]

cdef class ShortSeq128:                # 16 bytes (PyObject_HEAD)
    cdef uint128_t _packed             # 16 bytes
    cdef uint8_t _length               # 1 byte
                                       # Remain: 15 bytes pad/free

# Forward declaration to allow for cimport
cdef uint128_t _marshall_bytes_128(uint8_t* seq_bytes, uint8_t length, bint with_length=*) nogil