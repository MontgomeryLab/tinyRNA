from .short_seq_util cimport *

# Reusable buffer for unmarshalling
cdef char out_ascii_buffer_32[32]

cdef class ShortSeq64:               # 16 bytes (PyObject_HEAD)
    cdef uint64_t _packed_seq        # 8 bytes
    cdef uint8_t _length             # 1 byte
                                     # Total: 32 bytes (7 bytes pad/free)

cdef ShortSeq64 make_ShortSeq64(bytes sequence)