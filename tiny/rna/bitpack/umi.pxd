from .short_seq_util cimport *

from cpython.mem cimport PyObject_Calloc
from libc.stdlib cimport free

ctypedef object (*factory_method)(UMIFactory self, char* read, uint8_t length)
ctypedef uint64_t* llstr

cdef class UMIFactory:
    cdef factory_method _factory
    cdef uint8_t len_5p
    cdef uint8_t len_3p

    cpdef object from_bytes(self, bytes read)
    cdef inline object from_chars(self, char* read)
    cdef inline object _factory_5p(self, char * read, uint8_t length)
    cdef inline object _factory_3p(self, char * read, uint8_t length)
    cdef inline object _factory_both(self, char * read, uint8_t length)

cdef uint64_t pext_mask = 0x0606060606060606
cdef uint64_t* _marshall_bytes_256(uint8_t* seq_bytes, uint8_t length)

"""
============================================================================
!!! TODO (maybe). After writing the implementation... might be a better way.
    For now, the length is indicated by the attribute "length".
Seq:
    Encoded length: Lower 8 bits of the first uint64_t indicate length
    Max length: 252 bases (28 bases at seq[0]; 32 bases at seq[1:6])

                uint64_t                    uint64_t       uint64_t
    |--------------------------------| ==> |--------| ==> |--------| ==> ...
    |-- NT Bases -----------|-length-|     |NT Bases|     |NT Bases|

============================================================================

UMI (each):
    Encoded length: Lower 4 bits of each UMI indicate the UMI's length
    Max length: 14 bases

         uint32_t
    |----------------|
    |NT Bases |length|
    
============================================================================
"""

cdef class UMI:               # 16 bytes (PyObject head)
    cdef uint32_t umi[2]      # 8 bytes arr (fixed)
    cdef uint64_t *seq        # 8 bytes ptr (~stretchy~) + heap allocation
    cdef uint8_t seq_len      # 1 byte (for now, may eventually encode)
                              # Total: 32 bytes PLUS heap

    cpdef printit(self)

cdef class UMI5p(UMI):
    pass

cdef class UMI3p(UMI):
    pass

cdef class UMIboth(UMI):
    pass