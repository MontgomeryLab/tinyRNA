import cython
import sys

from cython.operator cimport dereference as deref, postincrement as inc

cdef class UMI:
    # All UMI subtypes hash by sequence
    # Deduplication then happens naturally via __eq__()
    def __hash__(self):
        return deref(self.seq)

    def __dealloc__(self):
        if self.seq is not NULL:
            PyObject_Free(<void *>self.seq)

cdef class UMI5p(UMI):
    def __eq__(self, UMI5p other):
        return self.length == other.length and \
               self.umi[0] == other.umi[0] and \
               deref(self.seq) == deref(other.seq)

    def __str__(self): return super().__str__()

cdef class UMI3p(UMI):
    def __eq__(self, UMI3p other):
        return self.length == other.length and \
               self.umi[1] == other.umi[1] and \
               deref(self.seq) == deref(other.seq)

cdef class UMIboth(UMI):
    def __eq__(self, UMIboth other):
        return self.length == other.length and \
               self.umi[0] == other.umi[0] and \
               self.umi[1] == other.umi[1] and \
               deref(self.seq) == deref(other.seq)


cdef class UMIFactory:
    def __init__(self, **kwargs):
        self.len_5p = len_5p = kwargs.get('len_5p', 0)  # Default init arg support is currently still an alpha feature
        self.len_3p = len_3p = kwargs.get('len_3p', 0)

        if len_5p and len_3p:
            self._factory = self._factory_both
        elif len_5p:
            self._factory = self._factory_5p
        elif len_3p:
            self._factory = self._factory_3p
        else:
            raise Exception("At least one UMI length is required.")

    cpdef object from_bytes(self, bytes read):
        cdef char* read_chars = deref(<PyBytesObject *> read).ob_sval
        cdef uint8_t length = Py_SIZE(read_chars)
        return self._factory(self, read_chars, length)

    cdef inline object from_chars(self, char* read):
        cdef uint8_t length = strlen(read) - 1  # omit newline char
        return self._factory(self, read, length)

    cdef inline object _factory_5p(self, char* read, uint8_t length):
        cdef UMI5p obj = UMI5p.__new__(UMI5p)
        obj.seq = _marshall_bytes_256(<uint8_t *>read, length)
        obj.umi[0] = 0
        return obj

    cdef inline object _factory_3p(self, char * read, uint8_t length):
        cdef UMI3p obj = UMI3p.__new__(UMI3p)
        obj.seq = _marshall_bytes_256(<uint8_t *>read, length)
        obj.umi[1] = 0
        return obj

    cdef inline object _factory_both(self, char * read, uint8_t length):
        cdef UMIboth obj = UMIboth.__new__(UMIboth)
        obj.seq = _marshall_bytes_256(<uint8_t *>read, length)
        obj.umi[0] = 0
        obj.umi[1] = 0
        return obj

cdef uint64_t pext_mask_64 = 0x0606060606060606
cdef inline uint64_t* _marshall_bytes_256(uint8_t* seq_bytes, uint8_t length):
    cdef:
        uint64_t* sequence = reinterpret_cast[llstr](seq_bytes)
        uint8_t num_blocks = <uint8_t>ceil(length / 32)
        uint8_t block_overhang = <uint8_t>ceil((length % 32) / 8)   # In units of pext (
        uint8_t pext_overhang = length % 8
        uint8_t i, j
        uint64_t* hashed = <llstr>PyObject_Calloc(num_blocks, 8)

    if hashed is NULL:
        PyErr_NoMemory()
        return NULL

    if num_blocks > 1:
        print("bang")

    for i in range(0, length // 8, 4):
        block = 0LL
        for j in reversed(range(i, (i + 4))):
            block = (block << 16) | _pext_u64(sequence[j], pext_mask_64)

        hashed[i // 4] = block

    if block_overhang:
        i += 4
        block_idx = i // 4
        for j in range(i, i + block_overhang):
            hashed[block_idx] = (hashed[block_idx] << 16) | _pext_u64(sequence[j], pext_mask_64)

    if pext_overhang:
        # Mask lower n bits (n = 2 * pext_overhang)
        # This is making me rethink the benefits of using pext for the last lil bit...
        hashed[block_idx] &= (1 << pext_overhang * 2) - 1

    for i in range(num_blocks):
        printbin(f"Block {i}: ", hashed[i], 64, 2)

    return hashed

cdef inline uint32_t _marshall_bytes_32(uint8_t* seq_bytes, uint8_t length):
    cdef uint32_t* sequence = reinterpret_cast[istr](seq_bytes)
    cdef uint32_t hashed = _pext_u32(sequence[0], pext_mask_32)

    # hashed >>= (16 - length) * 2
    hashed &= (1 << ((16 - length) * 2)) - 1
    hashed = (hashed << 4) | length

    return hashed

cdef inline unicode _unmarshall_bytes_256(uint64_t* enc_seq, uint8_t length):
    cdef:
        uint8_t full_blocks = <uint8_t> floor(length / 32)
        uint8_t i, j

    for i in range(full_blocks):
        # for j in range()
        pass