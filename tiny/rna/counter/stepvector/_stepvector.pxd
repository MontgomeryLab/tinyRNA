# distutils: language = c++
# cython: language_level = 3

from cpython.ref cimport PyObject
from libcpp.map cimport map as ccmap

ctypedef PyObject* PyPtr

cdef extern from "src/PyRef.h" namespace "PyPtr":
    cdef cppclass PyRef:
        PyRef() except +
        PyRef(PyObject* set_obj) except +
        PyObject* get() const
        PyRef& operator+(const PyRef &po)

cdef extern from "src/StepVector.cpp":
    pass

cdef extern from "src/StepVector.h" namespace "sparse_vectors":
    cdef cppclass _StepVector[T]:
        ctypedef ccmap[long int, T].const_iterator const_iterator

        _StepVector() except +

        long int min_index
        long int max_index

        T operator[](long int i) const
        void set_value(long int start, long int end, T value) except +
        void add_value(long int start, long int end, T value) except +
        void apply_to_values(long int start, long int end, void (*func)(T &val))
        long int num_values() const
        const_iterator get_values(long int start) const
        const_iterator begin() const
        const_iterator end() const

cdef class StepVector:
    cdef _StepVector[PyRef] *c_step
    cdef long int start
    cdef long int stop
    cdef char _typecode
    @staticmethod
    cdef void _construct_wrapper(StepVector obj)
    @staticmethod
    cdef void _set_wrapper(StepVector obj, _StepVector[PyRef] *ref)
