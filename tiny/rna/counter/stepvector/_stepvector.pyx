# distutils: language = c++
# cython: language_level = 3

import sys
import HTSeq

from cython.operator cimport dereference as deref, postincrement as cni

cdef class StepVector:

    @classmethod
    def create(cls, *args, **kwargs):
        cdef StepVector self = StepVector(*args, **kwargs)
        StepVector._construct_wrapper(self)
        return self

    def __cinit__(self, length = sys.maxsize, typecode = 'O', start_index = 0):
        self.c_step = NULL
        self._typecode = ord(typecode)
        self.start = start_index
        self.stop = start_index + length

    @staticmethod
    cdef void _construct_wrapper(StepVector obj):
        obj.c_step = new _StepVector[PyRef]()

    @staticmethod
    cdef void _set_wrapper(StepVector obj, _StepVector[PyRef] *ref):
        obj.c_step = ref

    def __setitem__(self, index, value):
        if isinstance(value, StepVector):
            if value.start == index.start and value.stop == index.stop:
                return
            else:
                raise NotImplemented("Stepvector-to-Stepvector assignment still missing")
        if isinstance(index, slice):
            if index.step is not None and index.step != 1:
                raise ValueError("Striding slices (i.e., step != 1) are not supported")
            if index.start is None:
                start = self.start
            else:
                if index.start < self.start:
                    raise IndexError("start too small")
                start = index.start
            if index.stop is None:
                stop = self.stop
            else:
                if index.stop > self.stop:
                    raise IndexError("stop too large")
                stop = index.stop
            # Note the "-1": The C++ object uses closed intervals, but we follow
            # Python convention here and use half-open ones.
            self.c_step.set_value(start, stop-1, PyRef(<PyPtr>value))
        else:
            self.c_step.set_value(index, index, PyRef(<PyPtr>value))

    def get_steps(self, values_only = False, merge_steps = True):
        cdef _StepVector[PyRef].const_iterator startvals
        cdef set value, prevval
        startvals = self.c_step.get_values(self.start)
        prevstart = self.start

        # Iterator bounds check until a better solution is written
        if startvals == self.c_step.end():
            yield (self.start, self.stop, set())
            return
        else:
            prevval = <set>deref(cni(startvals)).second.get()

        while startvals != self.c_step.end():
            pair = deref(cni(startvals))
            stepstart, value = <long int>pair.first, <set>pair.second.get()
            if merge_steps and value == prevval:
                continue
            if self.stop is not None and stepstart >= self.stop:
                if not values_only:
                    yield prevstart, self.stop, prevval
                else:
                    yield prevval
                return
            if not values_only:
                yield prevstart, stepstart, prevval
            else:
                yield prevval
            prevstart, prevval = stepstart, value
        else:
            if not values_only:
                yield prevstart, min(self.stop, self.c_step.max_index), prevval
            else:
                yield prevval

    def __getitem__(self, index):
        cdef StepVector res
        cdef int start
        cdef int stop

        if isinstance(index, slice):
            if index.step is not None and index.step != 1:
                raise ValueError("Striding slices (i.e., step != 1) are not supported")
            if index.start is None:
                start = self.start
            else:
                if index.start < self.start:
                    raise IndexError("start too small")
                start = index.start
            if index.stop is None:
                stop = self.stop
            else:
                if index.stop > self.stop:
                    raise IndexError("stop too large")
                stop = index.stop
            res = StepVector(stop - start, 'O', start)
            StepVector._set_wrapper(res, self.c_step)
            res.start = start
            res.stop = stop
            return res
        else:
            return <set>deref(self.c_step.get_values(index)).second.get()

    def __iter__(self):
        for start, stop, value in self.get_steps():
            for i in range(start, stop):
                yield value

    def __repr__(self):
        if self.start == -sys.maxsize - 1:
            start_s = "-inf"
        else:
            start_s = str(self.start)
        if self.stop == sys.maxsize:
            stop_s = "inf"
        else:
            stop_s = str(self.stop)
        return "<%s object, type '%s', index range %s:%s, %d step(s)>" % (
            self.__class__.__name__, self.typecode(), start_s,
            stop_s, self.num_steps())

    def typecode(self):
        return chr(self._typecode)

    @property
    def start(self):
        return self.start

    @property
    def stop(self):
        return self.stop

    @start.setter
    def start(self, value):
        self.start = value

    @stop.setter
    def stop(self, value):
        self.stop = value

    def __len__(self):
        return self.stop - self.start

    def num_steps(self):
        return self.c_step.num_values()

    def __iadd__(self, value):
        cdef PyRef ref = PyRef(<PyPtr>value)
        self.c_step.add_value(self.start, self.stop-1, ref)
        return self

    # Todo
    def __eq__(self, other): pass
    def __reduce__(self): pass
    def apply(self, func, start=None, stop=None): pass

    @property
    def __class__(self):
        """~~master_of_disguise.exe~~"""
        return HTSeq.StepVector.StepVector