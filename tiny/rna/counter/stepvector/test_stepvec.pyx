# distutils: language = c++
# cython: language_level = 3

from cpython.ref cimport PyObject
from cython.operator cimport dereference as deref, preincrement as inc, address

from ._stepvector cimport _StepVector, PyRef
ctypedef PyObject* PyPtr

def main():
    cdef PyRef r1, r2, r3, c1, c2

    vctre = new _StepVector[PyRef]()
    f1 = {"featA"}
    f2 = {"featB", "featC"}
    f3 = {"featD"}
    r1 = PyRef(<PyPtr>f1)
    r2 = PyRef(<PyPtr>f2)
    r3 = PyRef(<PyPtr>f3)

    vctre.add_value(11, 99, r1)
    vctre.add_value(22, 44, r2)
    vctre.add_value(33, 98, r3)

    it = vctre.begin()
    while it != vctre.end():
        step = deref(it)
        print(f"{<int> step.first}: {<object> (<PyRef>step.second).get()}")
        inc(it)

if __name__ == '__main__':
    main()