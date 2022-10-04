# distutils: language = c++
# cython: language_level = 3
import sys

from cpython.ref cimport PyObject
from libcpp.pair cimport pair
from cython.operator cimport dereference as deref, preincrement as inc, postincrement as pinc

from tiny.rna.counter.stepvector._stepvector cimport _StepVector, PyRef
ctypedef PyObject* PyPtr

def main():
    cdef PyRef p1, p2, p3

    sv = new _StepVector[PyRef]()

    c1 = {"featA", "featB"}
    c2 = {"featC"}
    c3 = {"featD"}
    p1 = PyRef(<PyPtr>c1)
    p2 = PyRef(<PyPtr>c2)
    p3 = PyRef(<PyPtr>c3)

    sv.add_value(0, 10, p1)
    sv.add_value(8, 15, p2)
    sv.add_value(9, 20, p3)

    it = sv.get_values(-1)
    expected = iter([
        ((-1 * sys.maxsize - 1), set()),
        (0, {'featA', 'featB'}),
        (8, {'featA', 'featB', 'featC'}),
        (9, {'featA', 'featB', 'featC', 'featD'}),
        (11, {'featC', 'featD'}),
        (16, {'featD'}),
        (21, set())
    ])

    while it != sv.end():
        step = deref(pinc(it))
        exp_pos, exp_set = next(expected)

        try:
            assert step.first == exp_pos
            assert <object>step.second.get() == exp_set
        except AssertionError as e:
            print(f"Expected {exp_pos}: {exp_set}")
            print(f"Received {step.first}: {<object>step.second.get()}")
            raise

    print("StepVector basic test passed.")

if __name__ == '__main__':
    main()