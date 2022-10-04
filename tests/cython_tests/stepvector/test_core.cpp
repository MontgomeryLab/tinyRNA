#include "StepVector.h"
#include "PyRef.h"

using sparse_vectors::_StepVector;
using PyPtr::PyRef;

int main(){
    Py_Initialize();

    _StepVector<PyRef> sv;

    PyObject *c1 = PySet_New(NULL); // {featA, featB}
    PyObject *c2 = PySet_New(NULL); // {featC}
    PyObject *c3 = PySet_New(NULL); // {featD}

    PyObject *sA = PyUnicode_FromString("featA");
    PyObject *sB = PyUnicode_FromString("featB");
    PyObject *sC = PyUnicode_FromString("featC");
    PyObject *sD = PyUnicode_FromString("featD");

    PySet_Add(c1, sA);
    PySet_Add(c1, sB);
    PySet_Add(c2, sC);
    PySet_Add(c3, sD);

    PyRef p1 = PyRef(c1);
    PyRef p2 = PyRef(c2);
    PyRef p3 = PyRef(c3);

    sv.add_value(0, 10, p1);
    sv.add_value(8, 15, p2);
    sv.add_value(9, 20, p3);

    std::cout << std::endl;

    for (auto it = sv.get_values(-1); it != sv.end(); ++it){
        std::cout << it->first << ": ";
        PyObject_Print(it->second.get(), stdout, 0);
        std::cout << std::endl;
    }

    /* Expect:
    -9223372036854775808: set()
    0: {'featA', 'featB'}
    8: {'featA', 'featC', 'featB'}
    9: {'featD', 'featA', 'featC', 'featB'}
    11: {'featD', 'featC'}
    16: {'featD'}
    21: set()
    */
}