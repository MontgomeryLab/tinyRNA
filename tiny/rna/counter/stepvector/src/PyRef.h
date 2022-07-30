#include "Python.h"
#include <iostream>

namespace PyPtr {
    class PyRef {
        PyObject* obj;

        public:
            PyObject* get() const {return obj;}

            // Default constructor
            PyRef() {
                obj = NULL;
            }

            // Obj constructor
            PyRef(PyObject* obj) {
                if (PySet_Check(obj)){
                    this->obj = obj;
                } else {
                    // PySet_New treats the argument as an iterable
                    // We don't want that, e.g. if adding a tuple
                    // PySet_Add simply adds without iteration
                    this->obj = PySet_New(NULL);
                    PySet_Add(this->obj, obj);
                }
                Py_XINCREF(this->obj);
            }

            // Destructor
            ~PyRef() {
                Py_XDECREF(obj);
            }

            // Copy constructor
            PyRef(const PyRef& other) {
                obj = PySet_New(other.obj);
                Py_XINCREF(obj);
                Py_XDECREF(other.obj);
            }

            // Move constructor
            PyRef(PyRef&& other) {
                obj = other.obj;
                Py_XDECREF(other.obj);
            }

            // Assignment
            PyRef& operator=(const PyRef& other) {
                Py_XDECREF(obj);
                obj = other.obj;
                Py_XINCREF(obj);
                return *this;
            }

            // Move assignment
            PyRef& operator=(PyRef&& other) {
                Py_XDECREF(obj);
                obj = PySet_New(other.obj);
                Py_XINCREF(obj);
                other.obj = NULL;
                return *this;
            }

            PyRef& operator+= (const PyRef& other) {
                if (PySet_Check(other.obj)){
                    _PySet_Update(this->obj, PyObject_GetIter(other.obj));
                } else {
                    PySet_Add(this->obj, other.obj);
                }
                return *this;
            }

            PyRef operator+ (const PyRef& other) {
                PyRef result(obj);
                if (PySet_Check(other.obj)){
                    _PySet_Update(result.obj, PyObject_GetIter(other.obj));
                } else {
                    PySet_Add(result.obj, other.obj);
                }
                Py_XINCREF(result.obj);
                return result;
            }

            bool operator== (const PyRef &other) const {
                return PyObject_RichCompareBool(obj, other.obj, Py_EQ);
            }
    };
}