#include "Python.h"
#include <iostream>
#include <cstring>

namespace PyPtr {
    class PyRef {
        PyObject* obj;
        bool debug = false;

        public:
            PyObject* get() const {return obj;}

            // Default constructor
            PyRef() {
                obj = PySet_New(NULL);
            }

            // Obj constructor
            PyRef(PyObject* payload) {
                if (PySet_Check(payload)){
                    obj = payload;
                    Py_XINCREF(obj);
                } else {
                    // PySet_New treats the argument as an iterable
                    // We don't want that, e.g. if adding a tuple
                    // PySet_Add simply adds without iteration
                    obj = PySet_New(NULL);
                    PySet_Add(obj, payload);
                }

                if (debug){
                    std::cout << "Object ctor: ";
                    PyObject_Print(payload, stdout, 0);
                    std::cout << " (obj: " << obj << ")" << std::endl;
                }
            }

            // Destructor
            ~PyRef() {
                Py_XDECREF(obj);
            }

            // Copy constructor
            PyRef(const PyRef& other) {
                obj = PySet_New(other.obj);

                if (debug) {
                    std::cout << "Copy ctor: ";
                    PyObject_Print(other.obj, stdout, 0);
                    std::cout << " (obj: " << other.obj << " -> " << obj << ")";
                    std::cout << std::endl << std::flush;
                }
            }

            // Move constructor
            PyRef(PyRef&& other) {
                if (debug) {
                    std::cout << "Move ctor: ";
                    PyObject_Print(other.obj, stdout, 0);
                    std::cout << " (obj: " << other.obj << ")";
                    std::cout << std::endl << std::flush;
                }

                obj = other.obj;
                Py_XDECREF(other.obj);
            }

            // Assignment
            PyRef& operator=(const PyRef& other) {
                if (debug) {
                    std::cout << "Assignment: ";
                    PyObject_Print(other.obj, stdout, 0);
                    std::cout << " (obj: " << other.obj << ")";
                    std::cout << std::endl << std::flush;
                }

                Py_XDECREF(obj);
                obj = PySet_New(other.obj);
                return *this;
            }

            // Move assignment
            PyRef& operator=(PyRef&& other) {
                if (debug) {
                    std::cout << "Move assignment: ";
                    PyObject_Print(other.obj, stdout, 0);
                    std::cout << " (obj: " << other.obj << ")";
                    std::cout << std::endl << std::flush;
                }

                Py_XDECREF(obj);
                obj = other.obj;
                // Py_XINCREF(obj);
                other.obj = NULL;
                return *this;
            }

            PyRef& operator+= (const PyRef& other) {
                if (debug) {
                    std::cout << "+=" << std::endl;
                }

                if (PySet_Check(other.obj)){
                    _PySet_Update(obj, PyObject_GetIter(other.obj));
                } else {
                    std::cerr << "Error: rhs value holds a non-set type PyObject." << std::endl;
                }
                return *this;
            }

            PyRef operator+ (const PyRef& rhs) {
                if (debug){
                    PyObject_Print(obj, stdout, 0);
                    std::cout << " + ";
                    PyObject_Print(rhs.obj, stdout, 0);
                    std::cout << std::endl;
                }

                PyRef result(obj);
                if (PySet_Check(rhs.obj)){
                    _PySet_Update(result.obj, PyObject_GetIter(rhs.obj));
                } else {
                    std::cerr << "Error: rhs value holds a non-set type PyObject." << std::endl;
                }
                
                return result;
            }

            bool operator== (const PyRef &other) const {
                return PyObject_RichCompareBool(obj, other.obj, Py_EQ);
            }
    };
}