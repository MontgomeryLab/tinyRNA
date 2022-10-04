#include "Python.h"
#include <iostream>
#include <cstring>

//#define DEBUG

namespace PyPtr {
    class PyRef {
        PyObject* obj;

        public:
            PyObject* get() const {
                #ifdef DEBUG
                    std::cout << "Retrieve: ";
                    PyObject_Print(obj, stdout, 0);
                    std::cout << " (obj: " << obj << ")" << std::endl;
                #endif
                return obj;
            }

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

                #ifdef DEBUG
                    std::cout << "Object ctor: ";
                    PyObject_Print(payload, stdout, 0);
                    std::cout << " (obj: " << obj << ")" << std::endl;
                #endif
            }

            // Destructor
            ~PyRef() {
                #ifdef DEBUG
                    std::cout << "Object destructor: ";
                    PyObject_Print(obj, stdout, 0);
                    std::cout << " (obj: " << obj << ")" << std::endl;
                #endif
                Py_XDECREF(obj);
            }

            // Copy constructor
            PyRef(const PyRef& other) {
                obj = PySet_New(other.obj);

                #ifdef DEBUG
                    std::cout << "Copy ctor: ";
                    PyObject_Print(other.obj, stdout, 0);
                    std::cout << " (obj: " << other.obj << " -> " << obj << ")";
                    std::cout << std::endl << std::flush;
                #endif
            }

            // Move constructor
            PyRef(PyRef&& other) {
                #ifdef DEBUG
                    std::cout << "Move ctor: ";
                    PyObject_Print(other.obj, stdout, 0);
                    std::cout << " (obj: " << other.obj << ")";
                    std::cout << std::endl << std::flush;
                #endif

                obj = other.obj;
                Py_XDECREF(other.obj);
            }

            // Assignment
            PyRef& operator=(const PyRef& other) {
                #ifdef DEBUG
                    std::cout << "Assignment: ";
                    PyObject_Print(other.obj, stdout, 0);
                    std::cout << " (obj: " << other.obj << ")";
                    std::cout << std::endl << std::flush;
                #endif

                Py_XDECREF(obj);
                obj = PySet_New(other.obj);
                return *this;
            }

            // Move assignment
            PyRef& operator=(PyRef&& other) {
                #ifdef DEBUG
                    std::cout << "Move assignment: ";
                    PyObject_Print(other.obj, stdout, 0);
                    std::cout << " (obj: " << other.obj << ")";
                    std::cout << std::endl << std::flush;
                #endif

                Py_XDECREF(obj);
                obj = other.obj;
                // Py_XINCREF(obj);
                other.obj = NULL;
                return *this;
            }

            PyRef& operator+= (const PyRef& other) {
                #ifdef DEBUG
                    std::cout << "+=" << std::endl;
                #endif

                if (PySet_Check(other.obj)){
                    _PySet_Update(obj, PyObject_GetIter(other.obj));
                } else {
                    std::cerr << "Error: rhs value holds a non-set type PyObject." << std::endl;
                }
                return *this;
            }

            PyRef operator+ (const PyRef& rhs) {
                #ifdef DEBUG
                    PyObject_Print(obj, stdout, 0);
                    std::cout << " + ";
                    PyObject_Print(rhs.obj, stdout, 0);
                    std::cout << std::endl;
                #endif

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