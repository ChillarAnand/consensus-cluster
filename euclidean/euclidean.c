// 
// Euclidean distance method for comparison of numpy arrays
// 
// 
// Copyright 2008 Michael Seiler
// Rutgers University
// miseiler@gmail.com
// 
// This file is part of ConsensusCluster.
// 
// ConsensusCluster is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// ConsensusCluster is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with ConsensusCluster.  If not, see <http://www.gnu.org/licenses/>.
// 

#include <Python.h>
#include "arrayobject.h"
	
static PyObject *
find_euclid_dist(PyObject *self, PyObject *args)
{
    
	PyArrayObject *python_array1, *python_array2;

    if (!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &python_array1, &PyArray_Type, &python_array2))
        return NULL;

	if (python_array1->nd != 1 || python_array2->nd != 1 || python_array1->descr->type_num != PyArray_DOUBLE || python_array2->descr->type_num != PyArray_DOUBLE) {
		PyErr_SetString(PyExc_ValueError, "Array inputs must be one-dimensional and of type float");
		return NULL;
	}

	int array1_size = python_array1->dimensions[0];
	int array2_size = python_array2->dimensions[0];

	if (array1_size != array2_size) {
		PyErr_SetString(PyExc_TypeError, "Arrays of unequal length");
		return NULL;
	}

	int i = 0;
	double result = 0;
	
	for( i=0; i<array1_size; i++ ) {
		/* data is a pointer to the start of the array data
		 * strides is the address offset between successive data elements in a contiguous array at dimension [i] */
		result += pow((*(double *)(python_array1->data + i*python_array1->strides[0]) - *(double *)(python_array2->data + i*python_array2->strides[0])), 2);
	}

	return PyFloat_FromDouble(sqrt(result));
}

static PyMethodDef EuclidMethods[] = {
    {"euclidean",  find_euclid_dist, METH_VARARGS,
     "Get the euclidean distance between two list vectors."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

void
initeuclidean(void)
{
    Py_InitModule3("euclidean", EuclidMethods, "Returns the euclidean distance between two 1-dim arrays");
	import_array();
}


