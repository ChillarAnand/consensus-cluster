// 
// Cost function for simulated annealing methods
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
sa_cost(PyObject *self, PyObject *args)
{
    
	PyObject *order_list;
	PyArrayObject *sim_matrix;

    if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &order_list, &PyArray_Type, &sim_matrix))
        return NULL;

	if (sim_matrix->nd != 2 || sim_matrix->descr->type_num != PyArray_DOUBLE) {
		PyErr_SetString(PyExc_ValueError, "Matrix input must be two-dimensional and of type float");
		return NULL;
	}

	int order_list_size = PyList_Size(order_list);

	if (order_list_size != sim_matrix->dimensions[0]) {
		PyErr_SetString(PyExc_TypeError, "Indices array and similarity matrix have different lengths!");
		return NULL;
	}

	/* FIXME: I hear boundary checking on the indices is a good idea... */

	int i;
   	long sample1_idx, sample2_idx;
	double result = 0;
	
	/* sum([ matrix[indices[i]][indices[i+1]] for i in range(len(matrix) - 1) ])**2 */

	for( i = 0; i < (order_list_size - 1); i++ ) {
		/* data is a pointer to the start of the array data
		 * strides is the address offset between successive data elements in a contiguous array at dimension [i] */
		
		sample1_idx = PyInt_AsLong(PyList_GetItem(order_list, i));
		sample2_idx = PyInt_AsLong(PyList_GetItem(order_list, i+1));

		result += *(double *)(sim_matrix->data + sample1_idx*sim_matrix->strides[0] + sample2_idx*sim_matrix->strides[1]);
	}

	return PyFloat_FromDouble(pow(result, 2));
}

static PyMethodDef SAMethods[] = {
    {"c_energy",  sa_cost, METH_VARARGS,
     "Implements the matrix reordering cost function.  USAGE: c_energy(indices_order, similarity_matrix)"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

void
initsa(void)
{
    Py_InitModule3("sa", SAMethods, "Methods for SA speed enhancement.  Currently only implements the energy/cost function.");
	import_array();
}


