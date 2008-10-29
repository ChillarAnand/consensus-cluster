#include <Python.h>
	
static PyObject *
find_euclid_dist(PyObject *self, PyObject *args)
{
    
	PyObject *python_list1, *python_list2, *listitem1, *listitem2;

    if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &python_list1, &PyList_Type, &python_list2))
        return NULL;

	int list1_size = PyList_Size(python_list1);
	int list2_size = PyList_Size(python_list2);

	if (list1_size != list2_size) {
		PyErr_SetString(PyExc_TypeError, "Lists of unequal length");
		return NULL;
	}

	double result = 0;
	int i = 0;

	for( i=0; i<list1_size; i++ ) {
		listitem1 = PyList_GetItem(python_list1, i);
		listitem2 = PyList_GetItem(python_list2, i);

		result += pow((PyFloat_AsDouble(listitem1) - PyFloat_AsDouble(listitem2)), 2);
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
    Py_InitModule3("euclidean", EuclidMethods, "Implements euclidean distance algorithm for two lists");
}

