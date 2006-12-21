#include "Python.h"
#include "hfk_class.cpp"

extern "C" void init_hfk(void);

static PyObject *ErrorObject;
static PyObject *callback;

int Progress(const char *message, int percent){
  PyObject *args, *py_result;
  int result;
  if (callback){
    args = Py_BuildValue("(si)", message, percent);
    py_result = PyEval_CallObject(callback, args);
    if (py_result == NULL)
      return -1;
    result = PyInt_AS_LONG(py_result);
    Py_DECREF(py_result);
    return result;
  }
  return 0;
}

static PyObject *HFKhat(PyObject *self, PyObject *args, PyObject *keywds){
  PyObject *py_Xlist;
  PyObject *py_Olist;
  PyObject *py_callback = NULL;
  PyObject *result;
  int quiet = 1;
  int gridsize;
  static char *kwlist[] = {"Xlist", "Olist", "progress", "quiet", NULL};

  if ( !PyArg_ParseTupleAndKeywords(args, keywds, "OO|Oi:HFKhat", kwlist,
	      &py_Xlist, &py_Olist, &py_callback, &quiet) )
    return NULL;

  if ( !PySequence_Check(py_Xlist) || !PySequence_Check(py_Olist) ){
    PyErr_SetString( ErrorObject, 
     "Arguments Xlist and Olist must support the sequence protocol.");
    return NULL;
  }

  gridsize = PySequence_Length(py_Xlist);
  if ( gridsize != PySequence_Length(py_Olist) ){
    PyErr_SetString( ErrorObject, 
     "The Xlist and Olist sequences must have the same length.");
    return NULL;
  }
    
  if ( gridsize > 16 ) {
    PyErr_SetString( ErrorObject,
      "Get real!  You will never get an answer with more than 16 arcs.");
    return NULL;
  }

  {
    int Xlist[gridsize];
    int Olist[gridsize];
    int i, m, a, failed = 0;
    PyObject *item;

    for (i=0; i<gridsize; i++) {
      item = PySequence_GetItem(py_Xlist, i);
      if (item == NULL || !PyInt_Check(item)) {
	failed = 1;
	break;
      }
      Xlist[i] = (int)PyInt_AS_LONG(item);
      Py_DECREF(item);
      item = PySequence_GetItem(py_Olist, i);
      if (item == NULL || !PyInt_Check(item)) {
	failed = 1;
	break;
      }
      Olist[i] = (int)PyInt_AS_LONG(item);
      Py_DECREF(item);
    }
    if (failed || !ValidGrid(gridsize, Xlist, Olist) ) {
      PyErr_SetString( ErrorObject,
      "Xlist and Olist must be permutations of {0,...,N}.");
    return NULL;
    }

    if (py_callback) {
      callback = py_callback;
      Py_XINCREF(py_callback);
    }
    else
      callback = NULL;

    Link link(gridsize, Xlist, Olist, &Progress, quiet);

    if (py_callback) {
      Py_XDECREF(py_callback);
    }

    result = PyDict_New();
    if ( result && !link.aborted ) {
      for(a = -link.HFK_Asize; a <= link.HFK_Asize; a++) {
	for(m = -link.HFK_Msize - gridsize; m <= link.HFK_Msize + gridsize; m++) {
	  if ( link.HFK_Rank(m,a) != 0 ) {
	    PyDict_SetItem(result, 
			   Py_BuildValue("(ii)", m, a),
			   PyInt_FromLong(link.HFK_Rank(m,a)));
	  }
	}
      }
    }
  }
  return result;
}

/* Documentation */
static char HFK_doc[]=
"HFKhat(Xlist, Olist[, progress], quiet=1)\n\n"
"The Xlist and Olist arguments are lists representing permutations of\n"
"{1,..,N} that give the columns of the X's and O's in a grid diagram\n"
"for a link.\n\n"
"The return value is a dictionary whose keys are tuples (m,a) such that\n"
"the mod 2 Heegaard Floer Homology group of the link, in Alexander\n"
"grading a and Maslov grading m has non-zero rank.  The corresponding\n"
"value in the dictionary is the rank of that group.\n\n"
"The optional progress argument is a python callback that will be called\n"
"from time to time with two arguments: a message string and an integer\n"
"between -1 and 100.  The message describes the current computational\n"
"task and the integer represents the percentage of the task that has\n"
"been completed.  A percent value of -1 indicates that the task has\n"
"finished.  The python function progress(msg, percent) should return\n"
"an integer value of 0 under normal circumstances.  A non-zero return\n"
"value aborts the computation.\n\n"
"If the function is called with quiet=0 then additional messages are\n"
"written to the standard output during the computation.\n\n"
;

/* List of functions defined in the module */

static PyMethodDef _hfk_methods[] = {
  {"HFKhat", (PyCFunction)HFKhat, METH_VARARGS|METH_KEYWORDS, HFK_doc},
  {NULL,		NULL}		/* sentinel */
};


/* Initialization function for the module */

DL_EXPORT(void)
init_hfk(void)
{
	PyObject *m, *d;

	/* Create the module and add the functions */
	m = Py_InitModule("_hfk", _hfk_methods);

	/* Add some symbolic constants to the module */
	d = PyModule_GetDict(m);
	ErrorObject = PyErr_NewException("hfk.error", NULL, NULL);
	PyDict_SetItemString(d, "HFK_Error", ErrorObject);
}
