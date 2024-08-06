#include <Python.h>
/* #include <Numeric/arrayobject.h> */
#include <numpy/arrayobject.h>
#include <stdio.h>
#include "neighbors.h"

typedef struct {
  PyObject_HEAD
  grid *g ;
} py_grid;

static void py_grid_dealloc(py_grid *self) {
  free_grid(self->g);
  self->ob_type->tp_free((PyObject *) self);
}

static PyObject *py_grid_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  py_grid *self;

  self = (py_grid *) type->tp_alloc(type, 0);
  if(self != NULL)
    self->g = new_grid();

  return (PyObject *) self;
}

static PyObject *py_grid_build_grid(py_grid *self, PyObject *args) {
  PyArrayObject *abox;
  PyArrayObject *acoords;
  REAL rcut;
  rvect box[3];
  rvect *coords;
  int atomnum;
  int i, j;

#ifdef DOUBLE
  if(!PyArg_ParseTuple(args, "O!O!d", &PyArray_Type, &abox, 
                                      &PyArray_Type, &acoords, &rcut))
    return NULL;

  if((abox->descr->type_num != PyArray_DOUBLE) ||
     (acoords->descr->type_num != PyArray_DOUBLE)) {
    PyErr_SetString(PyExc_TypeError, "module was compiled for double precision but received other type of matrices");
    return NULL;
  }
#else
  if(!PyArg_ParseTuple(args, "O!O!f", &PyArray_Type, &abox, 
                                      &PyArray_Type, &acoords, &rcut))
    return NULL;

  if((abox->descr->type_num != PyArray_FLOAT) ||
     (acoords->descr->type_num != PyArray_FLOAT)) {
    PyErr_SetString(PyExc_TypeError, "module was compiled for single precision but received other type of matrices");
    return NULL;
  }
#endif

  atomnum = acoords->dimensions[0];

  coords = (rvect *) malloc(sizeof(rvect)*atomnum);
  for(i=0; i<atomnum; i++)
    for(j=0; j<3; j++)
#ifdef DOUBLE
      coords[i][j] = *(double *)(acoords->data+i*acoords->strides[0]+j*acoords->strides[1]);
#else
      coords[i][j] = *(float *)(acoords->data+i*acoords->strides[0]+j*acoords->strides[1]);
#endif

  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
#ifdef DOUBLE
      box[i][j] = *(double *)(abox->data+i*abox->strides[0]+j*abox->strides[1]);
#else
      box[i][j] = *(float *)(abox->data+i*abox->strides[0]+j*abox->strides[1]);
#endif

  build_grid(box, rcut, atomnum, coords, self->g);

  free(coords);

  Py_XINCREF(Py_None);
  return Py_None;
}

static PyObject *py_grid_find_neighbors(py_grid* self, PyObject *args) {
  int atom;
  atomlist neighbors;
  npy_intp dimensions[2];
  PyArrayObject *aneighbors, *aimages, *adists;
  int i, j;

  if(!PyArg_ParseTuple(args, "i", &atom))
    return NULL;

  neighbors = find_neighbors(atom, self->g);
  if(neighbors.n == 0) {
    Py_XINCREF(Py_None);
    return Py_None;
  }

  dimensions[0] = neighbors.n;
  dimensions[1] = 3;

  aneighbors = (PyArrayObject *) PyArray_SimpleNew(1, dimensions, PyArray_INT);
  aimages = (PyArrayObject *) PyArray_SimpleNew(2, dimensions, PyArray_INT);
  for(i=0; i<neighbors.n; i++)
  {
    *(int *)(aneighbors->data + i*aneighbors->strides[0]) = neighbors.idx[i];
    for(j=0; j<3; j++)
      *(int *)(aimages->data + i*aimages->strides[0] + j*aimages->strides[1]) = neighbors.images[i][j];
  }
#ifdef DOUBLE
  adists = (PyArrayObject *) PyArray_SimpleNew(1, dimensions, PyArray_DOUBLE);
  for(i=0; i<neighbors.n; i++)
    *(double *)(adists->data + i*adists->strides[0]) = neighbors.dists[i];
#else
  adists = (PyArrayObject *) PyArray_SimpleNew(1, dimensions, PyArray_FLOAT);
  for(i=0; i<neighbors.n; i++)
    *(float *)(adists->data + i*adists->strides[0]) = neighbors.dists[i];
#endif

  return Py_BuildValue("(NNN)", aneighbors, aimages, adists);
}

static PyMethodDef py_grid_methods[] = {
    {"build_grid", (PyCFunction) py_grid_build_grid, METH_VARARGS,
     "Construct a grid for given box, coordinates and rcut"
    },
    {"find_neighbors", (PyCFunction) py_grid_find_neighbors, METH_VARARGS,
     "Find neighbors for a given atom",
    },
    {NULL}  /* Sentinel */
};

static PyTypeObject py_gridType = {
  PyObject_HEAD_INIT(NULL)
  0,                                /*ob_size*/
  "py_grid",                        /*tp_name*/
  sizeof(py_grid),                  /*tp_basicsize*/
  0,                                /*tp_itemsize*/
  (destructor) py_grid_dealloc,     /*tp_dealloc*/
  0,                                /*tp_print*/
  0,                                /*tp_getattr*/
  0,                                /*tp_setattr*/
  0,                                /*tp_compare*/
  0,                                /*tp_repr*/
  0,                                /*tp_as_number*/
  0,                                /*tp_as_sequence*/
  0,                                /*tp_as_mapping*/
  0,                                /*tp_hash */
  0,                                /*tp_call*/
  0,                                /*tp_str*/
  0,                                /*tp_getattro*/
  0,                                /*tp_setattro*/
  0,                                /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT,               /*tp_flags*/
  "wrapped grid objects",           /* tp_doc */
  0,                                /* tp_traverse */
  0,                                /* tp_clear */
  0,                                /* tp_richcompare */
  0,                                /* tp_weaklistoffset */
  0,                                /* tp_iter */
  0,                                /* tp_iternext */
  py_grid_methods,                  /* tp_methods */
  0,                                /* tp_members */
  0,                                /* tp_getset */
  0,                                /* tp_base */
  0,                                /* tp_dict */
  0,                                /* tp_descr_get */
  0,                                /* tp_descr_set */
  0,                                /* tp_dictoffset */
  0,                                /* tp_init */
  0,                                /* tp_alloc */
  (newfunc) py_grid_new             /* tp_new */
};

static PyMethodDef gridwrapper_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC initgridwrapper(void) {
  PyObject *m;

  if(PyType_Ready(&py_gridType) < 0)
    return;

  m = Py_InitModule3("gridwrapper", gridwrapper_methods,
                     "Wrapper for the \"grid *\" type.");

  Py_INCREF(&py_gridType);
  PyModule_AddObject(m, "py_grid", (PyObject *) &py_gridType);

  import_array();
}

