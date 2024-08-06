#include <Python.h>
#include <numpy/arrayobject.h>
/* #include <MMTK/core.h> */
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define MAXLINE 100

#define ARGERROR         "wrong arguments error"
#define FLOATERROR       "float number convertion error"
#define INTERROR         "integer number convertion error"
#define CHRERROR         "character string convertion error"
#define UNKNOWNTYPE      "unsupported data type error"
#define LISTAPPENDERROR  "list append failed"
#define READERROR        "file reading error"

static PyObject *
dataBlock(PyObject *self, PyObject *args)
{
  PyObject *_pyfile, *_record;
  PyObject *data, *item, *record, *fieldspec;
  int linelen, i, rsize, r;
  long numrecords, _numrecords, fieldwidth;
  char buf[MAXLINE], field[MAXLINE], *fieldtype;
  FILE *textfile;

  if (!PyArg_ParseTuple(args, "O!O!l",
			&PyFile_Type,&_pyfile,
                        &PyTuple_Type,&_record,
                        &_numrecords))
      return Py_BuildValue("s", ARGERROR);

  rsize = (int)PyTuple_Size(_record);

  textfile = PyFile_AsFile(_pyfile);
  data = (PyObject *)PyList_New(0);

  numrecords = 0;
  i = 0;
  linelen = 0;
  while (numrecords < _numrecords) {
    if (rsize > 1) { record = (PyObject *)PyTuple_New(rsize); }
    for (r=0;r<rsize;r++) {
      fieldspec = (PyObject *)PyTuple_GetItem(_record,r);
      fieldtype = (char *)PyString_AsString
	((PyObject *)PyTuple_GetItem(fieldspec,0));
      fieldwidth = (long)PyInt_AsLong
	((PyObject *)PyTuple_GetItem(fieldspec,1));
      if (i > (linelen - fieldwidth)) {
	if (fgets(buf,MAXLINE,textfile) == NULL)
	  return Py_BuildValue("s", READERROR);
	linelen = strlen(buf);
	i = 0;
      }
      strncpy(field,buf+i,fieldwidth);
      i += fieldwidth;
      switch(*fieldtype) {
      case 'F':
	item = PyFloat_FromDouble(atof(field));
	if (rsize > 1) {
	  Py_INCREF(item); /* SetItem steals reference... */
	  if (PyTuple_SetItem((PyObject *)record,r,item) < 0)
	    return Py_BuildValue("s", FLOATERROR); }
	else { 	
	  if (PyList_Append((PyObject *)data,item) < 0)
	    return Py_BuildValue("s", FLOATERROR); }
	Py_DECREF(item);
	break;
      case 'I':
	/* item = PyLong_FromLong(atoi(field)); */
	item = PyInt_FromLong(atoi(field));
	if (rsize > 1) {
	  Py_INCREF(item); /* SetItem steals reference... */
	  if (PyTuple_SetItem((PyObject *)record,r,item) < 0)
	    return Py_BuildValue("s", INTERROR); }
	else { 	
	  if (PyList_Append((PyObject *)data,item) < 0)
	    return Py_BuildValue("s", INTERROR); }
	Py_DECREF(item);
	break;
      case 'A':
	item = PyString_FromString(field);
	if (rsize > 1) {
	  Py_INCREF(item); /* SetItem steals reference... */
	  if (PyTuple_SetItem((PyObject *)record,r,item) < 0)
	    return Py_BuildValue("s", CHRERROR); }
	else {
	  if (PyList_Append((PyObject *)data,item) < 0)
	    return Py_BuildValue("s", CHRERROR); }
	Py_DECREF(item);
	break;
      default:
	return Py_BuildValue("s",UNKNOWNTYPE);
      }
    }
    if (rsize > 1) {
      if (PyList_Append((PyObject *)data,record) < 0)
	return Py_BuildValue("s", LISTAPPENDERROR);
      Py_DECREF(record);
    }
    numrecords++;
  }
  /* Py_INCREF(data); */
  return data;
}

static PyMethodDef Methods[] = {
  {"dataBlock",dataBlock,METH_VARARGS},
  {NULL,NULL}
};

void initmdconv_caux()
{
  PyObject *module;
  module = Py_InitModule("caux",Methods);
#ifdef import_array
  import_array();
#endif
  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module 'caux'");
}
