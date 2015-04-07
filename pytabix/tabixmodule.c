/*-
 * The MIT License
 *
 * Copyright (c) 2011 by Seoul National University.
 *               2014 by Kamil Slowikowski <slowikow@broadinstitute.org>
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*
 * Contact: Hyeshik Chang <hyeshik@snu.ac.kr>
 */

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "tabix.h"

static PyObject *TabixError;

typedef struct {
    PyObject_HEAD
    tabix_t *tb;
    char *fn;
} TabixObject;

typedef struct {
    PyObject_HEAD
    TabixObject *tbobj;
    ti_iter_t iter;
} TabixIteratorObject;

static PyTypeObject Tabix_Type, TabixIterator_Type;

/* --- TabixIterator --------------------------------------------------- */

static PyObject *
tabixiter_create(TabixObject *parentidx, ti_iter_t iter)
{
    TabixIteratorObject *self;

    self = PyObject_New(TabixIteratorObject, &TabixIterator_Type);
    if (self == NULL)
        return NULL;

    Py_INCREF(parentidx);
    self->tbobj = parentidx;
    self->iter = iter;

    return (PyObject *)self;
}

static void
tabixiter_dealloc(TabixIteratorObject *self)
{
    ti_iter_destroy(self->iter);
    Py_DECREF(self->tbobj);
    PyObject_Del(self);
}

static PyObject *
tabixiter_iter(PyObject *self)
{
    Py_INCREF(self);
    return self;
}

#if PY_MAJOR_VERSION < 3
# define PYOBJECT_FROM_STRING_AND_SIZE PyString_FromStringAndSize
#else
# define PYOBJECT_FROM_STRING_AND_SIZE PyUnicode_FromStringAndSize
#endif

static PyObject *
tabixiter_iternext(TabixIteratorObject *self)
{
    const char *chunk;
    int len, i;

    chunk = ti_read(self->tbobj->tb, self->iter, &len);
    if (chunk != NULL) {
        PyObject *ret, *column;
        Py_ssize_t colidx;
        const char *ptr, *begin;

        ret = PyList_New(0);
        if (ret == NULL)
            return NULL;

        colidx = 0;
        ptr = begin = chunk;
        for (i = len; i > 0; i--, ptr++)
            if (*ptr == '\t') {
                column = PYOBJECT_FROM_STRING_AND_SIZE(begin,
                                                       (Py_ssize_t)(ptr - begin));
                if (column == NULL || PyList_Append(ret, column) == -1) {
                    Py_DECREF(ret);
                    return NULL;
                }

                Py_DECREF(column);
                begin = ptr + 1;
                colidx++;
            }

        column = PYOBJECT_FROM_STRING_AND_SIZE(begin, (Py_ssize_t)(ptr - begin));
        if (column == NULL || PyList_Append(ret, column) == -1) {
            Py_DECREF(ret);
            return NULL;
        }
        Py_DECREF(column);

        return ret;
    }
    else
        return NULL;
}

static PyMethodDef tabixiter_methods[] = {
    {NULL, NULL} /* sentinel */
};

static PyTypeObject TabixIterator_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "tabix.iter",      /*tp_name*/
    sizeof(TabixIteratorObject), /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    /* methods */
    (destructor)tabixiter_dealloc,  /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    0,                          /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash*/
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,         /*tp_flags*/
    "An iterator for records returned by tabix.query().\n\n"
    "    >>> tb = tabix.open(\"http://localhost/file.bgz\")\n"
    "    >>> records = list(tb.query(\"1\", 1000, 2000))\n",                          /*tp_doc*/
    0,                          /*tp_traverse*/
    0,                          /*tp_clear*/
    0,                          /*tp_richcompare*/
    0,                          /*tp_weaklistoffset*/
    tabixiter_iter,             /*tp_iter*/
    (iternextfunc)tabixiter_iternext, /*tp_iternext*/
    tabixiter_methods,          /*tp_methods*/
    0,                          /*tp_members*/
    0,                          /*tp_getset*/
    0,                          /*tp_base*/
    0,                          /*tp_dict*/
    0,                          /*tp_descr_get*/
    0,                          /*tp_descr_set*/
    0,                          /*tp_dictoffset*/
    0,                          /*tp_init*/
    0,                          /*tp_alloc*/
    0,                          /*tp_new*/
    0,                          /*tp_free*/
    0,                          /*tp_is_gc*/
};


/* --- Tabix ----------------------------------------------------------- */

static PyObject *
tabix_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    TabixObject *self;
    const char *fn, *fnidx=NULL;
    static char *kwnames[]={"fn", "fnidx", NULL};
    tabix_t *tb;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|z:open",
                                     kwnames, &fn, &fnidx))
        return NULL;

    tb = ti_open(fn, fnidx);
    if (tb == NULL) {
        PyErr_SetString(TabixError, "Can't open the index file.");
        return NULL;
    }

    self = (TabixObject *)type->tp_alloc(type, 0);
    if (self == NULL)
        return NULL;

    self->tb = tb;
    self->fn = strdup(fn);

    return (PyObject *)self;
}

static void
tabix_dealloc(TabixObject *self)
{
    free(self->fn);
    ti_close(self->tb);
    PyObject_Del(self);
}

static PyObject *
tabix_query(TabixObject *self, PyObject *args)
{
    char *name;
    int begin, end;
    ti_iter_t result;

    if (!PyArg_ParseTuple(args, "sii:query", &name, &begin, &end))
        return NULL;

    result = ti_query(self->tb, name, begin - 1, end);
    if (result == NULL) {
        PyErr_SetString(TabixError, "query failed");
        return NULL;
    }

    return tabixiter_create(self, result);
}

static PyObject *
tabix_queryi(TabixObject *self, PyObject *args)
{
    int tid, begin, end;
    ti_iter_t result;

    if (!PyArg_ParseTuple(args, "iii:queryi", &tid, &begin, &end))
        return NULL;

    result = ti_queryi(self->tb, tid, begin - 1, end);
    if (result == NULL) {
        PyErr_SetString(TabixError, "query failed");
        return NULL;
    }

    return tabixiter_create(self, result);
}

static PyObject *
tabix_querys(TabixObject *self, PyObject *args)
{
    const char *reg;
    ti_iter_t result;

    if (!PyArg_ParseTuple(args, "s:querys", &reg))
        return NULL;

    result = ti_querys(self->tb, reg);
    if (result == NULL) {
        PyErr_SetString(TabixError, "query failed");
        return NULL;
    }

    return tabixiter_create(self, result);
}

static PyObject *
tabix_repr(TabixObject *self)
{
#if PY_MAJOR_VERSION < 3
    return PyString_FromFormat("<tabix fn=\"%s\">", self->fn);
#else
    return PyUnicode_FromFormat("<tabix fn=\"%s\">", self->fn);
#endif
}

static PyMethodDef tabix_methods[] = {
    {
        "query",
        (PyCFunction)tabix_query,
        METH_VARARGS,
        PyDoc_STR("Retrieve items within a region.\n\n"
                  "    >>> tb.query(\"chr1\", 1000, 2000)\n"
                  "    <tabix.iter at 0x17b86e50>\n\n"
                  "Parameters\n"
                  "----------\n"
                  "name : str\n"
                  "    Name of the sequence in the file.\n"
                  "start : int\n"
                  "    Start of the query region.\n"
                  "end : int\n"
                  "    End of the query region.\n")
    },
    {
        "queryi",
        (PyCFunction)tabix_queryi,
        METH_VARARGS,
        PyDoc_STR("Retrieve items within a region.\n\n"
                  "    >>> tb.queryi(0, 1000, 2000)\n"
                  "    <tabix.iter at 0x17b86e50>\n\n"
                  "Parameters\n"
                  "----------\n"
                  "tid : int\n"
                  "    Index of the sequence in the file (first is 0).\n"
                  "start : int\n"
                  "    Start of the query region.\n"
                  "end : int\n"
                  "    End of the query region.\n")
    },
    {
        "querys",
        (PyCFunction)tabix_querys,
        METH_VARARGS,
        PyDoc_STR("Retrieve items within a region.\n\n"
                  "    >>> tb.querys(\"chr1:1000-2000\")\n"
                  "    <tabix.iter at 0x17b86e50>\n\n"
                  "Parameters\n"
                  "----------\n"
                  "region : str\n"
                  "    Query string like \"seq:start-end\".\n")
    },
    /*
    {
        "header",
        (PyCFunction)tabix_header,
        METH_VARARGS,
        PyDoc_STR("Get the header for a file (VCF, SAM, GTF, etc.).\n\n"
                  "    >>> tb.header()\n")
    },
    */
    {NULL, NULL}           /* sentinel */
};
/*
	ti_iter_t ti_query(tabix_t *t, const char *name, int beg, int end);
	ti_iter_t ti_queryi(tabix_t *t, int tid, int beg, int end);
	ti_iter_t ti_querys(tabix_t *t, const char *reg);
*/

static PyTypeObject Tabix_Type = {
    /* The ob_type field must be initialized in the module init function
     * to be portable to Windows without using C++. */
    PyVarObject_HEAD_INIT(NULL, 0)
    "tabix.open",              /*tp_name*/
    sizeof(TabixObject),        /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    /* methods */
    (destructor)tabix_dealloc,  /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    (reprfunc)tabix_repr,       /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash*/
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,         /*tp_flags*/
    "Open a local or remote file compressed with bgzip and indexed with tabix.\n\n"
    "    >>> tb = tabix.open(\"http://localhost/file.bgz\")\n"
    "    >>> tb.query(\"1\", 1000, 2000)\n"
    "    <tabix.iter at 0x17b86e50>",                          /*tp_doc*/
    0,                          /*tp_traverse*/
    0,                          /*tp_clear*/
    0,                          /*tp_richcompare*/
    0,                          /*tp_weaklistoffset*/
    0,                          /*tp_iter*/
    0,                          /*tp_iternext*/
    tabix_methods,              /*tp_methods*/
    0,                          /*tp_members*/
    0,                          /*tp_getset*/
    0,                          /*tp_base*/
    0,                          /*tp_dict*/
    0,                          /*tp_descr_get*/
    0,                          /*tp_descr_set*/
    0,                          /*tp_dictoffset*/
    0,                          /*tp_init*/
    0,                          /*tp_alloc*/
    (newfunc)tabix_new,         /*tp_new*/
    0,                          /*tp_free*/
    0,                          /*tp_is_gc*/
};
/* --------------------------------------------------------------------- */

static PyMethodDef tabix_functions[] = {
    {NULL, NULL} /* sentinel */
};

PyDoc_STRVAR(module_doc,
"This module allows fast random access to files compressed with ``bgzip`` and\n"
"indexed by ``tabix``. It includes a C extension with code from klib_.\n"
"\n"
"Genomics data is often in a table where each row corresponds to a genomic\n"
"region::\n"
"\n"
"    chrom  pos      snp\n"
"    1      1000760  rs75316104\n"
"    1      1000894  rs114006445\n"
"    1      1000910  rs79750022\n"
"    1      1001177  rs4970401\n"
"    1      1001256  rs78650406\n"
"\n"
"With ``tabix``, you can quickly retrieve all rows in a genomic region by\n"
"specifying a query with a sequence name, start, and end::\n"
"\n"
"    import tabix\n"
"\n"
"    # Open a remote or local file.\n"
"    url = \"ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/\"\n"
"    url += \"ALL.2of4intersection.20100804.genotypes.vcf.gz\"\n"
"\n"
"    tb = tabix.open(url)\n"
"\n"
"    # These queries are identical. A query returns an iterator over the\n"
"    # results.\n"
"    records = tb.query(\"1\", 1000000, 1250000)\n"
"\n"
"    records = tb.queryi(0, 1000000, 1250000)\n"
"\n"
"    records = tb.querys(\"1:1000000-1250000\")\n"
"\n"
"    # Each record is a list of strings.\n"
"    for record in records:\n"
"        print record[:5]\n"
"            break\n"
"\n"
"The ``bgzip`` and ``tabix`` programs are provided in samtools_.\n"
"\n"
".. _klib: https://github.com/jmarshall/klib\n"
".. _samtools: http://samtools.sourceforge.net\n");

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef tabixmodule = { 
    PyModuleDef_HEAD_INIT,
    "tabix",
    module_doc,
    -1, 
    tabix_functions,
    NULL,
    NULL,
    NULL,
    NULL
};
#endif

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC inittabix(void)
#else
PyMODINIT_FUNC PyInit_tabix(void)
#endif
{
    PyObject *m;

    if (PyType_Ready(&Tabix_Type) < 0)
        goto fail;
    if (PyType_Ready(&TabixIterator_Type) < 0)
        goto fail;

#if PY_MAJOR_VERSION < 3
    m = Py_InitModule3("tabix", tabix_functions, module_doc);
#else
    m = PyModule_Create(&tabixmodule);
#endif
    if (m == NULL)
        goto fail;

    if (TabixError == NULL) {
        TabixError = PyErr_NewException("tabix.TabixError", NULL, NULL);
        if (TabixError == NULL)
            goto fail;
    }
    Py_INCREF(TabixError);
    PyModule_AddObject(m, "TabixError", TabixError);

    PyModule_AddObject(m, "open", (PyObject *)&Tabix_Type);
    PyModule_AddObject(m, "iter", (PyObject *)&TabixIterator_Type);

#if PY_MAJOR_VERSION >= 3
    return m;
#endif

 fail:
#if PY_MAJOR_VERSION < 3
    return;
#else
    return NULL;
#endif
}
