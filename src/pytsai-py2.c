/**
 * pytsai.c
 * Copyright (C) Jonathan Merritt 2004.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "Python.h"
#include "tsai/cal_main.h"
#include "errors.h"

/*************************************
 * Forward Declarations of Functions *
 *************************************/
PyMODINIT_FUNC initpytsai(void);
static double* parse_calibration_data(PyObject *pyobj, int *size);
static int parse_camera_mapping(PyObject *obj);
static PyObject* build_camera_mapping();
static PyObject* tsai_coplanar_calibration(PyObject *self, PyObject *args);
static PyObject* tsai_noncoplanar_calibration(PyObject *self, PyObject *args);
static PyObject* tsai_coplanar_calibration_fo(PyObject *self, PyObject *args);
static PyObject* tsai_noncoplanar_calibration_fo(PyObject *self, 
        PyObject *args);
static PyObject* tsai_wc2ic(PyObject *self, PyObject *args);
static PyObject* tsai_ic2wc(PyObject *self, PyObject *args);
static PyObject* tsai_cc2wc(PyObject *self, PyObject *args);
static PyObject* tsai_add_sensor_coord_distortion(PyObject *self, 
        PyObject *args);

/***********************
 * Module Method Table *
 ***********************/
static PyMethodDef TsaiMethods[] = 
{

        {"_pytsai_coplanar_calibration", tsai_coplanar_calibration,
         METH_VARARGS, "Low level coplanar calibration routine."},

        {"_pytsai_noncoplanar_calibration", tsai_noncoplanar_calibration,
         METH_VARARGS, "Low level non-coplanar calibration routine."},

        {"_pytsai_coplanar_calibration_fo", tsai_coplanar_calibration_fo,
         METH_VARARGS, 
         "Low level coplanar, with full optimization calibration routine."},

        {"_pytsai_noncoplanar_calibration_fo", tsai_noncoplanar_calibration_fo,
         METH_VARARGS,
         "Low level noncoplanar, with full optimization calibration routine."},

        {"_pytsai_wc2ic", tsai_wc2ic, METH_VARARGS,
         "Low level conversion of world coordinates to image coordinates."},

        {"_pytsai_ic2wc", tsai_ic2wc, METH_VARARGS,
         "Low level conversion of image coordinates to world coordinates."},

        {"_pytsai_cc2wc", tsai_cc2wc, METH_VARARGS,
         "Low level conversion of camera coordinates to world coordinates."},

        {"_pytsai_add_sensor_coord_distortion",
         tsai_add_sensor_coord_distortion, METH_VARARGS,
         "Low level conversion of undistorted to distorted image coordinates."},

        {NULL, NULL, 0, NULL}
        
};

/****************************
 * Function Implementations *
 ****************************/

PyMODINIT_FUNC initpytsai(void)
{
    (void) Py_InitModule("pytsai", TsaiMethods);
}

/**
 * Parses calibration data.  The calibration data should be in the form of a
 * sequence of sequences. eg:
 *    [
 *        [ xs, ys, zs, xi, yi ], ...
 *    ]
 * where (xs, ys, zs) are the coordinates of a point in 3D space, and (xi, yi)
 * are the coordinates of the corresponding point in an image.
 *
 * If the method fails, it returns NULL and raises an appropriate exception.  
 * On success, a flat array of doubles will be PyMem_Malloc()'d and filled 
 * with the calibration data:
 *    { xs1, ys1, zs1, xi1, yi1, xs2, ys2, zs2, xi2, yi2, ... }
 * size will be populated with the number of calibration points.
 */
static double* parse_calibration_data(PyObject *pyobj, int *size)
{
        PyObject *subseq = NULL, *sstuple = NULL;
        int i = 0, ncoords = 0, index = 0;
        double *array = NULL, xs, ys, zs, xi, yi;

        /* check that pyobj is a sequence, and find out its size */
        if (PySequence_Check(pyobj) == 0)
        {
                PyErr_SetString(PyExc_TypeError, 
                        "First argument must be a sequence of coordinate " \
                        "sequences.");
                return NULL;
        }
        ncoords = PySequence_Size(pyobj);
        *size = ncoords;

        /* allocate memory */
        array = PyMem_Malloc(sizeof(double) * ncoords * 5);
        if (array == NULL)
                return NULL;

        /* iterate over each sub-sequence within pyobj, performing appropriate
         * checks and fetching data. */
        for (i = 0; i < ncoords; i++)
        {
                subseq = PySequence_GetItem(pyobj, i);
                if (PySequence_Check(subseq) == 0)
                {
                        Py_DECREF(subseq);
                        PyMem_Free(array);
                        PyErr_SetString(PyExc_TypeError,
                                "First argument must be a sequence of " \
                                "coordinate sequences.");
                        return NULL;
                }
                sstuple = PySequence_Tuple(subseq);
                Py_DECREF(subseq);
                if (sstuple == NULL)
                {
                        Py_DECREF(sstuple);
                        PyMem_Free(array);
                        PyErr_SetString(PyExc_TypeError,
                                "First argument must be a sequence of " \
                                "coordinate sequences.");
                        return NULL;
                }
                
                if (!PyArg_ParseTuple(sstuple, "ddddd", &xs, &ys, &zs,
                        &xi, &yi))
                {
                        Py_DECREF(sstuple);
                        PyMem_Free(array);
                        PyErr_SetString(PyExc_TypeError,
                                "First argument's coordinate sequences must " \
                                "contain 5 elements each: [x,y,z,xi,yi]");
                        return NULL;
                }
                
                array[index++] = xs;
                array[index++] = ys;
                array[index++] = zs;
                array[index++] = xi;
                array[index++] = yi;

                /* printf("xs = %f, ys = %f, zs = %f, xi = %f, yi = %f\n",
                        xs, ys, zs, xi, yi); */

                Py_DECREF(sstuple);
        }

        return array;
}

/**
 * Parses a mapping containing both camera constants and camera parameters.
 * The mapping should contain mappings whose keys are all strings.  The
 * following keys are recognized:
 *  Ncx
 *  Nfx
 *  dx
 *  dy
 *  dpx
 *  dpy
 *  Cx
 *  Cy
 *  sx
 *  f
 *  kappa1
 *  p1
 *  p2
 *  Tx
 *  Ty
 *  Tz
 *  Rx
 *  Ry
 *  Rz
 *  r1
 *  r2
 *  r3
 *  r4
 *  r5
 *  r6
 *  r7
 *  r8
 *  r9
 * 
 * The method returns 1 on success and 0 on failure.  If the method fails, an
 * appropriate exception is raised.
 */
#define TSAI_PARSE_ITEM(strx,name) { \
        mo = PyMapping_GetItemString(obj, #name); \
        if (mo != NULL) \
        { \
                number = PyNumber_Float(mo); \
                Py_DECREF(mo); \
                if (number == NULL) \
                { \
                        PyErr_SetString(PyExc_TypeError, \
                                "Camera parameter \"" #name \
                                "\" should be a number."); \
                        return 0; \
                } \
                strx.name = PyFloat_AsDouble(number); \
                Py_DECREF(number); \
        } \
}
static int parse_camera_mapping(PyObject *obj)
{
        PyObject *mo = NULL;
        PyObject *number = NULL;

        /* check that we have been passed a mapping */
        if (PyMapping_Check(obj) == 0)
        {
                PyErr_SetString(PyExc_TypeError,
                        "Second argument must be a mapping for camera " \
                        "parameters.");
                return 0;
        }
        
        /* parse all known items */
        TSAI_PARSE_ITEM(tsai_cp,Ncx);
        TSAI_PARSE_ITEM(tsai_cp,Nfx);
        TSAI_PARSE_ITEM(tsai_cp,dx);
        TSAI_PARSE_ITEM(tsai_cp,dy);
        TSAI_PARSE_ITEM(tsai_cp,dpx);
        TSAI_PARSE_ITEM(tsai_cp,dpy);
        TSAI_PARSE_ITEM(tsai_cp,Cx);
        TSAI_PARSE_ITEM(tsai_cp,Cy);
        TSAI_PARSE_ITEM(tsai_cp,sx);
        TSAI_PARSE_ITEM(tsai_cc,f);
        TSAI_PARSE_ITEM(tsai_cc,kappa1);
        TSAI_PARSE_ITEM(tsai_cc,p1);
        TSAI_PARSE_ITEM(tsai_cc,p2);
        TSAI_PARSE_ITEM(tsai_cc,Tx);
        TSAI_PARSE_ITEM(tsai_cc,Ty);
        TSAI_PARSE_ITEM(tsai_cc,Tz);
        TSAI_PARSE_ITEM(tsai_cc,Rx);
        TSAI_PARSE_ITEM(tsai_cc,Ry);
        TSAI_PARSE_ITEM(tsai_cc,Rz);
        TSAI_PARSE_ITEM(tsai_cc,r1);
        TSAI_PARSE_ITEM(tsai_cc,r2);
        TSAI_PARSE_ITEM(tsai_cc,r3);
        TSAI_PARSE_ITEM(tsai_cc,r4);
        TSAI_PARSE_ITEM(tsai_cc,r5);
        TSAI_PARSE_ITEM(tsai_cc,r6);
        TSAI_PARSE_ITEM(tsai_cc,r7);
        TSAI_PARSE_ITEM(tsai_cc,r8);
        TSAI_PARSE_ITEM(tsai_cc,r9);

        /* clear any exceptions that may have occurred while trying to fetch
         * keys */
        PyErr_Clear();

        /* return success */
        return 1;
}
#undef TSAI_PARSE_ITEM

/**
 * Constructs a mapping containing all camera parameters.  For parameters that
 * are known, see the parse_camera_mapping() function.
 */
static PyObject* build_camera_mapping()
{
        /* 28 parameters */
        return Py_BuildValue(
                "{sdsdsdsdsdsdsdsdsdsd" \
                 "sdsdsdsdsdsdsdsdsdsd" \
                 "sdsdsdsdsdsdsdsd}",
                "Ncx", tsai_cp.Ncx, "Nfx", tsai_cp.Nfx, 
                "dx", tsai_cp.dx, "dy", tsai_cp.dy, 
                "dpx", tsai_cp.dpx, "dpy", tsai_cp.dpy,
                "Cx", tsai_cp.Cx, "Cy", tsai_cp.Cy, 
                "sx", tsai_cp.sx, 
                "f", tsai_cc.f, 
                "kappa1", tsai_cc.kappa1, 
                "p1", tsai_cc.p1, "p2", tsai_cc.p2, 
                "Tx", tsai_cc.Tx, "Ty", tsai_cc.Ty, "Tz", tsai_cc.Tz, 
                "Rx", tsai_cc.Rx, "Ry", tsai_cc.Ry, "Rz", tsai_cc.Rz, 
                "r1", tsai_cc.r1, "r2", tsai_cc.r2, "r3", tsai_cc.r3, 
                "r4", tsai_cc.r4, "r5", tsai_cc.r5, "r6", tsai_cc.r6, 
                "r7", tsai_cc.r7, "r8", tsai_cc.r8, "r9", tsai_cc.r9);
}

/**
 * Performs coplanar calibration with optimization of f, Tz and kappa1.
 * The arguments to the function are:
 *      1 - set of calibration coordinates.
 *      2 - dictionary of camera parameters.
 */
static PyObject* tsai_coplanar_calibration(PyObject *self, PyObject *args)
{
        int i, index;
        PyObject *calibration_data = NULL;
        PyObject *params = NULL;
        double *calibration_array = NULL;
        int ncalibration_coords = 0;

        /* clear any error flags */
        pytsai_clear();

        if (!PyArg_ParseTuple(args, "OO", &calibration_data, &params))
                return NULL;
                
        /* fetch the calibration data */
        calibration_array = parse_calibration_data(calibration_data,
                &ncalibration_coords);
        if (calibration_array == NULL)
                return NULL;
        index = 0;
        tsai_cd.point_count = ncalibration_coords;
        for (i = 0; i < ncalibration_coords; i++)
        {
                tsai_cd.xw[i] = calibration_array[index++];
                tsai_cd.yw[i] = calibration_array[index++];
                tsai_cd.zw[i] = calibration_array[index++];
                tsai_cd.Xf[i] = calibration_array[index++];
                tsai_cd.Yf[i] = calibration_array[index++];
        }
        PyMem_Free(calibration_array);

        /* fetch the camera parameters mapping */
        if (parse_camera_mapping(params) == 0)
                return NULL;

        /* perform the C call */
        coplanar_calibration();

        /* check for an error */
        if (pytsai_haserror()) {
                PyErr_SetString(PyExc_RuntimeError, pytsai_string);
                return NULL;
        } else {
                return build_camera_mapping();
        }
}


/**
 * Performs non-coplanar calibration with optimization of f, Tz and kappa1.
 * The arguments to the function are:
 *      1 - set of calibration coordinates.
 *      2 - dictionary of camera parameters.
 */
static PyObject* tsai_noncoplanar_calibration(PyObject *self, PyObject *args)
{
        int i, index;
        PyObject *calibration_data = NULL;
        PyObject *params = NULL;
        double *calibration_array = NULL;
        int ncalibration_coords = 0;

        /* clear any error flags */
        pytsai_clear();

        if (!PyArg_ParseTuple(args, "OO", &calibration_data, &params))
                return NULL;
                
        /* fetch the calibration data */
        calibration_array = parse_calibration_data(calibration_data,
                &ncalibration_coords);
        if (calibration_array == NULL)
                return NULL;
        index = 0;
        tsai_cd.point_count = ncalibration_coords;
        for (i = 0; i < ncalibration_coords; i++)
        {
                tsai_cd.xw[i] = calibration_array[index++];
                tsai_cd.yw[i] = calibration_array[index++];
                tsai_cd.zw[i] = calibration_array[index++];
                tsai_cd.Xf[i] = calibration_array[index++];
                tsai_cd.Yf[i] = calibration_array[index++];
        }
        PyMem_Free(calibration_array);

        /* fetch the camera parameters mapping */
        if (parse_camera_mapping(params) == 0)
                return NULL;

        /* perform the C call */
        noncoplanar_calibration();

        /* check for an error */
        if (pytsai_haserror()) {
                PyErr_SetString(PyExc_RuntimeError, pytsai_string);
                return NULL;
        } else {
                return build_camera_mapping();
        }
}


/**
 * Performs coplanar calibration with full optimization.
 * The arguments to the function are:
 *      1 - set of calibration coordinates.
 *      2 - dictionary of camera parameters.
 */
static PyObject* tsai_coplanar_calibration_fo(PyObject *self, PyObject *args)
{
        int i, index;
        PyObject *calibration_data = NULL;
        PyObject *params = NULL;
        double *calibration_array = NULL;
        int ncalibration_coords = 0;

        /* clear any error flags */
        pytsai_clear();

        if (!PyArg_ParseTuple(args, "OO", &calibration_data, &params))
                return NULL;
                
        /* fetch the calibration data */
        calibration_array = parse_calibration_data(calibration_data,
                &ncalibration_coords);
        if (calibration_array == NULL)
                return NULL;
        index = 0;
        tsai_cd.point_count = ncalibration_coords;
        for (i = 0; i < ncalibration_coords; i++)
        {
                tsai_cd.xw[i] = calibration_array[index++];
                tsai_cd.yw[i] = calibration_array[index++];
                tsai_cd.zw[i] = calibration_array[index++];
                tsai_cd.Xf[i] = calibration_array[index++];
                tsai_cd.Yf[i] = calibration_array[index++];
        }
        PyMem_Free(calibration_array);

        /* fetch the camera parameters mapping */
        if (parse_camera_mapping(params) == 0)
                return NULL;

        /* perform the C call */
        coplanar_calibration_with_full_optimization();

        /* check for an error */
        if (pytsai_haserror()) {
                PyErr_SetString(PyExc_RuntimeError, pytsai_string);
                return NULL;
        } else {
                return build_camera_mapping();
        }
}


/**
 * Performs non-coplanar calibration with full optimization.
 * The arguments to the function are:
 *      1 - set of calibration coordinates.
 *      2 - dictionary of camera parameters.
 */
static PyObject* tsai_noncoplanar_calibration_fo(PyObject *self, 
        PyObject *args)
{
        int i, index;
        PyObject *calibration_data = NULL;
        PyObject *params = NULL;
        double *calibration_array = NULL;
        int ncalibration_coords = 0;

        /* clear any error flags */
        pytsai_clear();

        if (!PyArg_ParseTuple(args, "OO", &calibration_data, &params))
                return NULL;
                
        /* fetch the calibration data */
        calibration_array = parse_calibration_data(calibration_data,
                &ncalibration_coords);
        if (calibration_array == NULL)
                return NULL;
        index = 0;
        tsai_cd.point_count = ncalibration_coords;
        for (i = 0; i < ncalibration_coords; i++)
        {
                tsai_cd.xw[i] = calibration_array[index++];
                tsai_cd.yw[i] = calibration_array[index++];
                tsai_cd.zw[i] = calibration_array[index++];
                tsai_cd.Xf[i] = calibration_array[index++];
                tsai_cd.Yf[i] = calibration_array[index++];
        }
        PyMem_Free(calibration_array);

        /* fetch the camera parameters mapping */
        if (parse_camera_mapping(params) == 0)
                return NULL;

        /* perform the C call */
        noncoplanar_calibration_with_full_optimization();

        /* check for an error */
        if (pytsai_haserror()) {
                PyErr_SetString(PyExc_RuntimeError, pytsai_string);
                return NULL;
        } else {
                return build_camera_mapping();
        }
}


/**
 * Converts from world coordinates to image coordinates.
 * The arguments to the function are:
 *      1 - world coordinates (3-tuple)
 *      2 - dictionary of camera parameters
 * It returns a 2-tuple containing (Xf, Yf).
 */
static PyObject* tsai_wc2ic(PyObject *self, PyObject *args)
{
        double xw, yw, zw, Xf, Yf;
        PyObject *wc = NULL, *params = NULL, *wc2 = NULL;

        /* parse arguments */
        if (!PyArg_ParseTuple(args, "OO", &wc, &params))
                return NULL;
        if (!PyTuple_Check(wc))
        {
                wc2 = PySequence_Tuple(wc);
                if (wc2 == NULL)
                {
                        PyErr_SetString(PyExc_TypeError,
                                "First argument must be a 3-member sequence.");
                        return NULL;
                }
        }
        else
        {
                wc2 = wc;
                Py_INCREF(wc2);
        }
        if (!PyArg_ParseTuple(wc2, "ddd", &xw, &yw, &zw))
        {
                PyErr_SetString(PyExc_TypeError,
                        "First argument must be a 3-member sequence.");
                return NULL;
        }
        Py_DECREF(wc2);

        /* fetch the camera parameter mapping */
        if (parse_camera_mapping(params) == 0)
                return NULL;

        /* perform the C call */
        world_coord_to_image_coord(xw, yw, zw, &Xf, &Yf);

        /* return the value (Xf, Yf) */
        return Py_BuildValue("dd", Xf, Yf);
}


/**
 * Converts from image coordinates (+depth) to world coordinates.
 * The arguments to the function are:
 *      1 - image coordinates (3-tuple) (Xf, Yf, z)
 *      2 - dictionary of camera parameters
 * It returns a 2-tuple containing (xw, yw, zw).
 */
static PyObject* tsai_ic2wc(PyObject *self, PyObject *args)
{
        double xw, yw, zw, Xf, Yf;
        PyObject *wc = NULL, *params = NULL, *wc2 = NULL;

        /* parse arguments */
        if (!PyArg_ParseTuple(args, "OO", &wc, &params))
                return NULL;
        if (!PyTuple_Check(wc))
        {
                wc2 = PySequence_Tuple(wc);
                if (wc2 == NULL)
                {
                        PyErr_SetString(PyExc_TypeError,
                                "First argument must be a 3-member sequence.");
                        return NULL;
                }
        }
        else
        {
                wc2 = wc;
                Py_INCREF(wc2);
        }
        if (!PyArg_ParseTuple(wc2, "ddd", &Xf, &Yf, &zw))
        {
                PyErr_SetString(PyExc_TypeError,
                        "First argument must be a 3-member sequence.");
                return NULL;
        }
        Py_DECREF(wc2);

        /* fetch the camera parameter mapping */
        if (parse_camera_mapping(params) == 0)
                return NULL;

        /* perform the C call */
        image_coord_to_world_coord(Xf, Yf, zw, &xw, &yw);

        /* return the value (xw, yw, zw) */
        return Py_BuildValue("ddd", xw, yw, zw);
}


/**
 * Converts from camera coordinates to world coordinates.
 * The arguments to the function are:
 *      1 - camera coordinates (3-tuple) (xc, yc, zc)
 *      2 - dictionary of camera parameters
 * It returns a 2-tuple containing (xw, yw, zw).
 */
static PyObject* tsai_cc2wc(PyObject *self, PyObject *args)
{
        double xw, yw, zw, xc, yc, zc;
        PyObject *cc = NULL, *params = NULL, *cc2 = NULL;

        /* parse arguments */
        if (!PyArg_ParseTuple(args, "OO", &cc, &params))
                return NULL;
        if (!PyTuple_Check(cc))
        {
                cc2 = PySequence_Tuple(cc);
                if (cc2 == NULL)
                {
                        PyErr_SetString(PyExc_TypeError,
                                "First argument must be a 3-member sequence.");
                        return NULL;
                }
        }
        else
        {
                cc2 = cc;
                Py_INCREF(cc2);
        }
        if (!PyArg_ParseTuple(cc2, "ddd", &xc, &yc, &zc))
        {
                PyErr_SetString(PyExc_TypeError,
                        "First argument must be a 3-member sequence.");
                return NULL;
        }
        Py_DECREF(cc2);

        /* fetch the camera parameter mapping */
        if (parse_camera_mapping(params) == 0)
                return NULL;

        /* perform the C call */
        camera_coord_to_world_coord(xc, yc, zc, &xw, &yw, &zw);

        /* return the value (xw, yw, zw) */
        return Py_BuildValue("ddd", xw, yw, zw);
}


/**
 * Adds distortion to un-distorted sensor coordinates, returning distorted
 * sensor coordinates.  The arguments to this function are:
 *      1 - Xu
 *      2 - Yu
 *      3 - dictionary of camera parameters
 * It returns a 2-tuple containing (Xd, Yd).
 */
static PyObject* tsai_add_sensor_coord_distortion(PyObject *self, 
        PyObject *args)
{
        double Xu, Yu, Xd, Yd;
        PyObject *params = NULL;

        /* parse arguments */
        if (!PyArg_ParseTuple(args, "ddO", &Xu, &Yu, &params))
                return NULL;

        /* fetch the camera parameter mapping */
        if (parse_camera_mapping(params) == 0)
                return NULL;

        /* perform the C call */
        undistorted_to_distorted_sensor_coord(Xu, Yu, &Xd, &Yd);

        /* return the value (Xd, Yd) */
        return Py_BuildValue("dd", Xd, Yd);
}

