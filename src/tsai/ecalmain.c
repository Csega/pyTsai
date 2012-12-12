/**
 * ecalmain.c
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

/****************************************************************************\
*                                                                            *
* This file contains routines for calibrating the extrinsic parameters of    *
* Tsai's 11 parameter camera model.  The inputs to the routines are a set of *
* precalibrated intrinsic camera parameters:                                 *
*                                                                            *
*       f       - effective focal length of the pin hole camera              *
*       kappa1  - 1st order radial lens distortion coefficient               *
*       Cx, Cy  - coordinates of center of radial lens distortion            *
*                 (also used as the piercing point of the camera coordinate  *
*                  frame's Z axis with the camera's sensor plane)            *
*       sx      - uncertainty factor for scale of horizontal scanline        *
*                                                                            *
* and a set of calibration data consisting of the (x,y,z) world coordinates  *
* of a feature point (in mm) and the corresponding coordinates (Xf,Yf) (in   *
* pixels) of the feature point in the image.  The outputs of the routines    *
* are the 6 external (also called extrinsic or exterior) camera parameters:  *
*                                                                            *
*       Rx, Ry, Rz, Tx, Ty, Tz  - rotational and translational components of *
*                                 the transform between the world's          *
*                                 coordinate frame and the camera's          *
*                                 coordinate frame.                          *
*                                                                            *
* describing the camera's pose.                                              *
*                                                                            *
* This file provides two routines:                                           *
*                                                                            *
*       coplanar_extrinsic_parameter_estimation ()                           *
* and                                                                        *
*       noncoplanar_extrinsic_parameter_estimation ()                        *
*                                                                            *
* which are used respectively for coplanar and non-coplanar calibration      *
* data.                                                                      *
*                                                                            *
* Initial estimates for the extrinsic camera parameters are determined using *
* a modification of the first stage of Tsai's algorithm.  These estimates    *
* are then refined using iterative non-linear optimization.                  *
*                                                                            *
*                                                                            *
* History                                                                    *
* -------                                                                    *
*                                                                            *
* 01-Nov-04  Jonathan Merritt (j.merritt@pgrad.unimelb.edu.au) at            *
*            The University of Melbourne.                                    *
*       o Released new changes under the LGPL.                               *
*       o Removed exit() calls.                                              *
*                                                                            *
* 15-Oct-95  Reg Willson (rgwillson@mmm.com) at 3M St. Paul, MN              *
*       Added routines to coplanar_extrinsic_parameter_estimation to pick    *
*       the correct rotation matrix solution from the two possible solutions.*
*       Bug tracked down by Pete Rander <Peter.Rander@IUS4.IUS.CS.CMU.EDU>.  *
*                                                                            *
* 20-May-95  Reg Willson (rgwillson@mmm.com) at 3M St. Paul, MN              *
*       Return the error to lmdif rather than the squared error.             *
*         lmdif calculates the squared error internally during optimization. *
*         Before this change calibration was essentially optimizing error^4. *
*                                                                            *
* 02-Apr-95  Reg Willson (rgwillson@mmm.com) at 3M St. Paul, MN              *
*       Rewrite memory allocation to avoid memory alignment problems         *
*       on some machines.                                                    *
*       Strip out IMSL code.  MINPACK seems to work fine.                    *
*       Filename changes for DOS port.                                       *
*                                                                            *
* 04-Jun-94  Reg Willson (rgwillson@mmm.com) at 3M St. Paul, MN              *
*       Added alternate macro definitions for less common math functions.    *
*                                                                            *
* 25-Mar-94  Torfi Thorhallsson (torfit@verk.hi.is) at the University of     *
*            Iceland                                                         *
*       Added a new version of the routine epe_optimize() which uses the     *
*       *public domain* MINPACK optimization library instead of IMSL.        *
*       To select the new routine, compile this file with the flag -DMINPACK *
*                                                                            *
* 11-Nov-93  Reg Willson (rgw@cs.cmu.edu) at Carnegie-Mellon University      *
*       Original implementation.                                             *
*                                                                            *
\****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include "../matrix/matrix.h"
#include "cal_main.h"
#include "../minpack/f2c.h"
#include "../minpack/minpack.h"
#include "../errors.h"


/***********************************************************************\
* Routines for coplanar extrinsic parameter estimation			*
\***********************************************************************/
/* pytsai: can fail; int return type is required. */
int cepe_compute_U (U)
    double    U[];
{
    dmat      M,
              a,
              b;

    double    Xd,
              Yd,
              Xu,
              Yu,
              distortion_factor;

    int       i;

    M = newdmat (0, (tsai_cd.point_count - 1), 0, 4, &errno);
    if (errno) {
	pytsai_raise("cepe compute U: unable to allocate matrix M");
	return 0;
    }

    a = newdmat (0, 4, 0, 0, &errno);
    if (errno) {
        freemat(M);
	pytsai_raise("cepe compute U: unable to allocate vector a");
	return 0;
    }

    b = newdmat (0, (tsai_cd.point_count - 1), 0, 0, &errno);
    if (errno) {
        freemat(M);
        freemat(a);
	pytsai_raise("cepe compute U: unable to allocate vector b");
	return 0;
    }

    for (i = 0; i < tsai_cd.point_count; i++) {
	/* convert from image coordinates to distorted sensor coordinates */
	Xd = tsai_cp.dpx * (tsai_cd.Xf[i] - tsai_cp.Cx) / tsai_cp.sx;
	Yd = tsai_cp.dpy * (tsai_cd.Yf[i] - tsai_cp.Cy);

	/* convert from distorted sensor coordinates to undistorted sensor coordinates */
	distortion_factor = 1 + tsai_cc.kappa1 * (SQR (Xd) + SQR (Yd));
	Xu = Xd * distortion_factor;
	Yu = Yd * distortion_factor;

	M.el[i][0] = Yu * tsai_cd.xw[i];
	M.el[i][1] = Yu * tsai_cd.yw[i];
	M.el[i][2] = Yu;
	M.el[i][3] = -Xu * tsai_cd.xw[i];
	M.el[i][4] = -Xu * tsai_cd.yw[i];
	b.el[i][0] = Xu;
    }

    if (solve_system (M, a, b)) {
        freemat(M);
        freemat(a);
        freemat(b);
	pytsai_raise("cepe compute U: unable to solve system  Ma=b");
	return 0;
    }

    U[0] = a.el[0][0];
    U[1] = a.el[1][0];
    U[2] = a.el[2][0];
    U[3] = a.el[3][0];
    U[4] = a.el[4][0];

    freemat (M);
    freemat (a);
    freemat (b);

    return 1;
}


/* pytsai: cannot fail; void return type is fine. */
void cepe_compute_Tx_and_Ty (U)
    double    U[];
{
    double    Tx,
              Ty,
              Ty_squared,
              x,
              y,
              Sr,
              r1p,
              r2p,
              r4p,
              r5p,
              r1,
              r2,
              r4,
              r5,
              distance,
              far_distance;

    int       i,
              far_point;

    r1p = U[0];
    r2p = U[1];
    r4p = U[3];
    r5p = U[4];

    /* first find the square of the magnitude of Ty */
    if ((fabs (r1p) < EPSILON) && (fabs (r2p) < EPSILON))
	Ty_squared = 1 / (SQR (r4p) + SQR (r5p));
    else if ((fabs (r4p) < EPSILON) && (fabs (r5p) < EPSILON))
	Ty_squared = 1 / (SQR (r1p) + SQR (r2p));
    else if ((fabs (r1p) < EPSILON) && (fabs (r4p) < EPSILON))
	Ty_squared = 1 / (SQR (r2p) + SQR (r5p));
    else if ((fabs (r2p) < EPSILON) && (fabs (r5p) < EPSILON))
	Ty_squared = 1 / (SQR (r1p) + SQR (r4p));
    else {
	Sr = SQR (r1p) + SQR (r2p) + SQR (r4p) + SQR (r5p);
	Ty_squared = (Sr - sqrt (SQR (Sr) - 4 * SQR (r1p * r5p - r4p * r2p))) /
	 (2 * SQR (r1p * r5p - r4p * r2p));
    }

    /* find a point that is far from the image center */
    far_distance = 0;
    far_point = 0;
    for (i = 0; i < tsai_cd.point_count; i++)
	if ((distance = SQR (tsai_cd.Xf[i] - tsai_cp.Cx) + SQR (tsai_cd.Yf[i] - tsai_cp.Cy)) > far_distance) {
	    far_point = i;
	    far_distance = distance;
	}

    /* now find the sign for Ty */
    /* start by assuming Ty > 0 */
    Ty = sqrt (Ty_squared);
    r1 = U[0] * Ty;
    r2 = U[1] * Ty;
    Tx = U[2] * Ty;
    r4 = U[3] * Ty;
    r5 = U[4] * Ty;
    x = r1 * tsai_cd.xw[far_point] + r2 * tsai_cd.yw[far_point] + Tx;
    y = r4 * tsai_cd.xw[far_point] + r5 * tsai_cd.yw[far_point] + Ty;

    /* flip Ty if we guessed wrong */
    if ((SIGNBIT (x) != SIGNBIT (tsai_cd.Xf[far_point] - tsai_cp.Cx)) ||
	(SIGNBIT (y) != SIGNBIT (tsai_cd.Yf[far_point] - tsai_cp.Cy)))
	Ty = -Ty;

    /* update the calibration constants */
    tsai_cc.Tx = U[2] * Ty;
    tsai_cc.Ty = Ty;
}


/* pytsai: cannot fail; void return type is fine. */
void cepe_compute_R (U)
    double    U[];
{
    double    r1,
              r2,
              r3,
              r4,
              r5,
              r6,
              r7,
              r8,
              r9;

    r1 = U[0] * tsai_cc.Ty;
    r2 = U[1] * tsai_cc.Ty;
    r3 = sqrt (1 - SQR (r1) - SQR (r2));

    r4 = U[3] * tsai_cc.Ty;
    r5 = U[4] * tsai_cc.Ty;
    r6 = sqrt (1 - SQR (r4) - SQR (r5));
    if (!SIGNBIT (r1 * r4 + r2 * r5))
	r6 = -r6;

    /* use the outer product of the first two rows to get the last row */
    r7 = r2 * r6 - r3 * r5;
    r8 = r3 * r4 - r1 * r6;
    r9 = r1 * r5 - r2 * r4;

    /* update the calibration constants */
    tsai_cc.r1 = r1;
    tsai_cc.r2 = r2;
    tsai_cc.r3 = r3;
    tsai_cc.r4 = r4;
    tsai_cc.r5 = r5;
    tsai_cc.r6 = r6;
    tsai_cc.r7 = r7;
    tsai_cc.r8 = r8;
    tsai_cc.r9 = r9;

    /* fill in tsai_cc.Rx, tsai_cc.Ry and tsai_cc.Rz */
    solve_RPY_transform ();
}


/* pytsai: can fail; int return type is required. */
int cepe_compute_approximate_f (f)
    double   *f;
{
    dmat      M,
              a,
              b;

    double    Yd;

    int       i;

    M = newdmat (0, (tsai_cd.point_count - 1), 0, 1, &errno);
    if (errno) {
	pytsai_raise("cepe compute apx: unable to allocate matrix M");
	return 0;
    }

    a = newdmat (0, 1, 0, 0, &errno);
    if (errno) {
        freemat(M);
	pytsai_raise("cepe compute apx: unable to allocate vector a");
	return 0;
    }

    b = newdmat (0, (tsai_cd.point_count - 1), 0, 0, &errno);
    if (errno) {
        freemat(M);
        freemat(a);
	pytsai_raise("cepe compute apx: unable to allocate vector b");
	return 0;
    }

    for (i = 0; i < tsai_cd.point_count; i++) {
	Yd = tsai_cp.dpy * (tsai_cd.Yf[i] - tsai_cp.Cy);

	M.el[i][0] = tsai_cc.r4 * tsai_cd.xw[i] + tsai_cc.r5 * tsai_cd.yw[i] + tsai_cc.Ty;
	M.el[i][1] = -Yd;
	b.el[i][0] = (tsai_cc.r7 * tsai_cd.xw[i] + tsai_cc.r8 * tsai_cd.yw[i]) * Yd;
    }

    if (solve_system (M, a, b)) {
        freemat(M);
        freemat(a);
        freemat(b);
	pytsai_raise("cepe compute apx: unable to solve system  Ma=b");
	return 0;
    }

    /* return the approximate effective focal length */
    *f = a.el[0][0];

    freemat (M);
    freemat (a);
    freemat (b);

    return 1;
}


/***********************************************************************\
* Routines for noncoplanar extrinsic parameter estimation		*
\***********************************************************************/
/* pytsai: can fail; int return type is required. */
int ncepe_compute_U (U)
    double    U[];
{
    dmat      M,
              a,
              b;

    double    Xu,
              Yu,
              Xd,
              Yd,
              distortion_factor;

    int       i;

    M = newdmat (0, (tsai_cd.point_count - 1), 0, 6, &errno);
    if (errno) {
	pytsai_raise("ncepe compute U: unable to allocate matrix M");
	return 0;
    }

    a = newdmat (0, 6, 0, 0, &errno);
    if (errno) {
        freemat(M);
	pytsai_raise("ncepe compute U: unable to allocate vector a");
	return 0;
    }

    b = newdmat (0, (tsai_cd.point_count - 1), 0, 0, &errno);
    if (errno) {
        freemat(M);
        freemat(a);
	pytsai_raise("ncepe compute U: unable to allocate vector b");
	return 0;
    }

    for (i = 0; i < tsai_cd.point_count; i++) {
	/* convert from image coordinates to distorted sensor coordinates */
	Xd = tsai_cp.dpx * (tsai_cd.Xf[i] - tsai_cp.Cx) / tsai_cp.sx;
	Yd = tsai_cp.dpy * (tsai_cd.Yf[i] - tsai_cp.Cy);

	/* convert from distorted sensor coordinates to undistorted sensor coordinates */
	distortion_factor = 1 + tsai_cc.kappa1 * (SQR (Xd) + SQR (Yd));
	Xu = Xd * distortion_factor;
	Yu = Yd * distortion_factor;

	M.el[i][0] = Yu * tsai_cd.xw[i];
	M.el[i][1] = Yu * tsai_cd.yw[i];
	M.el[i][2] = Yu * tsai_cd.zw[i];
	M.el[i][3] = Yu;
	M.el[i][4] = -Xu * tsai_cd.xw[i];
	M.el[i][5] = -Xu * tsai_cd.yw[i];
	M.el[i][6] = -Xu * tsai_cd.zw[i];
	b.el[i][0] = Xu;
    }

    if (solve_system (M, a, b)) {
        freemat(M);
        freemat(a);
        freemat(b);
	pytsai_raise("ncepe compute U: unable to solve system  Ma=b");
	return 0;
    }

    U[0] = a.el[0][0];
    U[1] = a.el[1][0];
    U[2] = a.el[2][0];
    U[3] = a.el[3][0];
    U[4] = a.el[4][0];
    U[5] = a.el[5][0];
    U[6] = a.el[6][0];

    freemat (M);
    freemat (a);
    freemat (b);

    return 1;
}


/* pytsai: cannot fail; void return type is fine. */
void ncepe_compute_Tx_and_Ty (U)
    double    U[];
{
    double    Tx,
              Ty,
              Ty_squared,
              x,
              y,
              r1,
              r2,
              r3,
              r4,
              r5,
              r6,
              distance,
              far_distance;

    int       i,
              far_point;

    /* first find the square of the magnitude of Ty */
    Ty_squared = 1 / (SQR (U[4]) + SQR (U[5]) + SQR (U[6]));

    /* find a point that is far from the image center */
    far_distance = 0;
    far_point = 0;
    for (i = 0; i < tsai_cd.point_count; i++)
	if ((distance = SQR (tsai_cd.Xf[i] - tsai_cp.Cx) + SQR (tsai_cd.Yf[i] - tsai_cp.Cy)) > far_distance) {
	    far_point = i;
	    far_distance = distance;
	}

    /* now find the sign for Ty */
    /* start by assuming Ty > 0 */
    Ty = sqrt (Ty_squared);
    r1 = U[0] * Ty;
    r2 = U[1] * Ty;
    r3 = U[2] * Ty;
    Tx = U[3] * Ty;
    r4 = U[4] * Ty;
    r5 = U[5] * Ty;
    r6 = U[6] * Ty;
    x = r1 * tsai_cd.xw[far_point] + r2 * tsai_cd.yw[far_point] + r3 * tsai_cd.zw[far_point] + Tx;
    y = r4 * tsai_cd.xw[far_point] + r5 * tsai_cd.yw[far_point] + r6 * tsai_cd.zw[far_point] + Ty;

    /* flip Ty if we guessed wrong */
    if ((SIGNBIT (x) != SIGNBIT (tsai_cd.Xf[far_point] - tsai_cp.Cx)) ||
	(SIGNBIT (y) != SIGNBIT (tsai_cd.Yf[far_point] - tsai_cp.Cy)))
	Ty = -Ty;

    /* update the calibration constants */
    tsai_cc.Tx = U[3] * Ty;
    tsai_cc.Ty = Ty;
}


/* pytsai: cannot fail; void return type is fine. */
void ncepe_compute_R (U)
    double    U[];
{
    double    r1,
              r2,
              r3,
              r4,
              r5,
              r6,
              r7,
              r8,
              r9;

    r1 = U[0] * tsai_cc.Ty;
    r2 = U[1] * tsai_cc.Ty;
    r3 = U[2] * tsai_cc.Ty;

    r4 = U[4] * tsai_cc.Ty;
    r5 = U[5] * tsai_cc.Ty;
    r6 = U[6] * tsai_cc.Ty;

    /* use the outer product of the first two rows to get the last row */
    r7 = r2 * r6 - r3 * r5;
    r8 = r3 * r4 - r1 * r6;
    r9 = r1 * r5 - r2 * r4;

    /* update the calibration constants */
    tsai_cc.r1 = r1;
    tsai_cc.r2 = r2;
    tsai_cc.r3 = r3;
    tsai_cc.r4 = r4;
    tsai_cc.r5 = r5;
    tsai_cc.r6 = r6;
    tsai_cc.r7 = r7;
    tsai_cc.r8 = r8;
    tsai_cc.r9 = r9;

    /* fill in tsai_cc.Rx, tsai_cc.Ry and tsai_cc.Rz */
    solve_RPY_transform ();
}


/************************************************************************/
/* pytsai: can fail; int return type is required. */
int epe_compute_Tx_Ty_Tz ()
{
    dmat      M,
              a,
              b;

    double    xk,
              yk,
              zk,
              Xu,
              Yu,
              Xd,
              Yd,
              distortion_factor;

    int       i,
              j;

    M = newdmat (0, (2 * tsai_cd.point_count - 1), 0, 2, &errno);
    if (errno) {
	pytsai_raise("epe compute Tx Ty Tz: unable to allocate matrix M");
	return 0;
    }

    a = newdmat (0, 2, 0, 0, &errno);
    if (errno) {
        freemat(M);
	pytsai_raise("epe compute Tx Ty Tz: unable to allocate vector a");
	return 0;
    }

    b = newdmat (0, (2 * tsai_cd.point_count - 1), 0, 0, &errno);
    if (errno) {
        freemat(M);
        freemat(a);
	pytsai_raise("epe compute Tx Ty Tz: unable to allocate vector b");
	return 0;
    }

    for (i = 0, j = tsai_cd.point_count; i < tsai_cd.point_count; i++, j++) {
	/* convert from world coordinates to untranslated camera coordinates */
	xk = tsai_cc.r1 * tsai_cd.xw[i] + tsai_cc.r2 * tsai_cd.yw[i] + tsai_cc.r3 * tsai_cd.zw[i];
	yk = tsai_cc.r4 * tsai_cd.xw[i] + tsai_cc.r5 * tsai_cd.yw[i] + tsai_cc.r6 * tsai_cd.zw[i];
	zk = tsai_cc.r7 * tsai_cd.xw[i] + tsai_cc.r8 * tsai_cd.yw[i] + tsai_cc.r9 * tsai_cd.zw[i];

	/* convert from image coordinates to distorted sensor coordinates */
	Xd = tsai_cp.dpx * (tsai_cd.Xf[i] - tsai_cp.Cx) / tsai_cp.sx;
	Yd = tsai_cp.dpy * (tsai_cd.Yf[i] - tsai_cp.Cy);

	/* convert from distorted sensor coordinates to undistorted sensor coordinates */
	distortion_factor = 1 + tsai_cc.kappa1 * (SQR (Xd) + SQR (Yd));
	Xu = Xd * distortion_factor;
	Yu = Yd * distortion_factor;

	M.el[i][0] = tsai_cc.f;
	M.el[i][1] = 0;
	M.el[i][2] = -Xu;
	b.el[i][0] = Xu * zk - tsai_cc.f * xk;

	M.el[j][0] = 0;
	M.el[j][1] = tsai_cc.f;
	M.el[j][2] = -Yu;
	b.el[j][0] = Yu * zk - tsai_cc.f * yk;
    }

    if (solve_system (M, a, b)) {
        freemat(M);
        freemat(a);
        freemat(b);
	pytsai_raise("epe compute Tx Ty Tz: unable to solve system  Ma=b");
	return 0;
    }

    tsai_cc.Tx = a.el[0][0];
    tsai_cc.Ty = a.el[1][0];
    tsai_cc.Tz = a.el[2][0];

    freemat (M);
    freemat (a);
    freemat (b);

    return 1;
}


/************************************************************************/
/* pytsai: cannot fail; void return type is fine. */
void epe_optimize_error (m_ptr, n_ptr, params, err)
    integer  *m_ptr;		/* pointer to number of points to fit */
    integer  *n_ptr;		/* pointer to number of parameters */
    doublereal *params;		/* vector of parameters */
    doublereal *err;		/* vector of error from data */
{
    int       i;

    double    xc,
              yc,
              zc,
              Xd,
              Yd,
              Xu_1,
              Yu_1,
              Xu_2,
              Yu_2,
              distortion_factor,
              Rx,
              Ry,
              Rz,
              Tx,
              Ty,
              Tz,
              r1,
              r2,
              r3,
              r4,
              r5,
              r6,
              r7,
              r8,
              r9,
              sa,
              sb,
              sg,
              ca,
              cb,
              cg;

    Rx = params[0];
    Ry = params[1];
    Rz = params[2];
    Tx = params[3];
    Ty = params[4];
    Tz = params[5];

    SINCOS (Rx, sa, ca);
    SINCOS (Ry, sb, cb);
    SINCOS (Rz, sg, cg);

    r1 = cb * cg;
    r2 = cg * sa * sb - ca * sg;
    r3 = sa * sg + ca * cg * sb;
    r4 = cb * sg;
    r5 = sa * sb * sg + ca * cg;
    r6 = ca * sb * sg - cg * sa;
    r7 = -sb;
    r8 = cb * sa;
    r9 = ca * cb;

    for (i = 0; i < tsai_cd.point_count; i++) {
	/* convert from world coordinates to camera coordinates */
	xc = r1 * tsai_cd.xw[i] + r2 * tsai_cd.yw[i] + r3 * tsai_cd.zw[i] + Tx;
	yc = r4 * tsai_cd.xw[i] + r5 * tsai_cd.yw[i] + r6 * tsai_cd.zw[i] + Ty;
	zc = r7 * tsai_cd.xw[i] + r8 * tsai_cd.yw[i] + r9 * tsai_cd.zw[i] + Tz;

	/* convert from camera coordinates to undistorted sensor plane coordinates */
	Xu_1 = tsai_cc.f * xc / zc;
	Yu_1 = tsai_cc.f * yc / zc;

	/* convert from 2D image coordinates to distorted sensor coordinates */
	Xd = tsai_cp.dpx * (tsai_cd.Xf[i] - tsai_cp.Cx) / tsai_cp.sx;
	Yd = tsai_cp.dpy * (tsai_cd.Yf[i] - tsai_cp.Cy);

	/* convert from distorted sensor coordinates to undistorted sensor plane coordinates */
	distortion_factor = 1 + tsai_cc.kappa1 * (SQR (Xd) + SQR (Yd));
	Xu_2 = Xd * distortion_factor;
	Yu_2 = Yd * distortion_factor;

	/* record the error in the undistorted sensor coordinates */
	err[i] = hypot (Xu_1 - Xu_2, Yu_1 - Yu_2);
    }
}


/* pytsai: can fail; int return type is required. */
int epe_optimize ()
{
#define NPARAMS 6

    int       i;

    /* Parameters needed by MINPACK's lmdif() */

    integer     m = tsai_cd.point_count;
    integer     n = NPARAMS;
    doublereal  x[NPARAMS];
    doublereal *fvec;
    doublereal  ftol = REL_SENSOR_TOLERANCE_ftol;
    doublereal  xtol = REL_PARAM_TOLERANCE_xtol;
    doublereal  gtol = ORTHO_TOLERANCE_gtol;
    integer     maxfev = MAXFEV;
    doublereal  epsfcn = EPSFCN;
    doublereal  diag[NPARAMS];
    integer     mode = MODE;
    doublereal  factor = FACTOR;
    integer     nprint = 0;
    integer     info;
    integer     nfev;
    doublereal *fjac;
    integer     ldfjac = m;
    integer     ipvt[NPARAMS];
    doublereal  qtf[NPARAMS];
    doublereal  wa1[NPARAMS];
    doublereal  wa2[NPARAMS];
    doublereal  wa3[NPARAMS];
    doublereal *wa4;

    /* allocate some workspace */
    if (( fvec = (doublereal *) PyMem_Malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       pytsai_raise("PyMem_Malloc: Cannot allocate workspace fvec");
       return 0;
    }
 
    if (( fjac = (doublereal *) PyMem_Malloc ((unsigned int) m*n * (unsigned int) sizeof(doublereal))) == NULL ) {
       PyMem_Free(fvec);
       pytsai_raise("PyMem_Malloc: Cannot allocate workspace fjac");
       return 0;
    }
 
    if (( wa4 = (doublereal *) PyMem_Malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       PyMem_Free(fvec);
       PyMem_Free(fjac);
       pytsai_raise("PyMem_Malloc: Cannot allocate workspace wa4");
       return 0;
    }

    /* use the current calibration and camera constants as a starting point */
    x[0] = tsai_cc.Rx;
    x[1] = tsai_cc.Ry;
    x[2] = tsai_cc.Rz;
    x[3] = tsai_cc.Tx;
    x[4] = tsai_cc.Ty;
    x[5] = tsai_cc.Tz;

    /* define optional scale factors for the parameters */
    if ( mode == 2 ) {
        for (i = 0; i < NPARAMS; i++)
            diag[i] = 1.0;             /* some user-defined values */
    }
       
    /* perform the optimization */
    lmdif_ (epe_optimize_error,
            &m, &n, x, fvec, &ftol, &xtol, &gtol, &maxfev, &epsfcn,
            diag, &mode, &factor, &nprint, &info, &nfev, fjac, &ldfjac,
            ipvt, qtf, wa1, wa2, wa3, wa4);
    /* TODO: Check for error conditions (into parameter). */

    /* update the calibration and camera constants */
    tsai_cc.Rx = x[0];
    tsai_cc.Ry = x[1];
    tsai_cc.Rz = x[2];
    apply_RPY_transform ();

    tsai_cc.Tx = x[3];
    tsai_cc.Ty = x[4];
    tsai_cc.Tz = x[5];

    /* release allocated workspace */
    PyMem_Free(fvec);
    PyMem_Free(fjac);
    PyMem_Free(wa4);

#ifdef DEBUG
    /* print the number of function calls during iteration */
    fprintf(stderr,"info: %d nfev: %d\n\n",info,nfev);
#endif

    return 1;

#undef NPARAMS
}


/************************************************************************/
/* pytsai: can fail; int return type is required. */
int coplanar_extrinsic_parameter_estimation ()
{
    double    trial_f,
              U[5];

    if (!cepe_compute_U(U))
        return 0;

    cepe_compute_Tx_and_Ty (U);

    cepe_compute_R (U);

    if (!cepe_compute_approximate_f(&trial_f))
        return 0;

    if (trial_f < 0) {
	/* try the other rotation matrix solution */
	tsai_cc.r3 = -tsai_cc.r3;
	tsai_cc.r6 = -tsai_cc.r6;
	tsai_cc.r7 = -tsai_cc.r7;
	tsai_cc.r8 = -tsai_cc.r8;
	solve_RPY_transform ();

	if (!cepe_compute_approximate_f(&trial_f))
                return 0;

	if (trial_f < 0) {
	    pytsai_raise("error - possible handedness problem with data");
	    return 0;
	}
    }

    if (!epe_compute_Tx_Ty_Tz())
        return 0;

    if (!epe_optimize())
        return 0;

    return 1;
}


/************************************************************************/
/* pytsai: can fail; int return type is required. */
int noncoplanar_extrinsic_parameter_estimation ()
{
    double    U[7];

    if (!ncepe_compute_U(U))
        return 0;

    ncepe_compute_Tx_and_Ty (U);

    ncepe_compute_R (U);

    if (!epe_compute_Tx_Ty_Tz())
        return 0;

    if (!epe_optimize())
        return 0;

    return 1;
}
