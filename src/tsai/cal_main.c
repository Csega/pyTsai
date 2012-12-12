/**
 * cal_main.c
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
* This file contains routines for calibrating Tsai's 11 parameter camera     *
* model. The camera model is based on the pin hole model of 3D-2D            *
* perspective projection with 1st order radial lens distortion.  The model   *
* consists of 5 internal (also called intrinsic or interior) camera          *
* parameters:                                                                *
*                                                                            *
*       f       - effective focal length of the pin hole camera              *
*       kappa1  - 1st order radial lens distortion coefficient               *
*       Cx, Cy  - coordinates of center of radial lens distortion            *
*		  (also used as the piercing point of the camera coordinate  *
*		   frame's Z axis with the camera's sensor plane)            *
*       sx      - uncertainty factor for scale of horizontal scanline        *
*                                                                            *
* and 6 external (also called extrinsic or exterior) camera parameters:      *
*                                                                            *
*       Rx, Ry, Rz, Tx, Ty, Tz  - rotational and translational components of *
*                                 the transform between the world's          *
*                                 coordinate frame and the camera's          *
*                                 coordinate frame.                          *
*                                                                            *
* Data for model calibration consists of the (x,y,z) world coordinates of a  *
* feature point (in mm) and the corresponding coordinates (Xf,Yf) (in        *
* pixels) of the feature point in the image.  Two types of calibration are   *
* available:                                                                 *
*                                                                            *
*       coplanar     - all of the calibration points lie in a single plane   *
*       non-coplanar - the calibration points do not lie in a single plane   *
*                                                                            *
* This file contains routines for two levels of calibration.  The first      *
* level of calibration is a direct implementation of Tsai's algorithm in     *
* which only the f, Tz and kappa1 parameters are optimized for.  The         *
* routines are:                                                              *
*                                                                            *
*       coplanar_calibration ()                                              *
*       noncoplanar_calibration ()                                           *
*                                                                            *
* The second level of calibration optimizes for everything.  This level is   *
* very slow but provides the most accurate calibration.  The routines are:   *
*                                                                            *
*       coplanar_calibration_with_full_optimization ()                       *
*       noncoplanar_calibration_with_full_optimization ()                    *
*                                                                            *
* Routines are also provided for initializing camera parameter variables     *
* for five of our camera/frame grabber systems.  These routines are:         *
*                                                                            *
*       initialize_photometrics_parms ()                                     *
*       initialize_general_imaging_mos5300_matrox_parms ()                   *
*       initialize_panasonic_gp_mf702_matrox_parms ()                        *
*       initialize_sony_xc75_matrox_parms ()                                 *
*       initialize_sony_xc77_matrox_parms ()                                 *
*       initialize_sony_xc57_androx_parms ()                                 *
*                                                                            *
*                                                                            *
* External routines                                                          *
* -----------------                                                          *
*                                                                            *
* Nonlinear optimization for these camera calibration routines is performed  *
* by the MINPACK lmdif subroutine.  lmdif uses a modified                    *
* Levenberg-Marquardt with a jacobian calculated by a forward-difference     *
* approximation. The MINPACK FORTRAN routines were translated into C         *
* generated using f2c.                                                       *
*                                                                            *
* Matrix operations (inversions, multiplications, etc.) are also provided by *
* external routines.                                                         *
*                                                                            *
*                                                                            *
* Extra notes                                                                *
* -----------                                                                *
*                                                                            *
* An explanation of the basic algorithms and description of the variables    *
* can be found in several publications, including:                           *
*                                                                            *
* "An Efficient and Accurate Camera Calibration Technique for 3D Machine     *
*  Vision", Roger Y. Tsai, Proceedings of IEEE Conference on Computer Vision *
*  and Pattern Recognition, Miami Beach, FL, 1986, pages 364-374.            *
*                                                                            *
*  and                                                                       *
*                                                                            *
* "A versatile Camera Calibration Technique for High-Accuracy 3D Machine     *
*  Vision Metrology Using Off-the-Shelf TV Cameras and Lenses", Roger Y.     *
*  Tsai, IEEE Journal of Robotics and Automation, Vol. RA-3, No. 4, August   *
*  1987, pages 323-344.                                                      *
*                                                                            *
*                                                                            *
* Notation                                                                   *
* --------                                                                   *
*                                                                            *
* The camera's X axis runs along increasing column coordinates in the        *
* image/frame.  The Y axis runs along increasing row coordinates.            *
*                                                                            *
* pix == image/frame grabber picture element                                 *
* sel == camera sensor element                                               *
*                                                                            *
* Internal routines starting with "cc_" are for coplanar calibration.        *
* Internal routines starting with "ncc_" are for noncoplanar calibration.    *
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
* 20-May-95  Reg Willson (rgwillson@mmm.com) at 3M St. Paul, MN              *
*       Return the error to lmdif rather than the squared error.             *
*         lmdif calculates the squared error internally during optimization. *
*         Before this change calibration was essentially optimizing error^4. *
*       Put transform and evaluation routines into separate files.           *
*                                                                            *
* 02-Apr-95  Reg Willson (rgwillson@mmm.com) at 3M St. Paul, MN              *
*       Rewrite memory allocation to avoid memory alignment problems         *
*       on some machines.                                                    *
*       Strip out IMSL code.  MINPACK seems to work fine.                    *
*       Filename changes for DOS port.                                       *
*                                                                            *
* 04-Jun-94  Reg Willson (rgwillson@mmm.com) at 3M St. Paul, MN              *
*       Replaced ncc_compute_Xdp_and_Ydp with                                *
*       ncc_compute_Xd_Yd_and_r_squared.                                     *
*         (effectively propagates the 22-Mar-94 to the non-coplanar          *
*         routines)                                                          *
*       Added alternate macro definitions for less common math functions.    *
*                                                                            *
* 25-Mar-94  Torfi Thorhallsson (torfit@verk.hi.is) at the University of     *
*            Iceland                                                         *
*       Added a new version of the routines:                                 *
*            cc_compute_exact_f_and_Tz ()                                    *
*            cc_five_parm_optimization_with_late_distortion_removal ()       *
*            cc_five_parm_optimization_with_early_distortion_removal ()      *
*            cc_nic_optimization ()                                          *
*            cc_full_optimization ()                                         *
*            ncc_compute_exact_f_and_Tz ()                                   *
*            ncc_nic_optimization ()                                         *
*            ncc_full_optimization ()                                        *
*                                                                            *
*       The new routines use the *public domain* MINPACK library for         *
*       optimization instead of the commercial IMSL library.                 *
*       To select the new routines, compile this file with the flag          *
*       -DMINPACK                                                            *
*                                                                            *
* 22-Mar-94  Torfi Thorhallsson (torfit@verk.hi.is) at the University of     *
*            Iceland                                                         *
*       Fixed a bug in cc_nic_optimization_error and                         *
*       cc_full_optimization_error.                                          *
*       A division by cp.sx was missing.                                     *
*                                                                            *
* 15-Feb-94  Reg Willson (rgw@cs.cmu.edu) at Carnegie-Mellon University      *
*       Included Frederic Devernay's (<Frederic.Devernay@sophia.inria.fr)    *
*	significantly improved routine for converting from undistorted to    *
*	distorted sensor coordinates.  Rather than iteratively solving a     *
*	system of two non-linear equations to perform the conversion, the    *
*	new routine algebraically solves a cubic polynomial in Rd (using     *
*	the Cardan method).                                                  *
*                                                                            *
* 14-Feb-94  Reg Willson (rgw@cs.cmu.edu) at Carnegie-Mellon University      *
*	Fixed a coding bug in ncc_compute_R and ncc_compute_better_R.        *
*       The r4, r5, and r6 terms should not be divided by cp.sx.             *
*       Bug reported by: Volker Rodehorst <vr@cs.tu-berlin.de>               *
*                                                                            *
* 04-Jul-93  Reg Willson (rgw@cs.cmu.edu) at Carnegie-Mellon University      *
*	Added new routines to evaluate the accuracy of camera calibration.   *
*                                                                            *
*       Added check for coordinate handedness problem in calibration data.   *
*                                                                            *
* 01-May-93  Reg Willson (rgw@cs.cmu.edu) at Carnegie-Mellon University      *
*	For efficiency the non-linear optimizations are now all based on     *
*	the minimization of squared error in undistorted image coordinates   *
*	instead of the squared error in distorted image coordinates.         *
*                                                                            *
*	New routine for inverse perspective projection.                      *
*                                                                            *
* 14-Feb-93  Reg Willson (rgw@cs.cmu.edu) at Carnegie-Mellon University      *
*       Bug fixes and speed ups.                                             *
*                                                                            *
* 07-Feb-93  Reg Willson (rgw@cs.cmu.edu) at Carnegie-Mellon University      *
*       Original implementation.                                             *
*                                                                            *
\****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "../matrix/matrix.h"
#include "cal_main.h"
#include "../minpack/f2c.h"
#include "../minpack/minpack.h"
#include "../errors.h"


/* Variables used by the subroutines for I/O (perhaps not the best way of 
 * doing this) */
struct camera_parameters tsai_cp;
struct calibration_data tsai_cd;
struct calibration_constants tsai_cc;

/* Local working storage */
static
double Xd[MAX_POINTS],
       Yd[MAX_POINTS],
       r_squared[MAX_POINTS],
       U[7];


/* char camera_type[256] = "unknown"; */


/** We don't need the initialization routines in the Python version; these
  * are provided by the high-level Python wrapper. */
#if 0

/************************************************************************/
void      initialize_photometrics_parms ()
{
    strcpy (camera_type, "Photometrics Star I");

    tsai_cp.Ncx = 576;		/* [sel]        */
    tsai_cp.Nfx = 576;		/* [pix]        */
    tsai_cp.dx = 0.023;		/* [mm/sel]     */
    tsai_cp.dy = 0.023;		/* [mm/sel]     */
    tsai_cp.dpx = tsai_cp.dx * tsai_cp.Ncx / tsai_cp.Nfx;	/* [mm/pix]     */
    tsai_cp.dpy = tsai_cp.dy;		/* [mm/pix]     */
    tsai_cp.Cx = 576 / 2;		/* [pix]        */
    tsai_cp.Cy = 384 / 2;		/* [pix]        */
    tsai_cp.sx = 1.0;		/* []		 */

    /* something a bit more realistic */
    tsai_cp.Cx = 258.0;
    tsai_cp.Cy = 204.0;
}


/************************************************************************/
void      initialize_general_imaging_mos5300_matrox_parms ()
{
    strcpy (camera_type, "General Imaging MOS5300 + Matrox");

    tsai_cp.Ncx = 649;		/* [sel]        */
    tsai_cp.Nfx = 512;		/* [pix]        */
    tsai_cp.dx = 0.015;		/* [mm/sel]     */
    tsai_cp.dy = 0.015;		/* [mm/sel]     */
    tsai_cp.dpx = tsai_cp.dx * tsai_cp.Ncx / tsai_cp.Nfx;	/* [mm/pix]     */
    tsai_cp.dpy = tsai_cp.dy;		/* [mm/pix]     */
    tsai_cp.Cx = 512 / 2;		/* [pix]        */
    tsai_cp.Cy = 480 / 2;		/* [pix]        */
    tsai_cp.sx = 1.0;		/* []		 */
}


/************************************************************************/
void      initialize_panasonic_gp_mf702_matrox_parms ()
{
    strcpy (camera_type, "Panasonic GP-MF702 + Matrox");

    tsai_cp.Ncx = 649;		/* [sel]        */
    tsai_cp.Nfx = 512;		/* [pix]        */
    tsai_cp.dx = 0.015;		/* [mm/sel]     */
    tsai_cp.dy = 0.015;		/* [mm/sel]     */
    tsai_cp.dpx = tsai_cp.dx * tsai_cp.Ncx / tsai_cp.Nfx;	/* [mm/pix]     */
    tsai_cp.dpy = tsai_cp.dy;		/* [mm/pix]     */

    tsai_cp.Cx = 268;		/* [pix]        */
    tsai_cp.Cy = 248;		/* [pix]        */

    tsai_cp.sx = 1.078647;		/* []           */
}


/************************************************************************/
void      initialize_sony_xc75_matrox_parms ()
{
    strcpy (camera_type, "Sony XC75 + Matrox");

    tsai_cp.Ncx = 768;		/* [sel]        */
    tsai_cp.Nfx = 512;		/* [pix]        */
    tsai_cp.dx = 0.0084;		/* [mm/sel]     */
    tsai_cp.dy = 0.0098;		/* [mm/sel]     */
    tsai_cp.dpx = tsai_cp.dx * tsai_cp.Ncx / tsai_cp.Nfx;	/* [mm/pix]     */
    tsai_cp.dpy = tsai_cp.dy;		/* [mm/pix]     */
    tsai_cp.Cx = 512 / 2;		/* [pix]        */
    tsai_cp.Cy = 480 / 2;		/* [pix]        */
    tsai_cp.sx = 1.0;		/* []           */
}


/************************************************************************/
void      initialize_sony_xc77_matrox_parms ()
{
    strcpy (camera_type, "Sony XC77 + Matrox");

    tsai_cp.Ncx = 768;		/* [sel]        */
    tsai_cp.Nfx = 512;		/* [pix]        */
    tsai_cp.dx = 0.011;		/* [mm/sel]     */
    tsai_cp.dy = 0.013;		/* [mm/sel]     */
    tsai_cp.dpx = tsai_cp.dx * tsai_cp.Ncx / tsai_cp.Nfx;	/* [mm/pix]     */
    tsai_cp.dpy = tsai_cp.dy;		/* [mm/pix]     */
    tsai_cp.Cx = 512 / 2;		/* [pix]        */
    tsai_cp.Cy = 480 / 2;		/* [pix]        */
    tsai_cp.sx = 1.0;		/* []           */
}


/************************************************************************/
void      initialize_sony_xc57_androx_parms ()
{
    strcpy (camera_type, "Sony XC57 + Androx");

    tsai_cp.Ncx = 510;               /* [sel]        */
    tsai_cp.Nfx = 512;               /* [pix]        */
    tsai_cp.dx = 0.017;              /* [mm/sel]     */
    tsai_cp.dy = 0.013;              /* [mm/sel]     */
    tsai_cp.dpx = tsai_cp.dx * tsai_cp.Ncx / tsai_cp.Nfx;   /* [mm/pix]     */
    tsai_cp.dpy = tsai_cp.dy;             /* [mm/pix]     */
    tsai_cp.Cx = 512 / 2;            /* [pix]        */
    tsai_cp.Cy = 480 / 2;            /* [pix]        */
    tsai_cp.sx = 1.107914;           /* []           */
}


/************************************************************************/
/* John Krumm, 9 November 1993                                          */
/* This assumes that every other row, starting with the second row from */
/* the top, has been removed.  The Xap Shot CCD only has 250 lines, and */
/* it inserts new rows in between the real rows to make a full size     */
/* picture.  Its algorithm for making these new rows isn't very good,   */
/* so it's best just to throw them away.                                */
/* The camera's specs give the CCD size as 6.4(V)x4.8(H) mm.            */
/* A call to the manufacturer found that the CCD has 250 rows           */
/* and 782 columns.  It is underscanned to 242 rows and 739 columns.    */
/************************************************************************/
void      initialize_xapshot_matrox_parms ()
{
    strcpy (camera_type, "Canon Xap Shot");

    tsai_cp.Ncx = 739;			/* [sel]        */
    tsai_cp.Nfx = 512;			/* [pix]        */
    tsai_cp.dx = 6.4 / 782.0;		/* [mm/sel]     */
    tsai_cp.dy = 4.8 / 250.0;		/* [mm/sel]     */
    tsai_cp.dpx = tsai_cp.dx * tsai_cp.Ncx / tsai_cp.Nfx;	/* [mm/pix]     */
    tsai_cp.dpy = tsai_cp.dy;			/* [mm/pix]     */
    tsai_cp.Cx = 512 / 2;			/* [pix]        */
    tsai_cp.Cy = 240 / 2;			/* [pix]        */
    tsai_cp.sx = 1.027753;	/* from noncoplanar calibration *//* []           */
}

#endif  /* finish excluding camera settings routines */


/***********************************************************************\
* This routine solves for the roll, pitch and yaw angles (in radians)	*
* for a given orthonormal rotation matrix (from Richard P. Paul,        *
* Robot Manipulators: Mathematics, Programming and Control, p70).       *
* Note 1, should the rotation matrix not be orthonormal these will not  *
* be the "best fit" roll, pitch and yaw angles.                         *
* Note 2, there are actually two possible solutions for the matrix.     *
* The second solution can be found by adding 180 degrees to Rz before   *
* Ry and Rx are calculated.                                             *
\***********************************************************************/
/* pytsai: cannot fail; void return type is fine. */
void solve_RPY_transform (TSAI_ARGS_DECL)
{
    double    sg,
              cg;

    tsai_cc.Rz = atan2 (tsai_cc.r4, tsai_cc.r1);
    SINCOS (tsai_cc.Rz, sg, cg);
    tsai_cc.Ry = atan2 (-tsai_cc.r7, tsai_cc.r1 * cg + tsai_cc.r4 * sg);
    tsai_cc.Rx = atan2 (tsai_cc.r3 * sg - tsai_cc.r6 * cg, tsai_cc.r5 * cg - tsai_cc.r2 * sg);
}


/***********************************************************************\
* This routine simply takes the roll, pitch and yaw angles and fills in	*
* the rotation matrix elements r1-r9.					*
\***********************************************************************/
/* pytsai: cannot fail; void return type is fine. */
void apply_RPY_transform (TSAI_ARGS_DECL)
{
    double    sa,
              ca,
              sb,
              cb,
              sg,
              cg;

    SINCOS (tsai_cc.Rx, sa, ca);
    SINCOS (tsai_cc.Ry, sb, cb);
    SINCOS (tsai_cc.Rz, sg, cg);

    tsai_cc.r1 = cb * cg;
    tsai_cc.r2 = cg * sa * sb - ca * sg;
    tsai_cc.r3 = sa * sg + ca * cg * sb;
    tsai_cc.r4 = cb * sg;
    tsai_cc.r5 = sa * sb * sg + ca * cg;
    tsai_cc.r6 = ca * sb * sg - cg * sa;
    tsai_cc.r7 = -sb;
    tsai_cc.r8 = cb * sa;
    tsai_cc.r9 = ca * cb;
}


/***********************************************************************\
* Routines for coplanar camera calibration	 			*
\***********************************************************************/
/* pytsai: cannot fail; void return type is fine. */
void cc_compute_Xd_Yd_and_r_squared ()
{
    int       i;

    double    Xd_,
              Yd_;

    for (i = 0; i < tsai_cd.point_count; i++) {
	Xd[i] = Xd_ = tsai_cp.dpx * (tsai_cd.Xf[i] - tsai_cp.Cx) / tsai_cp.sx;	/* [mm] */
	Yd[i] = Yd_ = tsai_cp.dpy * (tsai_cd.Yf[i] - tsai_cp.Cy);	        /* [mm] */
	r_squared[i] = SQR (Xd_) + SQR (Yd_);                   /* [mm^2] */
    }
}


/* pytsai: can fail: need int return type. */
int cc_compute_U ()
{
    int       i;

    dmat      M,
              a,
              b;

    M = newdmat (0, (tsai_cd.point_count - 1), 0, 4, &errno);
    if (errno) {
	pytsai_raise ("cc compute U: unable to allocate matrix M");
        return 0;
    }

    a = newdmat (0, 4, 0, 0, &errno);
    if (errno) {
        freemat(M);
	pytsai_raise ("cc compute U: unable to allocate vector a");
	return 0;
    }

    b = newdmat (0, (tsai_cd.point_count - 1), 0, 0, &errno);
    if (errno) {
        freemat(M);
        freemat(a);
	pytsai_raise ("cc compute U: unable to allocate vector b");
	return 0;
    }

    for (i = 0; i < tsai_cd.point_count; i++) {
	M.el[i][0] = Yd[i] * tsai_cd.xw[i];
	M.el[i][1] = Yd[i] * tsai_cd.yw[i];
	M.el[i][2] = Yd[i];
	M.el[i][3] = -Xd[i] * tsai_cd.xw[i];
	M.el[i][4] = -Xd[i] * tsai_cd.yw[i];
	b.el[i][0] = Xd[i];
    }

    if (solve_system (M, a, b)) {
        freemat(M);
        freemat(a);
        freemat(b);
	pytsai_raise ("cc compute U: unable to solve system  Ma=b");
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
void cc_compute_Tx_and_Ty ()
{
    int       i,
              far_point;

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
	if ((distance = r_squared[i]) > far_distance) {
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
    if ((SIGNBIT (x) != SIGNBIT (Xd[far_point])) ||
	(SIGNBIT (y) != SIGNBIT (Yd[far_point])))
	Ty = -Ty;

    /* update the calibration constants */
    tsai_cc.Tx = U[2] * Ty;
    tsai_cc.Ty = Ty;
}


/* pytsai: cannot fail; void return type is fine. */
void cc_compute_R ()
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

 
/* pytsai: can fail; need int return type */
int cc_compute_approximate_f_and_Tz ()
{
    int       i;

    dmat      M,
              a,
              b;

    M = newdmat (0, (tsai_cd.point_count - 1), 0, 1, &errno);
    if (errno) {
	pytsai_raise ("cc compute apx: unable to allocate matrix M");
	return 0;
    }

    a = newdmat (0, 1, 0, 0, &errno);
    if (errno) {
        freemat(M);
	pytsai_raise ("cc compute apx: unable to allocate vector a");
	return 0;
    }

    b = newdmat (0, (tsai_cd.point_count - 1), 0, 0, &errno);
    if (errno) {
        freemat(M);
        freemat(a);
	pytsai_raise ("cc compute apx: unable to allocate vector b");
	return 0;
    }

    for (i = 0; i < tsai_cd.point_count; i++) {
	M.el[i][0] = tsai_cc.r4 * tsai_cd.xw[i] + tsai_cc.r5 * tsai_cd.yw[i] + tsai_cc.Ty;
	M.el[i][1] = -Yd[i];
	b.el[i][0] = (tsai_cc.r7 * tsai_cd.xw[i] + tsai_cc.r8 * tsai_cd.yw[i]) * Yd[i];
    }

    if (solve_system (M, a, b)) {
        freemat(M);
        freemat(a);
        freemat(b);
	pytsai_raise ("cc compute apx: unable to solve system  Ma=b");
	return 0;
    }

    /* update the calibration constants */
    tsai_cc.f = a.el[0][0];
    tsai_cc.Tz = a.el[1][0];
    tsai_cc.kappa1 = 0.0;  /* this is the assumption that our calculation was 
                       * made under */

    freemat (M);
    freemat (a);
    freemat (b);

    return 1;
}


/************************************************************************/
/* pytsai: cannot fail; void return type is fine. */
void cc_compute_exact_f_and_Tz_error (
        integer *m_ptr,         /* pointer to number of points to fit */
        integer *n_ptr,         /* pointer to number of parameters */
        doublereal *params,     /* vector of parameters */
        doublereal *err         /* vector of error from data */
)
{
    int       i;

    double    f,
              Tz,
              kappa1,
              xc,
              yc,
              zc,
              Xu_1,
              Yu_1,
              Xu_2,
              Yu_2,
              distortion_factor;

    f = params[0];
    Tz = params[1];
    kappa1 = params[2];

    for (i = 0; i < tsai_cd.point_count; i++) {
	/* convert from world coordinates to camera coordinates */
	/* Note: zw is implicitly assumed to be zero for these (coplanar) 
	 * calculations */
	xc = tsai_cc.r1 * tsai_cd.xw[i] + tsai_cc.r2 * tsai_cd.yw[i] + tsai_cc.Tx;
	yc = tsai_cc.r4 * tsai_cd.xw[i] + tsai_cc.r5 * tsai_cd.yw[i] + tsai_cc.Ty;
	zc = tsai_cc.r7 * tsai_cd.xw[i] + tsai_cc.r8 * tsai_cd.yw[i] + Tz;

	/* convert from camera coordinates to undistorted sensor coordinates */
	Xu_1 = f * xc / zc;
	Yu_1 = f * yc / zc;

	/* convert from distorted sensor coordinates to undistorted sensor 
	 * coordinates */
	distortion_factor = 1 + kappa1 * (SQR (Xd[i]) + SQR (Yd[i]));
	Xu_2 = Xd[i] * distortion_factor;
	Yu_2 = Yd[i] * distortion_factor;

        /* record the error in the undistorted sensor coordinates */
        err[i] = hypot (Xu_1 - Xu_2, Yu_1 - Yu_2);
    }
}


/* pytsai: can fail; need int return type */
int cc_compute_exact_f_and_Tz ()
{
#define NPARAMS 3

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
    if (( fvec = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       pytsai_raise("malloc: Cannot allocate workspace fvec");
       return 0;
    }

    if (( fjac = (doublereal *) malloc ((unsigned int) m*n * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       pytsai_raise("malloc: Cannot allocate workspace fjac");
       return 0;
    }

    if (( wa4 = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       free(fjac);
       pytsai_raise("malloc: Cannot allocate workspace wa4");
       return 0;
    }

    /* use the current calibration constants as an initial guess */
    x[0] = tsai_cc.f;
    x[1] = tsai_cc.Tz;
    x[2] = tsai_cc.kappa1;

    /* define optional scale factors for the parameters */
    if ( mode == 2 ) {
        for (i = 0; i < NPARAMS; i++)
	    diag[i] = 1.0;             /* some user-defined values */
    }

    /* perform the optimization */ 
    lmdif_ (cc_compute_exact_f_and_Tz_error,
            &m, &n, x, fvec, &ftol, &xtol, &gtol, &maxfev, &epsfcn,
            diag, &mode, &factor, &nprint, &info, &nfev, fjac, &ldfjac,
            ipvt, qtf, wa1, wa2, wa3, wa4);

    /* update the calibration constants */
    tsai_cc.f = x[0];
    tsai_cc.Tz = x[1];
    tsai_cc.kappa1 = x[2];

    /* release allocated workspace */
    free(fvec);
    free(fjac);
    free(wa4);

#ifdef DEBUG
    /* print the number of function calls during iteration */
    fprintf(stderr,"info: %d nfev: %d\n\n",info,nfev);
#endif

    return 1;

#undef NPARAMS
}


/************************************************************************/
/* pytsai: can fail; int return type is required. */
int cc_three_parm_optimization ()
{
    int       i;

    for (i = 0; i < tsai_cd.point_count; i++)
	if (tsai_cd.zw[i]) {
	    pytsai_raise("error - coplanar calibration tried with data outside of Z plane");
	    return 0;
	}

    cc_compute_Xd_Yd_and_r_squared ();

    if (!cc_compute_U())
        return 0;

    cc_compute_Tx_and_Ty ();

    cc_compute_R ();

    if (!cc_compute_approximate_f_and_Tz())
        return 0;

    if (tsai_cc.f < 0) {
        /* try the other solution for the orthonormal matrix */
	tsai_cc.r3 = -tsai_cc.r3;
	tsai_cc.r6 = -tsai_cc.r6;
	tsai_cc.r7 = -tsai_cc.r7;
	tsai_cc.r8 = -tsai_cc.r8;
	solve_RPY_transform ();

	if (!cc_compute_approximate_f_and_Tz())
                return 0;

	if (tsai_cc.f < 0) {
	    pytsai_raise("error - possible handedness problem with data");
	    return 0;
	}
    }

    if (!cc_compute_exact_f_and_Tz())
        return 0;

    return 1;
}


/************************************************************************/
/* pytsai: cannot fail; void return type is fine. */
void cc_remove_sensor_plane_distortion_from_Xd_and_Yd ()
{
    int       i;
    double    Xu,
              Yu;

    for (i = 0; i < tsai_cd.point_count; i++) {
	distorted_to_undistorted_sensor_coord (Xd[i], Yd[i], &Xu, &Yu);
	Xd[i] = Xu;
	Yd[i] = Yu;
	r_squared[i] = SQR (Xu) + SQR (Yu);
    }
}


/************************************************************************/
/* pytsai: can fail.  This method is called from within the lmdif_ routine.
 * To indicate failure, set *iflag = -1. */
void cc_five_parm_optimization_with_late_distortion_removal_error (m_ptr, n_ptr, params, err, iflag)
    integer  *m_ptr;		/* pointer to number of points to fit */
    integer  *n_ptr;		/* pointer to number of parameters */
    doublereal *params;		/* vector of parameters */
    doublereal *err;		/* vector of error from data */
    integer  *iflag;            /* flag to indicate error to caller */
{
    int       i;

    double    f,
              Tz,
              kappa1,
              xc,
              yc,
              zc,
              Xu_1,
              Yu_1,
              Xu_2,
              Yu_2,
              distortion_factor;

    /* in this routine radial lens distortion is only taken into account */
    /* after the rotation and translation constants have been determined */

    f = params[0];
    Tz = params[1];
    kappa1 = params[2];

    tsai_cp.Cx = params[3];
    tsai_cp.Cy = params[4];

    cc_compute_Xd_Yd_and_r_squared ();

    if (!cc_compute_U())
    {
        *iflag = -1;
        return;
    }

    cc_compute_Tx_and_Ty ();

    cc_compute_R ();

    if (!cc_compute_approximate_f_and_Tz())
    {
        *iflag = -1;
        return;
    }

    if (tsai_cc.f < 0) {
        /* try the other solution for the orthonormal matrix */
	tsai_cc.r3 = -tsai_cc.r3;
	tsai_cc.r6 = -tsai_cc.r6;
	tsai_cc.r7 = -tsai_cc.r7;
	tsai_cc.r8 = -tsai_cc.r8;
	solve_RPY_transform ();

        if (!cc_compute_approximate_f_and_Tz())
        {
                *iflag = -1;
                return;
        }

        if (tsai_cc.f < 0) {
            pytsai_raise("error - possible handedness problem with data");
            *iflag = -1;
            return;
	}
    }

    for (i = 0; i < tsai_cd.point_count; i++) {
	/* convert from world coordinates to camera coordinates */
	/* Note: zw is implicitly assumed to be zero for these (coplanar) calculations */
	xc = tsai_cc.r1 * tsai_cd.xw[i] + tsai_cc.r2 * tsai_cd.yw[i] + tsai_cc.Tx;
	yc = tsai_cc.r4 * tsai_cd.xw[i] + tsai_cc.r5 * tsai_cd.yw[i] + tsai_cc.Ty;
	zc = tsai_cc.r7 * tsai_cd.xw[i] + tsai_cc.r8 * tsai_cd.yw[i] + Tz;

	/* convert from camera coordinates to undistorted sensor coordinates */
	Xu_1 = f * xc / zc;
	Yu_1 = f * yc / zc;

	/* convert from distorted sensor coordinates to undistorted sensor coordinates */
	distortion_factor = 1 + kappa1 * r_squared[i];
	Xu_2 = Xd[i] * distortion_factor;
	Yu_2 = Yd[i] * distortion_factor;

        /* record the error in the undistorted sensor coordinates */
        err[i] = hypot (Xu_1 - Xu_2, Yu_1 - Yu_2);
    }

}


/* pytsai: can fail; need int return type. */
int cc_five_parm_optimization_with_late_distortion_removal ()
{
#define NPARAMS 5

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
    if (( fvec = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       pytsai_raise("malloc: Cannot allocate workspace fvec");
       return 0;
    }
 
    if (( fjac = (doublereal *) malloc ((unsigned int) m*n * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       pytsai_raise("malloc: Cannot allocate workspace fjac");
       return 0;
    }
 
    if (( wa4 = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       free(fjac);
       pytsai_raise("malloc: Cannot allocate workspace wa4");
       return 0;
    }
 
    /* use the current calibration constants as an initial guess */
    x[0] = tsai_cc.f;
    x[1] = tsai_cc.Tz;
    x[2] = tsai_cc.kappa1;
    x[3] = tsai_cp.Cx;
    x[4] = tsai_cp.Cy;

    /* define optional scale factors for the parameters */
    if ( mode == 2 ) {
        for (i = 0; i < NPARAMS; i++)
            diag[i] = 1.0;             /* some user-defined values */
    }
 
    /* perform the optimization */
    lmdif_ (cc_five_parm_optimization_with_late_distortion_removal_error,
            &m, &n, x, fvec, &ftol, &xtol, &gtol, &maxfev, &epsfcn,
            diag, &mode, &factor, &nprint, &info, &nfev, fjac, &ldfjac,
            ipvt, qtf, wa1, wa2, wa3, wa4);
    /* check for pytsai error condition */
    if (info == -1)
    {
        free(fvec);
        free(fjac);
        free(wa4);
        return 0;
    }
    /* TODO: Check for and translate other error conditions.  See
     * lmdif.c for possible values of the info parameter. */
 
    /* update the calibration and camera constants */
    tsai_cc.f = x[0];
    tsai_cc.Tz = x[1];
    tsai_cc.kappa1 = x[2];
    tsai_cp.Cx = x[3];
    tsai_cp.Cy = x[4];

    /* release allocated workspace */
    free(fvec);
    free(fjac);
    free(wa4);

#ifdef DEBUG
    /* print the number of function calls during iteration */
    fprintf(stderr,"info: %d nfev: %d\n\n",info,nfev);
#endif
 
    return 1;

#undef NPARAMS
}


/************************************************************************/
/* pytsai: can fail.  This method is called from within the lmdif_ routine.
 * To indicate failure, set *iflag = -1. */
void cc_five_parm_optimization_with_early_distortion_removal_error (m_ptr, n_ptr, params, err, iflag)
    integer  *m_ptr;		/* pointer to number of points to fit */
    integer  *n_ptr;		/* pointer to number of parameters */
    doublereal *params;		/* vector of parameters */
    doublereal *err;		/* vector of error from data */
    integer  *iflag;            /* flag to indicate error to caller */
{
    int       i;

    double    f,
              Tz,
              xc,
              yc,
              zc,
              Xu_1,
              Yu_1,
              Xu_2,
              Yu_2;

    /* in this routine radial lens distortion is taken into account */
    /* before the rotation and translation constants are determined */
    /* (this assumes we have the distortion reasonably modelled)    */

    f = params[0];
    Tz = params[1];

    tsai_cc.kappa1 = params[2];
    tsai_cp.Cx = params[3];
    tsai_cp.Cy = params[4];

    cc_compute_Xd_Yd_and_r_squared ();

    /* remove the sensor distortion before computing the translation and rotation stuff */
    cc_remove_sensor_plane_distortion_from_Xd_and_Yd ();

    if (!cc_compute_U())
    {
        *iflag = -1;
        return;
    }

    cc_compute_Tx_and_Ty ();

    cc_compute_R ();

    /* we need to do this just to see if we have to flip the rotation matrix */
    if (!cc_compute_approximate_f_and_Tz())
    {
        *iflag = -1;
        return;
    }

    if (tsai_cc.f < 0) {
        /* try the other solution for the orthonormal matrix */
	tsai_cc.r3 = -tsai_cc.r3;
	tsai_cc.r6 = -tsai_cc.r6;
	tsai_cc.r7 = -tsai_cc.r7;
	tsai_cc.r8 = -tsai_cc.r8;
	solve_RPY_transform ();

        if (!cc_compute_approximate_f_and_Tz())
        {
                *iflag = -1;
                return;
        }

        if (tsai_cc.f < 0) 
        {
            pytsai_raise("error - possible handedness problem with data");
            *iflag = -1;
            return;
	}
    }

    /* now calculate the squared error assuming zero distortion */
    for (i = 0; i < tsai_cd.point_count; i++) {
	/* convert from world coordinates to camera coordinates */
	/* Note: zw is implicitly assumed to be zero for these (coplanar) calculations */
	xc = tsai_cc.r1 * tsai_cd.xw[i] + tsai_cc.r2 * tsai_cd.yw[i] + tsai_cc.Tx;
	yc = tsai_cc.r4 * tsai_cd.xw[i] + tsai_cc.r5 * tsai_cd.yw[i] + tsai_cc.Ty;
	zc = tsai_cc.r7 * tsai_cd.xw[i] + tsai_cc.r8 * tsai_cd.yw[i] + Tz;

	/* convert from camera coordinates to undistorted sensor coordinates */
	Xu_1 = f * xc / zc;
	Yu_1 = f * yc / zc;

	/* convert from distorted sensor coordinates to undistorted sensor coordinates  */
	/* (already done, actually)							 */
	Xu_2 = Xd[i];
	Yu_2 = Yd[i];

        /* record the error in the undistorted sensor coordinates */
        err[i] = hypot (Xu_1 - Xu_2, Yu_1 - Yu_2);
    }

    return;
}

 
/* pytsai: can fail; int return type is required. */
int cc_five_parm_optimization_with_early_distortion_removal ()
{
#define NPARAMS 5

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
    if (( fvec = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       pytsai_raise("malloc: Cannot allocate workspace fvec");
       return 0;
    }
 
    if (( fjac = (doublereal *) malloc ((unsigned int) m*n * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       pytsai_raise("malloc: Cannot allocate workspace fjac");
       return 0;
    }
 
    if (( wa4 = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       free(fjac);
       pytsai_raise("malloc: Cannot allocate workspace wa4");
       return 0;
    }

    /* use the current calibration and camera constants as a starting point */
    x[0] = tsai_cc.f;
    x[1] = tsai_cc.Tz;
    x[2] = tsai_cc.kappa1;
    x[3] = tsai_cp.Cx;
    x[4] = tsai_cp.Cy;

    /* define optional scale factors for the parameters */
    if ( mode == 2 ) {
        for (i = 0; i < NPARAMS; i++)
            diag[i] = 1.0;             /* some user-defined values */
    }
 
    /* perform the optimization */
    lmdif_ (cc_five_parm_optimization_with_early_distortion_removal_error,
            &m, &n, x, fvec, &ftol, &xtol, &gtol, &maxfev, &epsfcn,
            diag, &mode, &factor, &nprint, &info, &nfev, fjac, &ldfjac,
            ipvt, qtf, wa1, wa2, wa3, wa4);
    /* check for pytsai error condition */
    if (info == -1)
    {
        free(fvec);
        free(fjac);
        free(wa4);
        return 0;
    }
    /* TODO: Check for other error conditions.  See lmdif.c for possible
     * values of the info parameter. */

    /* update the calibration and camera constants */
    tsai_cc.f = x[0];
    tsai_cc.Tz = x[1];
    tsai_cc.kappa1 = x[2];
    tsai_cp.Cx = x[3];
    tsai_cp.Cy = x[4];

    /* release allocated workspace */
    free(fvec);
    free(fjac);
    free(wa4);

#ifdef DEBUG
    /* print the number of function calls during iteration */
    fprintf(stderr,"info: %d nfev: %d\n\n",info,nfev);
#endif

    return 1;

#undef NPARAMS
}


/************************************************************************/
/* pytsai: cannot fail; void return type is fine. */
void cc_nic_optimization_error (m_ptr, n_ptr, params, err)
    integer  *m_ptr;		/* pointer to number of points to fit */
    integer  *n_ptr;		/* pointer to number of parameters */
    doublereal *params;		/* vector of parameters */
    doublereal *err;		/* vector of error from data */
{
    int       i;

    double    xc,
              yc,
              zc,
              Xd_,
              Yd_,
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
              kappa1,
              f,
              r1,
              r2,
              r4,
              r5,
              r7,
              r8,
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
    kappa1 = params[6];
    f = params[7];

    SINCOS (Rx, sa, ca);
    SINCOS (Ry, sb, cb);
    SINCOS (Rz, sg, cg);
    r1 = cb * cg;
    r2 = cg * sa * sb - ca * sg;
    r4 = cb * sg;
    r5 = sa * sb * sg + ca * cg;
    r7 = -sb;
    r8 = cb * sa;

    for (i = 0; i < tsai_cd.point_count; i++) {
	/* convert from world coordinates to camera coordinates */
	/* Note: zw is implicitly assumed to be zero for these (coplanar) calculations */
	xc = r1 * tsai_cd.xw[i] + r2 * tsai_cd.yw[i] + Tx;
	yc = r4 * tsai_cd.xw[i] + r5 * tsai_cd.yw[i] + Ty;
	zc = r7 * tsai_cd.xw[i] + r8 * tsai_cd.yw[i] + Tz;

	/* convert from camera coordinates to undistorted sensor plane coordinates */
	Xu_1 = f * xc / zc;
	Yu_1 = f * yc / zc;

	/* convert from 2D image coordinates to distorted sensor coordinates */
	Xd_ = tsai_cp.dpx * (tsai_cd.Xf[i] - tsai_cp.Cx) / tsai_cp.sx; 
	Yd_ = tsai_cp.dpy * (tsai_cd.Yf[i] - tsai_cp.Cy);

	/* convert from distorted sensor coordinates to undistorted sensor plane coordinates */
	distortion_factor = 1 + kappa1 * (SQR (Xd_) + SQR (Yd_));
	Xu_2 = Xd_ * distortion_factor;
	Yu_2 = Yd_ * distortion_factor;

        /* record the error in the undistorted sensor coordinates */
        err[i] = hypot (Xu_1 - Xu_2, Yu_1 - Yu_2);
    }
}


/* pytsai: can fail; int return type is required. */
int cc_nic_optimization ()
{
#define NPARAMS 8

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
    if (( fvec = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       pytsai_raise("malloc: Cannot allocate workspace fvec");
       return 0;
    }

    if (( fjac = (doublereal *) malloc ((unsigned int) m*n * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       pytsai_raise("malloc: Cannot allocate workspace fjac");
       return 0;
    }
 
    if (( wa4 = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       free(fjac);
       pytsai_raise("malloc: Cannot allocate workspace wa4");
       return 0;
    }

    /* use the current calibration and camera constants as a starting point */
    x[0] = tsai_cc.Rx;
    x[1] = tsai_cc.Ry;
    x[2] = tsai_cc.Rz;
    x[3] = tsai_cc.Tx;
    x[4] = tsai_cc.Ty;
    x[5] = tsai_cc.Tz;
    x[6] = tsai_cc.kappa1;
    x[7] = tsai_cc.f;

    /* define optional scale factors for the parameters */
    if ( mode == 2 ) {
        for (i = 0; i < NPARAMS; i++)
            diag[i] = 1.0;             /* some user-defined values */
    }
 
    /* perform the optimization */
    lmdif_ (cc_nic_optimization_error,
            &m, &n, x, fvec, &ftol, &xtol, &gtol, &maxfev, &epsfcn,
            diag, &mode, &factor, &nprint, &info, &nfev, fjac, &ldfjac,
            ipvt, qtf, wa1, wa2, wa3, wa4);
    /* TODO: check for error conditions in info. */

    /* update the calibration and camera constants */
    tsai_cc.Rx = x[0];
    tsai_cc.Ry = x[1];
    tsai_cc.Rz = x[2];
    apply_RPY_transform ();

    tsai_cc.Tx = x[3];
    tsai_cc.Ty = x[4];
    tsai_cc.Tz = x[5];
    tsai_cc.kappa1 = x[6];
    tsai_cc.f = x[7];

    /* release allocated workspace */
    free(fvec);
    free(fjac);
    free(wa4);

#ifdef DEBUG
    /* print the number of function calls during iteration */
    fprintf(stderr,"info: %d nfev: %d\n\n",info,nfev);
#endif

    return 1;

#undef NPARAMS
}


/************************************************************************/
/* pytsai: cannot fail; void return type is fine. */
void cc_full_optimization_error (m_ptr, n_ptr, params, err)
    integer  *m_ptr;		/* pointer to number of points to fit */
    integer  *n_ptr;		/* pointer to number of parameters */
    doublereal *params;		/* vector of parameters */
    doublereal *err;		/* vector of error from data */
{
    int       i;

    double    xc,
              yc,
              zc,
              Xd_,
              Yd_,
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
              kappa1,
              f,
              Cx,
              Cy,
              r1,
              r2,
              r4,
              r5,
              r7,
              r8,
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
    kappa1 = params[6];
    f = params[7];
    Cx = params[8];
    Cy = params[9];

    SINCOS (Rx, sa, ca);
    SINCOS (Ry, sb, cb);
    SINCOS (Rz, sg, cg);
    r1 = cb * cg;
    r2 = cg * sa * sb - ca * sg;
    r4 = cb * sg;
    r5 = sa * sb * sg + ca * cg;
    r7 = -sb;
    r8 = cb * sa;

    for (i = 0; i < tsai_cd.point_count; i++) {
	/* convert from world coordinates to camera coordinates */
	/* Note: zw is implicitly assumed to be zero for these (coplanar) calculations */
	xc = r1 * tsai_cd.xw[i] + r2 * tsai_cd.yw[i] + Tx;
	yc = r4 * tsai_cd.xw[i] + r5 * tsai_cd.yw[i] + Ty;
	zc = r7 * tsai_cd.xw[i] + r8 * tsai_cd.yw[i] + Tz;

	/* convert from camera coordinates to undistorted sensor plane coordinates */
	Xu_1 = f * xc / zc;
	Yu_1 = f * yc / zc;

	/* convert from 2D image coordinates to distorted sensor coordinates */
	Xd_ = tsai_cp.dpx * (tsai_cd.Xf[i] - Cx) / tsai_cp.sx; 
	Yd_ = tsai_cp.dpy * (tsai_cd.Yf[i] - Cy);

	/* convert from distorted sensor coordinates to undistorted sensor plane coordinates */
	distortion_factor = 1 + kappa1 * (SQR (Xd_) + SQR (Yd_));
	Xu_2 = Xd_ * distortion_factor;
	Yu_2 = Yd_ * distortion_factor;

        /* record the error in the undistorted sensor coordinates */
        err[i] = hypot (Xu_1 - Xu_2, Yu_1 - Yu_2);
    }
}


/* pytsai: can fail; int return type is required. */
int cc_full_optimization ()
{
#define NPARAMS 10

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
    if (( fvec = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       pytsai_raise("malloc: Cannot allocate workspace fvec");
       return 0;
    }

    if (( fjac = (doublereal *) malloc ((unsigned int) m*n * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       pytsai_raise("malloc: Cannot allocate workspace fjac");
       return 0;
    }

    if (( wa4 = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       free(fjac);
       pytsai_raise("malloc: Cannot allocate workspace wa4");
       return 0;
    }

    /* use the current calibration and camera constants as a starting point */
    x[0] = tsai_cc.Rx;
    x[1] = tsai_cc.Ry;
    x[2] = tsai_cc.Rz;
    x[3] = tsai_cc.Tx;
    x[4] = tsai_cc.Ty;
    x[5] = tsai_cc.Tz;
    x[6] = tsai_cc.kappa1;
    x[7] = tsai_cc.f;
    x[8] = tsai_cp.Cx;
    x[9] = tsai_cp.Cy;

    /* define optional scale factors for the parameters */
    if ( mode == 2 ) {
        for (i = 0; i < NPARAMS; i++)
            diag[i] = 1.0;             /* some user-defined values */
    }

    /* perform the optimization */
    lmdif_ (cc_full_optimization_error,
            &m, &n, x, fvec, &ftol, &xtol, &gtol, &maxfev, &epsfcn,
            diag, &mode, &factor, &nprint, &info, &nfev, fjac, &ldfjac,
            ipvt, qtf, wa1, wa2, wa3, wa4);
    /* TODO: Check for error conditions (info parameter).  See lmdif.c for
     * possible values. */

    /* update the calibration and camera constants */
    tsai_cc.Rx = x[0];
    tsai_cc.Ry = x[1];
    tsai_cc.Rz = x[2];
    apply_RPY_transform ();

    tsai_cc.Tx = x[3];
    tsai_cc.Ty = x[4];
    tsai_cc.Tz = x[5];
    tsai_cc.kappa1 = x[6];
    tsai_cc.f = x[7];
    tsai_cp.Cx = x[8];
    tsai_cp.Cy = x[9];

    /* release allocated workspace */
    free(fvec);
    free(fjac);
    free(wa4);

#ifdef DEBUG
    /* print the number of function calls during iteration */
    fprintf(stderr,"info: %d nfev: %d\n\n",info,nfev);
#endif

    return 1;

#undef NPARAMS
}


/***********************************************************************\
* Routines for noncoplanar camera calibration	 			*
\***********************************************************************/
/* pytsai: cannot fail; void return type is fine. */
void ncc_compute_Xd_Yd_and_r_squared ()
{
    int       i;

    double    Xd_,
              Yd_;

    for (i = 0; i < tsai_cd.point_count; i++) {
	Xd[i] = Xd_ = tsai_cp.dpx * (tsai_cd.Xf[i] - tsai_cp.Cx) / tsai_cp.sx;      /* [mm] */
	Yd[i] = Yd_ = tsai_cp.dpy * (tsai_cd.Yf[i] - tsai_cp.Cy);              /* [mm] */
	r_squared[i] = SQR (Xd_) + SQR (Yd_);                   /* [mm^2] */
    }
}


/* pytsai: can fail; int return type required. */
int ncc_compute_U ()
{
    int       i;

    dmat      M,
              a,
              b;

    M = newdmat (0, (tsai_cd.point_count - 1), 0, 6, &errno);
    if (errno) {
	pytsai_raise("ncc compute U: unable to allocate matrix M");
	return 0;
    }

    a = newdmat (0, 6, 0, 0, &errno);
    if (errno) {
        freemat(M);
	pytsai_raise("ncc compute U: unable to allocate vector a");
	return 0;
    }

    b = newdmat (0, (tsai_cd.point_count - 1), 0, 0, &errno);
    if (errno) {
        freemat(M);
        freemat(a);
	pytsai_raise("ncc compute U: unable to allocate vector b");
	return 0;
    }

    for (i = 0; i < tsai_cd.point_count; i++) {
	M.el[i][0] = Yd[i] * tsai_cd.xw[i];
	M.el[i][1] = Yd[i] * tsai_cd.yw[i];
	M.el[i][2] = Yd[i] * tsai_cd.zw[i];
	M.el[i][3] = Yd[i];
	M.el[i][4] = -Xd[i] * tsai_cd.xw[i];
	M.el[i][5] = -Xd[i] * tsai_cd.yw[i];
	M.el[i][6] = -Xd[i] * tsai_cd.zw[i];
	b.el[i][0] = Xd[i];
    }

    if (solve_system (M, a, b)) {
        freemat(M);
        freemat(a);
        freemat(b);
	pytsai_raise("ncc compute U: error - non-coplanar calibration tried with data which may possibly be coplanar");
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
void ncc_compute_Tx_and_Ty ()
{
    int       i,
              far_point;

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

    /* first find the square of the magnitude of Ty */
    Ty_squared = 1 / (SQR (U[4]) + SQR (U[5]) + SQR (U[6]));

    /* find a point that is far from the image center */
    far_distance = 0;
    far_point = 0;
    for (i = 0; i < tsai_cd.point_count; i++)
	if ((distance = r_squared[i]) > far_distance) {
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
    if ((SIGNBIT (x) != SIGNBIT (Xd[far_point])) ||
	(SIGNBIT (y) != SIGNBIT (Yd[far_point])))
	Ty = -Ty;

    /* update the calibration constants */
    tsai_cc.Tx = U[3] * Ty;
    tsai_cc.Ty = Ty;
}


/* pytsai: cannot fail; void return type is fine. */
void ncc_compute_sx ()
{
    tsai_cp.sx = sqrt (SQR (U[0]) + SQR (U[1]) + SQR (U[2])) * fabs (tsai_cc.Ty);
}


/* pytsai: cannot fail; void return type is fine. */
void ncc_compute_R ()
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

    r1 = U[0] * tsai_cc.Ty / tsai_cp.sx;
    r2 = U[1] * tsai_cc.Ty / tsai_cp.sx;
    r3 = U[2] * tsai_cc.Ty / tsai_cp.sx;

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


/* pytsai: cannot fail; void return type is fine. */
void ncc_compute_better_R ()
{
    double    r1,
              r2,
              r3,
              r4,
              r5,
              r6,
              r7,
              sa,
              ca,
              sb,
              cb,
              sg,
              cg;

    r1 = U[0] * tsai_cc.Ty / tsai_cp.sx;
    r2 = U[1] * tsai_cc.Ty / tsai_cp.sx;
    r3 = U[2] * tsai_cc.Ty / tsai_cp.sx;

    r4 = U[4] * tsai_cc.Ty;
    r5 = U[5] * tsai_cc.Ty;
    r6 = U[6] * tsai_cc.Ty;

    /* use the outer product of the first two rows to get the last row */
    r7 = r2 * r6 - r3 * r5;

    /* now find the RPY angles corresponding to the estimated rotation matrix */
    tsai_cc.Rz = atan2 (r4, r1);

    SINCOS (tsai_cc.Rz, sg, cg);

    tsai_cc.Ry = atan2 (-r7, r1 * cg + r4 * sg);

    tsai_cc.Rx = atan2 (r3 * sg - r6 * cg, r5 * cg - r2 * sg);

    SINCOS (tsai_cc.Rx, sa, ca);

    SINCOS (tsai_cc.Ry, sb, cb);

    /* now generate a more orthonormal rotation matrix from the RPY angles */
    tsai_cc.r1 = cb * cg;
    tsai_cc.r2 = cg * sa * sb - ca * sg;
    tsai_cc.r3 = sa * sg + ca * cg * sb;
    tsai_cc.r4 = cb * sg;
    tsai_cc.r5 = sa * sb * sg + ca * cg;
    tsai_cc.r6 = ca * sb * sg - cg * sa;
    tsai_cc.r7 = -sb;
    tsai_cc.r8 = cb * sa;
    tsai_cc.r9 = ca * cb;
}


/* pytsai: can fail; int return type required. */
int ncc_compute_approximate_f_and_Tz ()
{
    int       i;

    dmat      M,
              a,
              b;

    M = newdmat (0, (tsai_cd.point_count - 1), 0, 1, &errno);
    if (errno) {
	pytsai_raise("ncc compute apx: unable to allocate matrix M");
	return 0;
    }

    a = newdmat (0, 1, 0, 0, &errno);
    if (errno) {
        freemat(M);
	pytsai_raise("ncc compute apx: unable to allocate vector a");
	return 0;
    }

    b = newdmat (0, (tsai_cd.point_count - 1), 0, 0, &errno);
    if (errno) {
        freemat(M);
        freemat(a);
	pytsai_raise("ncc compute apx: unable to allocate vector b");
	return 0;
    }

    for (i = 0; i < tsai_cd.point_count; i++) {
	M.el[i][0] = tsai_cc.r4 * tsai_cd.xw[i] + tsai_cc.r5 * tsai_cd.yw[i] + tsai_cc.r6 * tsai_cd.zw[i] + tsai_cc.Ty;
	M.el[i][1] = -Yd[i];
	b.el[i][0] = (tsai_cc.r7 * tsai_cd.xw[i] + tsai_cc.r8 * tsai_cd.yw[i] + tsai_cc.r9 * tsai_cd.zw[i]) * Yd[i];
    }

    if (solve_system (M, a, b)) {
        freemat(M);
        freemat(a);
        freemat(b);
	pytsai_raise("ncc compute apx: unable to solve system  Ma=b");
	return 0;
    }

    /* update the calibration constants */
    tsai_cc.f = a.el[0][0];
    tsai_cc.Tz = a.el[1][0];
    tsai_cc.kappa1 = 0.0;		/* this is the assumption that our calculation was made under */

    freemat (M);
    freemat (a);
    freemat (b);

    return 1;
}


/************************************************************************/
/* pytsai: cannot fail; void return type is fine. */
void ncc_compute_exact_f_and_Tz_error (m_ptr, n_ptr, params, err)
    integer  *m_ptr;		/* pointer to number of points to fit */
    integer  *n_ptr;		/* pointer to number of parameters */
    doublereal *params;		/* vector of parameters */
    doublereal *err;		/* vector of error from data */
{
    int       i;

    double    xc,
              yc,
              zc,
              Xu_1,
              Yu_1,
              Xu_2,
              Yu_2,
              distortion_factor,
              f,
              Tz,
              kappa1;

    f = params[0];
    Tz = params[1];
    kappa1 = params[2];

    for (i = 0; i < tsai_cd.point_count; i++) {
	/* convert from world coordinates to camera coordinates */
	xc = tsai_cc.r1 * tsai_cd.xw[i] + tsai_cc.r2 * tsai_cd.yw[i] + tsai_cc.r3 * tsai_cd.zw[i] + tsai_cc.Tx;
	yc = tsai_cc.r4 * tsai_cd.xw[i] + tsai_cc.r5 * tsai_cd.yw[i] + tsai_cc.r6 * tsai_cd.zw[i] + tsai_cc.Ty;
	zc = tsai_cc.r7 * tsai_cd.xw[i] + tsai_cc.r8 * tsai_cd.yw[i] + tsai_cc.r9 * tsai_cd.zw[i] + Tz;

	/* convert from camera coordinates to undistorted sensor coordinates */
	Xu_1 = f * xc / zc;
	Yu_1 = f * yc / zc;

	/* convert from distorted sensor coordinates to undistorted sensor coordinates */
	distortion_factor = 1 + kappa1 * r_squared[i];
	Xu_2 = Xd[i] * distortion_factor;
	Yu_2 = Yd[i] * distortion_factor;

        /* record the error in the undistorted sensor coordinates */
        err[i] = hypot (Xu_1 - Xu_2, Yu_1 - Yu_2);
    }
}


/* pytsai: can fail; int return type is required. */
int ncc_compute_exact_f_and_Tz ()
{
#define NPARAMS 3

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
    if (( fvec = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       pytsai_raise("malloc: Cannot allocate workspace fvec");
       return 0;
    }

    if (( fjac = (doublereal *) malloc ((unsigned int) m*n * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       pytsai_raise("malloc: Cannot allocate workspace fjac");
       return 0;
    }

    if (( wa4 = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       free(fjac);
       pytsai_raise("malloc: Cannot allocate workspace wa4");
       return 0;
    }

    /* use the current calibration constants as an initial guess */
    x[0] = tsai_cc.f;
    x[1] = tsai_cc.Tz;
    x[2] = tsai_cc.kappa1;

    /* define optional scale factors for the parameters */
    if ( mode == 2 ) {
        for (i = 0; i < NPARAMS; i++)
            diag[i] = 1.0;             /* some user-defined values */
    }
 
    /* perform the optimization */
    lmdif_ (ncc_compute_exact_f_and_Tz_error,
            &m, &n, x, fvec, &ftol, &xtol, &gtol, &maxfev, &epsfcn,
            diag, &mode, &factor, &nprint, &info, &nfev, fjac, &ldfjac,
            ipvt, qtf, wa1, wa2, wa3, wa4);
    /* TODO: Check for error conditions (info parameter).  See lmbif.c
     * for possible values. */

    /* update the calibration constants */
    tsai_cc.f = x[0];
    tsai_cc.Tz = x[1];
    tsai_cc.kappa1 = x[2];

    /* release allocated workspace */
    free(fvec);
    free(fjac);
    free(wa4);

#ifdef DEBUG
    /* print the number of function calls during iteration */
    fprintf(stderr,"info: %d nfev: %d\n\n",info,nfev);
#endif

    return 1;

#undef NPARAMS
}


/************************************************************************/
/* pytsai: can fail; int return type is required. */
int ncc_three_parm_optimization ()
{
    ncc_compute_Xd_Yd_and_r_squared ();

    if (!ncc_compute_U())
        return 0;

    ncc_compute_Tx_and_Ty ();

    ncc_compute_sx ();

    ncc_compute_Xd_Yd_and_r_squared ();

    ncc_compute_better_R ();

    if (!ncc_compute_approximate_f_and_Tz())
        return 0;

    if (tsai_cc.f < 0) {
	/* try the other solution for the orthonormal matrix */
	tsai_cc.r3 = -tsai_cc.r3;
	tsai_cc.r6 = -tsai_cc.r6;
	tsai_cc.r7 = -tsai_cc.r7;
	tsai_cc.r8 = -tsai_cc.r8;
	solve_RPY_transform ();

	if (!ncc_compute_approximate_f_and_Tz())
                return 0;

        if (tsai_cc.f < 0) {
            pytsai_raise("error - possible handedness problem with data");
            return 0;
	}
    }

    if (!ncc_compute_exact_f_and_Tz())
        return 0;

    return 1;
}


/************************************************************************/
/* pytsai: cannot fail; void return type is fine. */
void ncc_nic_optimization_error (m_ptr, n_ptr, params, err)
    integer  *m_ptr;		/* pointer to number of points to fit */
    integer  *n_ptr;		/* pointer to number of parameters */
    doublereal *params;		/* vector of parameters */
    doublereal *err;		/* vector of error from data */
{
    int       i;

    double    xc,
              yc,
              zc,
              Xd_,
              Yd_,
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
              kappa1,
              sx,
              f,
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
    kappa1 = params[6];
    f = params[7];
    sx = params[8];

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
	Xu_1 = f * xc / zc;
	Yu_1 = f * yc / zc;

	/* convert from 2D image coordinates to distorted sensor coordinates */
	Xd_ = tsai_cp.dpx * (tsai_cd.Xf[i] - tsai_cp.Cx) / sx;
	Yd_ = tsai_cp.dpy * (tsai_cd.Yf[i] - tsai_cp.Cy);

	/* convert from distorted sensor coordinates to undistorted sensor plane coordinates */
	distortion_factor = 1 + kappa1 * (SQR (Xd_) + SQR (Yd_));
	Xu_2 = Xd_ * distortion_factor;
	Yu_2 = Yd_ * distortion_factor;

        /* record the error in the undistorted sensor coordinates */
        err[i] = hypot (Xu_1 - Xu_2, Yu_1 - Yu_2);
    }
}


/* pytsai: can fail; int return type is required. */
int ncc_nic_optimization ()
{
#define NPARAMS 9

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
    if (( fvec = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       pytsai_raise("malloc: Cannot allocate workspace fvec");
       return 0;
    }

    if (( fjac = (doublereal *) malloc ((unsigned int) m*n * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       pytsai_raise("malloc: Cannot allocate workspace fjac");
       return 0;
    }

    if (( wa4 = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       free(fjac);
       pytsai_raise("malloc: Cannot allocate workspace wa4");
       return 0;
    }

    /* use the current calibration and camera constants as a starting point */
    x[0] = tsai_cc.Rx;
    x[1] = tsai_cc.Ry;
    x[2] = tsai_cc.Rz;
    x[3] = tsai_cc.Tx;
    x[4] = tsai_cc.Ty;
    x[5] = tsai_cc.Tz;
    x[6] = tsai_cc.kappa1;
    x[7] = tsai_cc.f;
    x[8] = tsai_cp.sx;

    /* define optional scale factors for the parameters */
    if ( mode == 2 ) {
        for (i = 0; i < NPARAMS; i++)
            diag[i] = 1.0;             /* some user-defined values */
    }

    /* perform the optimization */
    lmdif_ (ncc_nic_optimization_error,
            &m, &n, x, fvec, &ftol, &xtol, &gtol, &maxfev, &epsfcn,
            diag, &mode, &factor, &nprint, &info, &nfev, fjac, &ldfjac,
            ipvt, qtf, wa1, wa2, wa3, wa4);
    /* TODO: Check for error conditions (info parameter).  See lmdif.c for
     * possible values of the info parameter. */

    /* update the calibration and camera constants */
    tsai_cc.Rx = x[0];
    tsai_cc.Ry = x[1];
    tsai_cc.Rz = x[2];
    apply_RPY_transform ();

    tsai_cc.Tx = x[3];
    tsai_cc.Ty = x[4];
    tsai_cc.Tz = x[5];
    tsai_cc.kappa1 = x[6];
    tsai_cc.f = x[7];
    tsai_cp.sx = x[8];

    /* release allocated workspace */
    free(fvec);
    free(fjac);
    free(wa4);

#ifdef DEBUG
    /* print the number of function calls during iteration */
    fprintf(stderr,"info: %d nfev: %d\n\n",info,nfev);
#endif

    return 1;

#undef NPARAMS
}


/************************************************************************/
/* pytsai: cannot fail; void return type is fine. */
void ncc_full_optimization_error (m_ptr, n_ptr, params, err)
    integer  *m_ptr;		/* pointer to number of points to fit */
    integer  *n_ptr;		/* pointer to number of parameters */
    doublereal *params;		/* vector of parameters */
    doublereal *err;		/* vector of error from data */
{
    int       i;

    double    xc,
              yc,
              zc,
              Xd_,
              Yd_,
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
              kappa1,
              sx,
              f,
              Cx,
              Cy,
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
    kappa1 = params[6];
    f = params[7];
    sx = params[8];
    Cx = params[9];
    Cy = params[10];

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
	Xu_1 = f * xc / zc;
	Yu_1 = f * yc / zc;

	/* convert from 2D image coordinates to distorted sensor coordinates */
	Xd_ = tsai_cp.dpx * (tsai_cd.Xf[i] - Cx) / sx;
	Yd_ = tsai_cp.dpy * (tsai_cd.Yf[i] - Cy);

	/* convert from distorted sensor coordinates to undistorted sensor plane coordinates */
	distortion_factor = 1 + kappa1 * (SQR (Xd_) + SQR (Yd_));
	Xu_2 = Xd_ * distortion_factor;
	Yu_2 = Yd_ * distortion_factor;

        /* record the error in the undistorted sensor coordinates */
        err[i] = hypot (Xu_1 - Xu_2, Yu_1 - Yu_2);
    }
}


/* pytsai: can fail; int return type is required. */
int ncc_full_optimization ()
{
#define NPARAMS 11

    int     i;

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
    if (( fvec = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       pytsai_raise("malloc: Cannot allocate workspace fvec");
       return 0;
    }
 
    if (( fjac = (doublereal *) malloc ((unsigned int) m*n * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       pytsai_raise("malloc: Cannot allocate workspace fjac");
       return 0;
    }
 
    if (( wa4 = (doublereal *) malloc ((unsigned int) m * (unsigned int) sizeof(doublereal))) == NULL ) {
       free(fvec);
       free(fjac);
       pytsai_raise("malloc: Cannot allocate workspace wa4");
       return 0;
    }

    /* use the current calibration and camera constants as a starting point */
    x[0] = tsai_cc.Rx;
    x[1] = tsai_cc.Ry;
    x[2] = tsai_cc.Rz;
    x[3] = tsai_cc.Tx;
    x[4] = tsai_cc.Ty;
    x[5] = tsai_cc.Tz;
    x[6] = tsai_cc.kappa1;
    x[7] = tsai_cc.f;
    x[8] = tsai_cp.sx;
    x[9] = tsai_cp.Cx;
    x[10] = tsai_cp.Cy;

    /* define optional scale factors for the parameters */
    if ( mode == 2 ) {
        for (i = 0; i < NPARAMS; i++)
            diag[i] = 1.0;             /* some user-defined values */
    }
 
    /* perform the optimization */
    lmdif_ (ncc_full_optimization_error,
            &m, &n, x, fvec, &ftol, &xtol, &gtol, &maxfev, &epsfcn,
            diag, &mode, &factor, &nprint, &info, &nfev, fjac, &ldfjac,
            ipvt, qtf, wa1, wa2, wa3, wa4);
    /* TODO: Check for error conditions (info parameter).  See lmdif.c for
     * possible values. */

    /* update the calibration and camera constants */
    tsai_cc.Rx = x[0];
    tsai_cc.Ry = x[1];
    tsai_cc.Rz = x[2];
    apply_RPY_transform ();

    tsai_cc.Tx = x[3];
    tsai_cc.Ty = x[4];
    tsai_cc.Tz = x[5];
    tsai_cc.kappa1 = x[6];
    tsai_cc.f = x[7];
    tsai_cp.sx = x[8];
    tsai_cp.Cx = x[9];
    tsai_cp.Cy = x[10];

    /* release allocated workspace */
    free(fvec);
    free(fjac);
    free(wa4);

#ifdef DEBUG
    /* print the number of function calls during iteration */
    fprintf(stderr,"info: %d nfev: %d\n\n",info,nfev);
#endif

    return 1;

#undef NPARAMS
}


/************************************************************************/

/* pytsai: can fail; int return type is required. */
int coplanar_calibration ()
{
    /* just do the basic 3 parameter (Tz, f, kappa1) optimization */
    return cc_three_parm_optimization ();
}

 
/* pytsai: can fail; int return type is required. */
int coplanar_calibration_with_full_optimization ()
{
    /* start with a 3 parameter (Tz, f, kappa1) optimization */
    if (!cc_three_parm_optimization())
        return 0;

    /* do a 5 parameter (Tz, f, kappa1, Cx, Cy) optimization */
    if (!cc_five_parm_optimization_with_late_distortion_removal())
        return 0;

    /* do a better 5 parameter (Tz, f, kappa1, Cx, Cy) optimization */
    if (!cc_five_parm_optimization_with_early_distortion_removal())
        return 0;

    /* do a full optimization minus the image center */
    if (!cc_nic_optimization())
        return 0;

    /* do a full optimization including the image center */
    if (!cc_full_optimization())
        return 0;

    return 1;
}


/* pytsai: can fail; int return type is required. */
int noncoplanar_calibration ()
{
    /* just do the basic 3 parameter (Tz, f, kappa1) optimization */
    return ncc_three_parm_optimization ();
}

 
/* pytsai: can fail; int return type is required. */
int noncoplanar_calibration_with_full_optimization ()
{
    /* start with a 3 parameter (Tz, f, kappa1) optimization */
    if (!ncc_three_parm_optimization())
        return 0;

    /* do a full optimization minus the image center */
    if (!ncc_nic_optimization())
        return 0;

    /* do a full optimization including the image center */
    if (!ncc_full_optimization())
        return 0;

    return 1;
}
