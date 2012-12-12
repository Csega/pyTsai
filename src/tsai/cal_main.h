/**
 * cal_main.h
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
* This file contains operating constants, data structure definitions and I/O *
* variable declarations for the calibration routines contained in the file   *
* cal_main.c                                                                 *
*                                                                            *
* An explanation of the basic algorithms and description of the variables    *
* can be found in "An Efficient and Accurate Camera Calibration Technique    *
* for 3D Machine Vision", Roger Y. Tsai, Proceedings of IEEE Conference on   *
* Computer Vision and Pattern Recognition, 1986, pages 364-374.  This work   *
* also appears in several other publications under Roger Tsai's name.        *
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
*                                                                            *
* History                                                                    *
* -------                                                                    *
*                                                                            *
* 03-Nov-2004 Jonathan Merritt (j.merritt@pgrad.unimelb.edu.au) at           *
*             The University of Melbourne.                                   *
*                                                                            *
* 01-Apr-95  Reg Willson (rgwillson@mmm.com) at 3M St. Paul, MN              *
*       Filename changes for DOS port.                                       *
*                                                                            *
* 04-Jun-94  Reg Willson (rgwillson@mmm.com) at 3M St. Paul, MN              *
*       Added alternate macro definitions for less common math functions.    *
*                                                                            *
* 25-Mar-94  Torfi Thorhallsson (torfit@verk.hi.is) at the University of     *
*            Iceland                                                         *
*       Added definitions of parameters controlling MINPACK's lmdif()        *
*       optimization routine.                                                *
*                                                                            *
* 30-Nov-93  Reg Willson (rgw@cs.cmu.edu) at Carnegie-Mellon University      *
*       Added declaration for camera type.                                   *
*                                                                            *
* 07-Feb-93  Reg Willson (rgw@cs.cmu.edu) at Carnegie-Mellon University      *
*       Original implementation.                                             *
*                                                                            *
*                                                                            *
\****************************************************************************/

#ifndef CAL_MAIN_H
#define CAL_MAIN_H

/* Maximum number of data points allowed */
#define MAX_POINTS	500

/* An arbitrary tolerance factor */
#define EPSILON		1.0E-8

/* Commonly used macros */
#define ABS(a)          (((a) > 0) ? (a) : -(a))
#define SIGNBIT(a)      (((a) > 0) ? 0 : 1)
#define SQR(a)          ((a) * (a))
#define CUB(a)          ((a)*(a)*(a))
#define MAX(a,b)        (((a) > (b)) ? (a) : (b))
#define MIN(a,b)        (((a) < (b)) ? (a) : (b))

/* Some math libraries may not have the sincos() routine */
#ifdef _SINCOS_
void sincos();
#define SINCOS(x,s,c)   sincos(x,&s,&c)
#else
double sin(), cos();
#define SINCOS(x,s,c)   s=sin(x);c=cos(x)
#endif


/* Parameters controlling MINPACK's lmdif() optimization routine. */
/* See the file lmdif.f for definitions of each parameter.        */
#define REL_SENSOR_TOLERANCE_ftol    1.0E-5           /* [pix] */
#define REL_PARAM_TOLERANCE_xtol     1.0E-7
#define ORTHO_TOLERANCE_gtol         0.0
#define MAXFEV                       (1000*n)
#define EPSFCN                       1.0E-16          /* Do not set to zero! */
#define MODE                         1   /* variables are scalled internally */
#define FACTOR                       100.0 

/****************************************************************************\
*                                                                            *
* Camera parameters are usually the fixed parameters of the given camera     *
* system, typically obtained from manufacturers specifications.              *
*                                                                            *
* Cy and Cy (the center of radial lens distortion), may be calibrated        *
* separately or as part of the coplanar/noncoplanar calibration.             *
* The same with sx (x scale uncertainty factor).                             *
*                                                                            *
\****************************************************************************/
struct camera_parameters {
    double    Ncx;		/* [sel]     Number of sensor elements in camera's x direction */
    double    Nfx;		/* [pix]     Number of pixels in frame grabber's x direction   */
    double    dx;		/* [mm/sel]  X dimension of camera's sensor element (in mm)    */
    double    dy;		/* [mm/sel]  Y dimension of camera's sensor element (in mm)    */
    double    dpx;		/* [mm/pix]  Effective X dimension of pixel in frame grabber   */
    double    dpy;		/* [mm/pix]  Effective Y dimension of pixel in frame grabber   */
    double    Cx;		/* [pix]     Z axis intercept of camera coordinate system      */
    double    Cy;		/* [pix]     Z axis intercept of camera coordinate system      */
    double    sx;		/* []        Scale factor to compensate for any error in dpx   */
};

/****************************************************************************\
*                                                                            *
* Calibration data consists of the (x,y,z) world coordinates of a feature    *
* point	(in mm) and the corresponding coordinates (Xf,Yf) (in pixels) of the *
* feature point in the image.  Point count is the number of points in the    *
* data set.                                                                  *
*                                                                            *
*                                                                            *
* Coplanar calibration:                                                      *
*                                                                            *
* For coplanar calibration the z world coordinate of the data must be zero.  *
* In addition, for numerical stability the coordinates of the feature points *
* should be placed at some arbitrary distance from the origin (0,0,0) of the *
* world coordinate system.  Finally, the plane containing the feature points *
* should not be parallel to the imaging plane.  A relative angle of 30       *
* degrees is recommended.                                                    *
*                                                                            *
*                                                                            *
* Noncoplanar calibration:                                                   *
*                                                                            *
* For noncoplanar calibration the data must not lie in a single plane.       *
*                                                                            *
\****************************************************************************/
struct calibration_data {
    int       point_count;	/* [points] 	 */
    double    xw[MAX_POINTS];	/* [mm]          */
    double    yw[MAX_POINTS];	/* [mm]          */
    double    zw[MAX_POINTS];	/* [mm]          */
    double    Xf[MAX_POINTS];	/* [pix]         */
    double    Yf[MAX_POINTS];	/* [pix]         */
};


/****************************************************************************\
*                                                                            *
* Calibration constants are the model constants that are determined from the *
* calibration data.                                                          *
*                                                                            *
\****************************************************************************/
struct calibration_constants {
    double    f;		/* [mm]          */
    double    kappa1;		/* [1/mm^2]      */
    double    p1;		/* [1/mm]        */
    double    p2;		/* [1/mm]        */
    double    Tx;		/* [mm]          */
    double    Ty;		/* [mm]          */
    double    Tz;		/* [mm]          */
    double    Rx;		/* [rad]	 */
    double    Ry;		/* [rad]	 */
    double    Rz;		/* [rad]	 */
    double    r1;		/* []            */
    double    r2;		/* []            */
    double    r3;		/* []            */
    double    r4;		/* []            */
    double    r5;		/* []            */
    double    r6;		/* []            */
    double    r7;		/* []            */
    double    r8;		/* []            */
    double    r9;		/* []            */
};

/* External declarations for variables used by the subroutines for I/O */
extern struct camera_parameters tsai_cp;
extern struct calibration_data tsai_cd;
extern struct calibration_constants tsai_cc;
//extern char   camera_type[];

/* Forward declarations for the calibration routines */
#if 0   /* not needed for the Python version. */
void  initialize_photometrics_parms ();
void  initialize_general_imaging_mos5300_matrox_parms ();
void  initialize_panasonic_gp_mf702_matrox_parms ();
void  initialize_sony_xc77_matrox_parms ();
void  initialize_sony_xc57_androx_parms ();
#endif

int   coplanar_calibration ();
int   coplanar_calibration_with_full_optimization ();

int   noncoplanar_calibration ();
int   noncoplanar_calibration_with_full_optimization ();

int   coplanar_extrinsic_parameter_estimation ();
int   noncoplanar_extrinsic_parameter_estimation ();

void  world_coord_to_image_coord ();
void  image_coord_to_world_coord ();
void  world_coord_to_camera_coord ();
void  camera_coord_to_world_coord ();
void  distorted_to_undistorted_sensor_coord ();
void  undistorted_to_distorted_sensor_coord ();
void  distorted_to_undistorted_image_coord ();
void  undistorted_to_distorted_image_coord ();

void  distorted_image_plane_error_stats ();
void  undistorted_image_plane_error_stats ();
void  object_space_error_stats ();
void  normalized_calibration_error ();

void  solve_RPY_transform ();
void  apply_RPY_transform ();

#endif /* CAL_MAIN_H */

