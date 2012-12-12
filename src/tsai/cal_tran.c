/**
 * cal_tran.c
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
* This file contains routines for transforming between the different         *
* coordinate systems used in Tsai's perspective projection camera model.     *
* The routines are:                                                          *
*                                                                            *
*       world_coord_to_image_coord ()                                        *
*       image_coord_to_world_coord ()                                        *
*       world_coord_to_camera_coord ()                                       *
*       camera_coord_to_world_coord ()                                       *
*       distorted_to_undistorted_sensor_coord ()                             *
*       undistorted_to_distorted_sensor_coord ()                             *
*       distorted_to_undistorted_image_coord ()                              *
*       undistorted_to_distorted_image_coord ()                              *
*                                                                            *
* The routines make use of the calibrated camera parameters and calibration  *
* constants contained in the two external data structures tsai_cp and tsai_cc.         *
*                                                                            *
* Notation                                                                   *
* --------                                                                   *
*                                                                            *
* The camera's X axis runs along increasing column coordinates in the        *
* image/frame.  The Y axis runs along increasing row coordinates.            *
* All 3D coordinates are right-handed.                                       *
*                                                                            *
* pix == image/frame grabber picture element                                 *
* sel == camera sensor element                                               *
*                                                                            *
*                                                                            *
* History                                                                    *
* -------                                                                    *
*                                                                            *
* 18-Oct-95  Reg Willson (rgwillson@mmm.com) at 3M St. Paul, MN              *
*       Added check in undistorted_to_distorted_sensor_coord () for          *
*       situation where an undistorted sensor point maps to the maximum      *
*       barrel distortion radius.                                            *
*                                                                            *
* 18-May-95  Reg Willson (rgwillson@mmm.com) at 3M St. Paul, MN              *
*       Split out from the cal_main.c file.                                  *
*       Fix the CBRT routine so it handles negative arguments properly.      *
*                                                                            *
\****************************************************************************/
#include <stdio.h>
#include <math.h>
#include "cal_main.h"


extern struct camera_parameters tsai_cp;
extern struct calibration_constants tsai_cc;


#define SQRT(x) sqrt(fabs(x))


/************************************************************************/
/* This cube root routine handles negative arguments (unlike cbrt).     */

double    CBRT (x)
    double    x;
{
    double    pow ();

    if (x == 0)
	return (0);
    else if (x > 0)
	return (pow (x, (double) 1.0 / 3.0));
    else
	return (-pow (-x, (double) 1.0 / 3.0));
} 


/************************************************************************/
/*
       This routine converts from undistorted to distorted sensor coordinates.
       The technique involves algebraically solving the cubic polynomial

            Ru = Rd * (1 + kappa1 * Rd**2)

       using the Cardan method.

       Note: for kappa1 < 0 the distorted sensor plane extends out to a 
             maximum barrel distortion radius of  Rd = sqrt (-1/(3 * kappa1)).

	     To see the logic used in this routine try graphing the above 
	     polynomial for positive and negative kappa1's
*/

void      undistorted_to_distorted_sensor_coord (Xu, Yu, Xd, Yd)
    double    Xu,
              Yu,
             *Xd,
             *Yd;
{
#define SQRT3   1.732050807568877293527446341505872366943

    double    Ru,
              Rd,
              lambda,
              c,
              d,
              Q,
              R,
              D,
              S,
              T,
              sinT,
              cosT;

    if (((Xu == 0) && (Yu == 0)) || (tsai_cc.kappa1 == 0)) {
	*Xd = Xu;
	*Yd = Yu;
	return;
    }

    Ru = hypot (Xu, Yu);	/* SQRT(Xu*Xu+Yu*Yu) */

    c = 1 / tsai_cc.kappa1;
    d = -c * Ru;

    Q = c / 3;
    R = -d / 2;
    D = CUB (Q) + SQR (R);

    if (D >= 0) {		/* one real root */
	D = SQRT (D);
	S = CBRT (R + D);
	T = CBRT (R - D);
	Rd = S + T;

	if (Rd < 0) {
	    Rd = SQRT (-1 / (3 * tsai_cc.kappa1));
	    fprintf (stderr, "\nWarning: undistorted image point to distorted image point mapping limited by\n");
	    fprintf (stderr, "         maximum barrel distortion radius of %lf\n", Rd);
	    fprintf (stderr, "         (Xu = %lf, Yu = %lf) -> (Xd = %lf, Yd = %lf)\n\n",
		     Xu, Yu, Xu * Rd / Ru, Yu * Rd / Ru);
	}
    } else {			/* three real roots */
	D = SQRT (-D);
	S = CBRT (hypot (R, D));
	T = atan2 (D, R) / 3;
	SINCOS (T, sinT, cosT);

	/* the larger positive root is    2*S*cos(T)                   */
	/* the smaller positive root is   -S*cos(T) + SQRT(3)*S*sin(T) */
	/* the negative root is           -S*cos(T) - SQRT(3)*S*sin(T) */

	Rd = -S * cosT + SQRT3 * S * sinT;	/* use the smaller positive root */
    }

    lambda = Rd / Ru;

    *Xd = Xu * lambda;
    *Yd = Yu * lambda;
}


/************************************************************************/
void      distorted_to_undistorted_sensor_coord (Xd, Yd, Xu, Yu)
    double    Xd,
              Yd,
             *Xu,
             *Yu;
{
    double    distortion_factor;

    /* convert from distorted to undistorted sensor plane coordinates */
    distortion_factor = 1 + tsai_cc.kappa1 * (SQR (Xd) + SQR (Yd));
    *Xu = Xd * distortion_factor;
    *Yu = Yd * distortion_factor;
}


/************************************************************************/
void      undistorted_to_distorted_image_coord (Xfu, Yfu, Xfd, Yfd)
    double    Xfu,
              Yfu,
             *Xfd,
             *Yfd;
{
    double    Xu,
              Yu,
              Xd,
              Yd;

    /* convert from image to sensor coordinates */
    Xu = tsai_cp.dpx * (Xfu - tsai_cp.Cx) / tsai_cp.sx;
    Yu = tsai_cp.dpy * (Yfu - tsai_cp.Cy);

    /* convert from undistorted sensor to distorted sensor plane coordinates */
    undistorted_to_distorted_sensor_coord (Xu, Yu, &Xd, &Yd);

    /* convert from sensor to image coordinates */
    *Xfd = Xd * tsai_cp.sx / tsai_cp.dpx + tsai_cp.Cx;
    *Yfd = Yd / tsai_cp.dpy + tsai_cp.Cy;
}


/************************************************************************/
void      distorted_to_undistorted_image_coord (Xfd, Yfd, Xfu, Yfu)
    double    Xfd,
              Yfd,
             *Xfu,
             *Yfu;
{
    double    Xd,
              Yd,
              Xu,
              Yu;

    /* convert from image to sensor coordinates */
    Xd = tsai_cp.dpx * (Xfd - tsai_cp.Cx) / tsai_cp.sx;
    Yd = tsai_cp.dpy * (Yfd - tsai_cp.Cy);

    /* convert from distorted sensor to undistorted sensor plane coordinates */
    distorted_to_undistorted_sensor_coord (Xd, Yd, &Xu, &Yu);

    /* convert from sensor to image coordinates */
    *Xfu = Xu * tsai_cp.sx / tsai_cp.dpx + tsai_cp.Cx;
    *Yfu = Yu / tsai_cp.dpy + tsai_cp.Cy;
}


/***********************************************************************\
* This routine takes the position of a point in world coordinates [mm]	*
* and determines the position of its image in image coordinates [pix].	*
\***********************************************************************/
void      world_coord_to_image_coord (xw, yw, zw, Xf, Yf)
    double    xw,
              yw,
              zw,
             *Xf,
             *Yf;
{
    double    xc,
              yc,
              zc,
              Xu,
              Yu,
              Xd,
              Yd;

    /* convert from world coordinates to camera coordinates */
    xc = tsai_cc.r1 * xw + tsai_cc.r2 * yw + tsai_cc.r3 * zw + tsai_cc.Tx;
    yc = tsai_cc.r4 * xw + tsai_cc.r5 * yw + tsai_cc.r6 * zw + tsai_cc.Ty;
    zc = tsai_cc.r7 * xw + tsai_cc.r8 * yw + tsai_cc.r9 * zw + tsai_cc.Tz;

    /* convert from camera coordinates to undistorted sensor plane coordinates */
    Xu = tsai_cc.f * xc / zc;
    Yu = tsai_cc.f * yc / zc;

    /* convert from undistorted to distorted sensor plane coordinates */
    undistorted_to_distorted_sensor_coord (Xu, Yu, &Xd, &Yd);

    /* convert from distorted sensor plane coordinates to image coordinates */
    *Xf = Xd * tsai_cp.sx / tsai_cp.dpx + tsai_cp.Cx;
    *Yf = Yd / tsai_cp.dpy + tsai_cp.Cy;
}


/***********************************************************************\
* This routine performs an inverse perspective projection to determine	*
* the position of a point in world coordinates that corresponds to a 	*
* given position in image coordinates.  To constrain the inverse	*
* projection to a single point the routine requires a Z world	 	*
* coordinate for the point in addition to the X and Y image coordinates.* 
\***********************************************************************/
void      image_coord_to_world_coord (Xfd, Yfd, zw, xw, yw)
    double    Xfd,
              Yfd, 
              zw,
             *xw,
             *yw;
{
    double    Xd,
              Yd,
              Xu,
              Yu,
              common_denominator;

    /* convert from image to distorted sensor coordinates */
    Xd = tsai_cp.dpx * (Xfd - tsai_cp.Cx) / tsai_cp.sx;
    Yd = tsai_cp.dpy * (Yfd - tsai_cp.Cy);

    /* convert from distorted sensor to undistorted sensor plane coordinates */
    distorted_to_undistorted_sensor_coord (Xd, Yd, &Xu, &Yu);

    /* calculate the corresponding xw and yw world coordinates	 */
    /* (these equations were derived by simply inverting	 */
    /* the perspective projection equations using Macsyma)	 */
    common_denominator = ((tsai_cc.r1 * tsai_cc.r8 - tsai_cc.r2 * tsai_cc.r7) * Yu +
			  (tsai_cc.r5 * tsai_cc.r7 - tsai_cc.r4 * tsai_cc.r8) * Xu -
			  tsai_cc.f * tsai_cc.r1 * tsai_cc.r5 + tsai_cc.f * tsai_cc.r2 * tsai_cc.r4);

    *xw = (((tsai_cc.r2 * tsai_cc.r9 - tsai_cc.r3 * tsai_cc.r8) * Yu +
	    (tsai_cc.r6 * tsai_cc.r8 - tsai_cc.r5 * tsai_cc.r9) * Xu -
	    tsai_cc.f * tsai_cc.r2 * tsai_cc.r6 + tsai_cc.f * tsai_cc.r3 * tsai_cc.r5) * zw +
	   (tsai_cc.r2 * tsai_cc.Tz - tsai_cc.r8 * tsai_cc.Tx) * Yu +
	   (tsai_cc.r8 * tsai_cc.Ty - tsai_cc.r5 * tsai_cc.Tz) * Xu -
	   tsai_cc.f * tsai_cc.r2 * tsai_cc.Ty + tsai_cc.f * tsai_cc.r5 * tsai_cc.Tx) / common_denominator;

    *yw = -(((tsai_cc.r1 * tsai_cc.r9 - tsai_cc.r3 * tsai_cc.r7) * Yu +
	     (tsai_cc.r6 * tsai_cc.r7 - tsai_cc.r4 * tsai_cc.r9) * Xu -
	     tsai_cc.f * tsai_cc.r1 * tsai_cc.r6 + tsai_cc.f * tsai_cc.r3 * tsai_cc.r4) * zw +
	    (tsai_cc.r1 * tsai_cc.Tz - tsai_cc.r7 * tsai_cc.Tx) * Yu +
	    (tsai_cc.r7 * tsai_cc.Ty - tsai_cc.r4 * tsai_cc.Tz) * Xu -
	    tsai_cc.f * tsai_cc.r1 * tsai_cc.Ty + tsai_cc.f * tsai_cc.r4 * tsai_cc.Tx) / common_denominator;
}


/***********************************************************************\
* This routine takes the position of a point in world coordinates [mm]	*
* and determines its position in camera coordinates [mm].		*
\***********************************************************************/
void      world_coord_to_camera_coord (xw, yw, zw, xc, yc, zc)
    double    xw,
              yw,
              zw,
             *xc,
             *yc,
	     *zc;
{
    *xc = tsai_cc.r1 * xw + tsai_cc.r2 * yw + tsai_cc.r3 * zw + tsai_cc.Tx;
    *yc = tsai_cc.r4 * xw + tsai_cc.r5 * yw + tsai_cc.r6 * zw + tsai_cc.Ty;
    *zc = tsai_cc.r7 * xw + tsai_cc.r8 * yw + tsai_cc.r9 * zw + tsai_cc.Tz;
}


/***********************************************************************\
* This routine takes the position of a point in camera coordinates [mm]	*
* and determines its position in world coordinates [mm].		*
\***********************************************************************/
void      camera_coord_to_world_coord (xc, yc, zc, xw, yw, zw)
    double    xc,
              yc,
              zc,
             *xw,
             *yw,
	     *zw;
{
    double    common_denominator;

    /* these equations were found by simply inverting the previous routine using Macsyma */

    common_denominator = ((tsai_cc.r1 * tsai_cc.r5 - tsai_cc.r2 * tsai_cc.r4) * tsai_cc.r9 +
			  (tsai_cc.r3 * tsai_cc.r4 - tsai_cc.r1 * tsai_cc.r6) * tsai_cc.r8 +
			  (tsai_cc.r2 * tsai_cc.r6 - tsai_cc.r3 * tsai_cc.r5) * tsai_cc.r7);

    *xw = ((tsai_cc.r2 * tsai_cc.r6 - tsai_cc.r3 * tsai_cc.r5) * zc +
	   (tsai_cc.r3 * tsai_cc.r8 - tsai_cc.r2 * tsai_cc.r9) * yc +
	   (tsai_cc.r5 * tsai_cc.r9 - tsai_cc.r6 * tsai_cc.r8) * xc +
	   (tsai_cc.r3 * tsai_cc.r5 - tsai_cc.r2 * tsai_cc.r6) * tsai_cc.Tz +
	   (tsai_cc.r2 * tsai_cc.r9 - tsai_cc.r3 * tsai_cc.r8) * tsai_cc.Ty +
	   (tsai_cc.r6 * tsai_cc.r8 - tsai_cc.r5 * tsai_cc.r9) * tsai_cc.Tx) / common_denominator;

    *yw = -((tsai_cc.r1 * tsai_cc.r6 - tsai_cc.r3 * tsai_cc.r4) * zc +
	    (tsai_cc.r3 * tsai_cc.r7 - tsai_cc.r1 * tsai_cc.r9) * yc +
	    (tsai_cc.r4 * tsai_cc.r9 - tsai_cc.r6 * tsai_cc.r7) * xc +
	    (tsai_cc.r3 * tsai_cc.r4 - tsai_cc.r1 * tsai_cc.r6) * tsai_cc.Tz +
	    (tsai_cc.r1 * tsai_cc.r9 - tsai_cc.r3 * tsai_cc.r7) * tsai_cc.Ty +
	    (tsai_cc.r6 * tsai_cc.r7 - tsai_cc.r4 * tsai_cc.r9) * tsai_cc.Tx) / common_denominator;

    *zw = ((tsai_cc.r1 * tsai_cc.r5 - tsai_cc.r2 * tsai_cc.r4) * zc +
	   (tsai_cc.r2 * tsai_cc.r7 - tsai_cc.r1 * tsai_cc.r8) * yc +
	   (tsai_cc.r4 * tsai_cc.r8 - tsai_cc.r5 * tsai_cc.r7) * xc +
	   (tsai_cc.r2 * tsai_cc.r4 - tsai_cc.r1 * tsai_cc.r5) * tsai_cc.Tz +
	   (tsai_cc.r1 * tsai_cc.r8 - tsai_cc.r2 * tsai_cc.r7) * tsai_cc.Ty +
	   (tsai_cc.r5 * tsai_cc.r7 - tsai_cc.r4 * tsai_cc.r8) * tsai_cc.Tx) / common_denominator;
}
