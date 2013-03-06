"""
Routines for calibrating a camera using the method of Tsai.
"""

import math, string

import pytsai


class CalibrationError(Exception):
        """
        Exception used to indicate that an error has occurred during
        calibration.
        """
        def __init__(self, value):
                self.value = value
        def __str__(self):
                return str(self.value)


class CameraParameters:

        """
        Utility class for camera parameters.
        
        This class can be used as an input and output to the L{calibrate}
        method.  It functions as a mapping type between camera parameters
        (stored as strings), and their values (numbers).
        """
        
        def __init__(self, existing=None, **keywords):
        
                """
                @param existing: A mapping type to copy values from.
                @keyword model: Specifies an existing (known) camera model.  
                                This can be any of:
                                        - 'photometrics-star-I'
                                        - 'general-imaging-mos5400-matrox'
                                        - 'panasonic-GP-MF702-matrox'
                                        - 'sony-xc75-matrox'
                                        - 'sony-xc77-matrox'
                                        - 'sony-xc57-androx'
                                        - 'xapshot-matrox'
                @keyword image_dim: Specifies dimensions of an image to derive
                                a camera model from.  The image_dim should be
                                a tuple containing M{(width, height)} of the
                                image.  The derived camera model presumes:
                                        - M{Ncx = Nfx = width}
                                        - M{dx = dy = dpx = dpy = 1.0}
                                        - M{Cx = width / 2}
                                        - M{Cy = height / 2}
                                        - M{sx = 1.0}
                """
                
                # initially create all parameters, setting them to 0.0
                parms = ['Ncx', 'Nfx', 'dx', 'dy', 'dpx', 'dpy', 'Cx', 'Cy',
                         'sx', 'f', 'kappa1', 'p1', 'p2', 'Tx', 'Ty', 'Tz',
                         'Rx', 'Ry', 'Rz', 'r1', 'r2', 'r3', 'r4', 'r5',
                         'r6', 'r7', 'r8', 'r9']
                def createParam(x): self[x]=0.0
                list(map(createParam, parms)) #map() function returns iterator in python 3. This is a workaround.

                # copy an existing map if one was provided
                if existing is not None:
                        for (key,item) in existing.items():
                                self[key] = item
                
                # set the camera model if specified
                if 'model' in keywords.keys():
                        model = keywords['model']
                        if model == 'photometrics-star-I':
                                self.Ncx = 576
                                self.Nfx = 576
                                self.dx  = 0.023
                                self.dy  = 0.023
                                self.dpx = self.dx * self.Ncx / self.Nfx
                                self.dpy = self.dy
                                self.Cx  = 258.0
                                self.Cy  = 204.0
                                self.sx  = 1.0
                        elif model == 'general-imaging-mos5400-matrox':
                                self.Ncx = 649
                                self.Nfx = 512
                                self.dx  = 0.015
                                self.dy  = 0.015
                                self.dpx = self.dx * self.Ncx / self.Nfx
                                self.dpy = self.dy
                                self.Cx  = 512/2
                                self.Cy  = 480/2
                                self.sx  = 1.0
                        elif model == 'panasonic-GP-MF702-matrox':
                                self.Ncx = 649
                                self.Nfx = 512
                                self.dx  = 0.015
                                self.dy  = 0.015
                                self.dpx = self.dx * self.Ncx / self.Nfx
                                self.dpy = self.dy
                                self.Cx  = 268
                                self.Cy  = 248
                                self.sx  = 1.078647
                        elif model == 'sony-xc75-matrox':
                                self.Ncx = 768
                                self.Nfx = 512
                                self.dx  = 0.0084
                                self.dy  = 0.0098
                                self.dpx = self.dx * self.Ncx / self.Nfx
                                self.dpy = self.dy
                                self.Cx  = 512 / 2
                                self.Cy  = 480 / 2
                                self.sx  = 1.0
                        elif model == 'sony-xc77-matrox':
                                self.Ncx = 768
                                self.Nfx = 512
                                self.dx  = 0.011
                                self.dy  = 0.013
                                self.dpx = self.dx * self.Ncx / self.Nfx
                                self.dpy = self.dy
                                self.Cx  = 512 / 2
                                self.Cy  = 480 / 2
                                self.sx  = 1.0
                        elif model == 'sony-xc57-androx':
                                self.Ncx = 510
                                self.Nfx = 512
                                self.dx  = 0.017
                                self.dy  = 0.013
                                self.dpx = self.dx * self.Ncx / self.Nfx
                                self.dpy = self.dy
                                self.Cx  = 512 / 2
                                self.Cy  = 480 / 2
                                self.sx  = 1.107914
                        elif model == 'xapshot-matrox':
                                self.Ncx = 739
                                self.Nfx = 512
                                self.dx  = 6.4 / 782.0
                                self.dy  = 4.8 / 250.0
                                self.dpx = self.dx * self.Ncx / self.Nfx
                                self.dpy = self.dy
                                self.Cx  = 512 / 2
                                self.Cy  = 240 / 2
                                self.sx  = 1.027753
                                
                # construct an artificial camera if image_dim has been
                #  provided
                elif 'image_dim' in keywords.keys():
                        image_dim = keywords['image_dim']
                        self.Ncx = image_dim[0]
                        self.Nfx = self.Ncx
                        self.dx  = 1.0
                        self.dy  = self.dx
                        self.dpx = self.dx
                        self.dpy = self.dx
                        self.Cx  = image_dim[0] / 2.0
                        self.Cy  = image_dim[1] / 2.0
                        self.sx  = 1.0
        
        def getPosition(self):
                """
                Returns the camera's position as a 3-tuple: (Tx, Ty, Tz).
                """
                return (self.Tx, self.Ty, self.Tz)
        
        def getEulerRotation(self):
                """
                Returns the camera's Euler rotation as a 3-tuple:
                (Rx, Ry, Rz).
                """
                return (self.Rx, self.Ry, self.Rz)

        def getRotationMatrix(self):
                """
                Returns the camera's rotation matrix::
                        [
                                [ r1, r2, r3 ],
                                [ r4, r5, r6 ],
                                [ r7, r8, r9 ]
                        ]
                """
                return [ [ self.r1, self.r2, self.r3 ], 
                         [ self.r4, self.r5, self.r6 ], 
                         [ self.r7, self.r8, self.r9 ] ]

        def getFOVx(self):
                """
                Returns the camera's x field-of-view angle (in radians).

                This angle is defined as: M{fovx = 2 * atan2(Ncx * dx, 2*f)}.
                """
                return 2.0 * math.atan2(self.Ncx * self.dx, 2.0 * self.f);

        def world2image(self, coord):
                """
                Converts from world coordinates to image coordinates.  The
                conversion is based upon the camera's current parameters.

                @param coord: A 3-member sequence containing the world
                        coordinates M{(xw, yw, zw)}.

                @return: A 2-member sequence containing the image coordinates
                        M{(xi, yi)}.
                """
                return pytsai._pytsai_wc2ic(coord, self)

        def image2world(self, coord):
                """
                Converts from image coordinates to world coordinates.  The
                conversion is based upon the camera's current parameters.

                @param coord: A 3-member sequence containing the image
                        coordinates and the depth coordinate in the world
                        M{(xi, yi, zw)}.

                @return: A 3-member sequence containing the world coordinates
                        M{(xw, yw, zw)}.
                """
                return pytsai._pytsai_ic2wc(coord, self)

        def world2camera(self, coord):
                """
                Converts from world coordinates to camera coordinates.  The
                conversion is based upon the camera's current parameters.

                @param coord: A 3-member sequence containing the world
                        coordinates M{(xw, yw, zw)}.

                @return: A 3-member sequence containing the camera coordinates
                        M{(xc, yc, zc)}.
                """
                xw, yw, zw = coord
                xc = self.r1 * xw + self.r2 * yw + self.r3 * zw + self.Tx
                yc = self.r4 * xw + self.r5 * yw + self.r6 * zw + self.Ty
                zc = self.r7 * xw + self.r8 * yw + self.r9 * zw + self.Tz
                return (xc, yc, zc)

        def camera2world(self, coord):
                """
                Converts from camera coordinates to world coordinates.  The
                conversion is based upon the camera's current parameters.

                @param coord: A 3-member sequence containing the camera
                        coordinates M{(xc, yc, zc)}.

                @return: A 3-member sequence containing the world
                        coordinates M{(xw, yw, zw)}.
                """
                return pytsai._pytsai_cc2wc(coord, self)

        def removeRadialDistortion(self, coord, type='sensor'):
                """
                Removes distortion from either image or sensor coordinates.
                The conversion is based upon the camera's current
                parameters.

                @param coord: A 2-member sequence containing the distorted
                        image or sensor coordinates.
                @param type: The type of the coordinates ('sensor' or
                        'image').
                @return: Un-distorted sensor or image coordinates (x,y).
                """
                if type == 'sensor':
                        Xd, Yd = coord
                        dfac = 1.0 + self.kappa1 * (Xd*Xd + Yd*Yd)
                        return (Xd*dfac, Yd*dfac)
                elif type == 'image':
                        Xfd,Yfd = coord
                        Xd = self.dpx * (Xfd - self.Cx) / self.sx
                        Yd = self.dpy * (Yfd - self.Cy)
                        Xu,Yu = removeDistortion((Xd,Yd), 'sensor')
                        Xfu = Xu*self.sx/self.dpx + self.Cx
                        Yfu = Yu/self.dpy + self.Cy
                        return (Xfu, Yfu)
                else:
                        raise TypeError('Unknown coordinate type: %s' % \
                                type)

        def addRadialDistortion(self, coord, type='sensor'):
                """
                Adds distortion to either image or sensor coordinates.
                The conversion is based upon the camera's current
                parameters.

                @param coord: A 2-member sequence containing the un-distorted
                        image or sensor coordinates.
                @param type: The type of the coordinates ('sensor' or
                        'image').
                @return: Distorted sensor or image coordinates (x,y).
                """
                if type == 'sensor':
                        Xu,Yu = coord
                        return pytsai._pytsai_add_sensor_coord_distortion(
                                Xu, Yu, self)
                elif type == 'image':
                        Xfu,Yfu = coord
                        Xu = self.dpx * (Xfu - self.Cx) / self.sx
                        Yu = self.dpy * (Yfu - self.Cy)
                        Xd,Yd = addDistortion((Xu, Yu), 'sensor')
                        Xfd = Xd*self.sx/self.dpx + self.Cx
                        Yfd = Yd/self.dpy + self.Cy
                        return (Xfd,Yfd)
                else:
                        raise TypeError('Unknown coordinate type: %s' % \
                                type)

        def setBlenderCamera(self, camobj, xres, yres):
                """
                This function is now fully compatible with Blender 2.65.
                (Haven't tested yet.)
                """
                import mathutils
                w2c = mathutils.Matrix((
                        (self.r1, self.r4, self.r7, 0.0),
                        (self.r2, self.r5, self.r8, 0.0),
                        (self.r3, self.r6, self.r9, 0.0),
                        (self.Tx, self.Ty, self.Tz, 1.0)
                ))
                rot180x = mathutils.Matrix.Rotation(180, 4, 'X')
                c2w = mathutils.Matrix.copy(w2c * rot180x)
                c2w.invert()
                camobj.matrix_world = c2w
                cam = camobj.data
                #asp = float(yres) / float(yres) #isn't it a semantic error???
                asp = float(xres) / float(yres) #test with this - maybe its yres/xres
                cam.lens = 16.0 / (asp * math.tan(self.getFOVx() / 2.0))

        def iterkeys(self):
                """
                Iterates over the keys of the mapping.
                """
                return self.__dict__.iterkeys()
        def __getitem__(self, key):
                return self.__dict__[key]
        def __setitem__(self, key, value):
                self.__dict__[key] = value
        def __delitem__(self, key):
                del self.__dict__[key]
        def __iter__(self):
                return self.iterkeys()
        def __contains__(self, item):
                return (item in self.__dict__)
        def __str__(self):
                str = 'CameraParameters:\n'
                parms = ['Ncx', 'Nfx', 'dx', 'dy', 'dpx', 'dpy', 'Cx', 'Cy',
                         'sx', 'f', 'kappa1', 'p1', 'p2', 'Tx', 'Ty', 'Tz',
                         'Rx', 'Ry', 'Rz', 'r1', 'r2', 'r3', 'r4', 'r5',
                         'r6', 'r7', 'r8', 'r9']
                for p in parms:
                        #s = string.ljust(p, 7) + ("= %f\n" % self[p]) #old, python 2 syntax
                        s = p.ljust(7) + ("= %f\n" % self[p])
                        str += s
                return str
        def __repr__(self):
                return str(self)


def calibrate(target_type, optimization_type, calibration_data, camera_params,
        origin_offset=(0.0,0.0,0.0)):
        """
        Calibrates a camera.

        @param target_type: The type of the target used for calibration.  This
                can be either:
                        - 'coplanar', in which all z-values for the calibration
                          points are zero.
                        - 'noncoplanar', in which some z-values for the
                          calibration points must be non-zero.

        @param optimization_type: The type of optimization to perform.  This
                can be either:
                        - 'three-param', for optimization of only M{f}, 
                          M{Tz} and M{kappa1}.
                        - 'full', for full optimization.
                          
        @param calibration_data: A sequence of sequences containing
                calibration points.  The sequence should consist of::
                        [
                                [ xs1, ys1, zs1, xi1, yi1 ],
                                [ xs2, ys2, zs2, xi2, yi2 ],
                                  ...  ...  ...  ...  ...
                                [ xsN, ysN, zsN, xiN, yiN ]
                        ]
                where:
                        - M{(xs, ys, zs)} are 3D space coordinates of the
                          calibration points.
                        - M{(xi, yi)} are corresponding 2D image space 
                          coordinates of the calibration points.
        @param camera_params: A dictionary mapping camera parameter names
                (stored as strings) to their values (which should be
                numbers).  The class L{CameraParameters} is a utility class
                designed to be used in this position.

        @param origin_offset: An artificial offset that is added to the origin
                of the calibration data coordinates.  This offset is later
                removed from the camera position as determined by the
                calibration. Shifting the origin may be useful since the
                Tsai method fails if the world space origin is near to the
                camera space origin or the camera space y axis.
        """

        # add an origin offset to the camera position
        #xo,yo,zo = origin_offset
        xo,yo,zo = 0.0,0.0,0.0
        def addOfs(c):
                return (c[0]+xo, c[1]+yo, c[2]+zo, c[3], c[4])
        ofsCalData = list(map(addOfs, calibration_data))

        # perform camera calibration
        if target_type == 'coplanar' and optimization_type == 'three-param':
                try:
                        cp = pytsai._pytsai_coplanar_calibration(
                                ofsCalData, camera_params)
                except RuntimeError as runtimeError:
                        raise CalibrationError(str(runtimeError))
                        
        elif target_type == 'noncoplanar' and \
             optimization_type == 'three-param':
                try:
                        cp = pytsai._pytsai_noncoplanar_calibration(
                                ofsCalData, camera_params)
                except RuntimeError as runtimeError:
                        raise CalibrationError(str(runtimeError))

        elif target_type == 'coplanar' and optimization_type == 'full':
                try:
                        cp = pytsai._pytsai_coplanar_calibration_fo(
                                ofsCalData, camera_params)
                except RuntimeError as runtimeError:
                        raise CalibrationError(str(runtimeError))

        elif target_type == 'noncoplanar' and optimization_type == 'full':
                try:
                        cp = pytsai._pytsai_noncoplanar_calibration_fo(
                                ofsCalData, camera_params)
                except RuntimeError as runtimeError:
                        raise CalibrationError(str(runtimeError))

        else:
                errstr = 'Unknown combination of target_type=\'%s\' and ' \
                         'optimization_type=\'%s\'' % \
                         (target_type, optimization_type)
                raise CalibrationError(errstr)

        # remove the origin offset from the camera position
        ccp = CameraParameters(cp)
        #camorigin = ccp.world2camera((xo,yo,zo))
        #ccp.Tx -= camorigin[0]
        #ccp.Ty -= camorigin[1]
        #ccp.Tz -= camorigin[2]

        # return the calculated camera parameters
        return ccp
        
