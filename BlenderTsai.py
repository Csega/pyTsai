#!BPY
#
# BlenderTsai.py
# Copyright (C) Jonathan Merritt 2004.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
#



################
# IMPORT BLOCK #
################

import Blender, math
from Blender import BGL, Draw, Image, Object, Scene, Text, Window, Registry

def menuError(str):
	"""
	Displays an error string in a pop-up menu.
	"""
	Draw.PupMenu('ERROR%%t|%s' % str)

try:
	import Tsai
except ImportError:
	print 'This script requires the Tsai extension'
	print 'module to be installed.  This module could'
	print 'not be imported.  Please make sure that it'
	print 'is installed.'
	menuError('Could not import Tsai extension module.')



####################
# GLOBAL CONSTANTS #
####################

BUTTON_QUIT      =  1
BUTTON_LOAD      =  2
BUTTON_ZOOM      =  3
BUTTON_ADD       =  4
BUTTON_DELETE    =  5
BUTTON_FULLOPT   =  6
BUTTON_COPLANAR  =  7
BUTTON_CALIBRATE =  8
BUTTON_OFSX      =  9
BUTTON_OFSY      = 10
BUTTON_OFSZ      = 11
ZOOM_MIN         = 0.1
ZOOM_MAX         = 10.0
ZOOM_DELTA       = 0.075
OFS_MIN          = -100.0
OFS_MAX          = 100.0
MODE_NORMAL      = 1
MODE_ADD         = 2
MODE_GRAB        = 3
MAX_SELECT_DIST  = 30
SCRIPTLINK_NAME  = '__TsaiCC_REDRAW'

# These colors should later be set by the Theme module
# once it is implemented.
COLOR_BACKGROUND = (130/255.0, 130/255.0, 130/255.0, 1.0)
COLOR_CROSSHAIRS = (1.0, 0.5, 0.5, 0.4)
COLOR_MINIMAP    = (0.0, 0.0, 0.0, 1.0)
COLOR_VERTSEL    = (255/255.0, 255/255.0, 112/255.0, 1.0)
COLOR_VERTUNSEL  = (255/255.0, 112/255.0, 255/255.0, 1.0)
COLOR_TEXT       = (1.0, 1.0, 1.0, 1.0)
POINT_SIZE       = 4


####################
# GLOBAL VARIABLES #
####################

class ADict:
	pass
G = ADict()
G.buttons = ADict()
G.buttons.quit      = None
G.buttons.load      = None
G.buttons.zoom      = None
G.buttons.fullopt   = None
G.buttons.coplanar  = None
G.buttons.ofsx      = None
G.buttons.ofsy      = None
G.buttons.ofsz      = None
G.buttons.calibrate = None
G.buttons.add       = None
G.buttons.delete    = None
G.zoom        = 1.0
G.fullopt     = 0
G.coplanar    = 1
G.ofsx        = 0.0
G.ofsy        = 0.0
G.ofsz        = 0.0
G.imagename   = None
G.image       = None
#G.ibuf        = None
G.iw          = 0
G.ih          = 0
G.imgpos      = [10,10]
G.mousepos    = None
G.grabpos     = None
G.mode        = MODE_NORMAL
G.selection   = []
G.last_camera = None
G.curempty    = None
# coordmap is a map between the names of Blender
#   empties and 2D image coordinates (stored as
#   tuples).
G.coordmap = {}
	


###########
# METHODS #
###########


def addCustomScriptLink():

	"""
	Adds a custom script link to the scene so that all script
	windows are redrawn whenever the scene is redrawn.
	"""

	txt = Text.New(SCRIPTLINK_NAME)
	txt.write('import Blender\n')
	txt.write('from Blender import Window\n')
	txt.write('Window.Redraw(Window.Types.SCRIPT)\n')
	Scene.GetCurrent().addScriptLink(SCRIPTLINK_NAME, 'Redraw')

	# do a full scene update
	Scene.GetCurrent().update(1)


def removeCustomScriptLink():

	"""
	Removes the custom script link (if it is present) from the
	scene.
	"""
	
	try:
		txt = Text.Get(SCRIPTLINK_NAME)
	except NameError:    # script link not found
		return
	
	scene = Scene.GetCurrent()
	slinks = scene.getScriptLinks('Redraw')
	if slinks is not None:
		scene.clearScriptLinks()
		if SCRIPTLINK_NAME in slinks:
			slinks.remove(SCRIPTLINK_NAME)
		for link in slinks:
			scene.addScriptLink(link, 'Redraw')
	
	Text.unlink(txt)

	# do a full scene update
	Scene.GetCurrent().update(1)


def getWinRect():

	"""
	Returns the rectangle of the current script window in
	screen coordinates.
	
	@return: Script window rectangle: [x, y, w, h]
	"""
	
	winrect = BGL.Buffer(BGL.GL_FLOAT, 4)
	BGL.glGetFloatv(BGL.GL_SCISSOR_BOX, winrect)
	return winrect.list
		

def getWMCoords():

	"""
	Returns the coordinates of the mouse relative to the
	current script window: (x,y)
	
	@return: Current script window mouse coordinates: (x,y).
	"""
	
	xm,ym = Window.GetMouseCoords()
	xw,yw,ww,hw = getWinRect()
	return (xm-xw, ym-yw)


def wc2ic(p):

	"""
	Converts window coordinates to image coordinates.
	
	@param p: Window coordinates (x,y).
	
	@return: Image coordinates (x,y).
	"""
	
	xi = p[0] / G.zoom - int(G.imgpos[0] / G.zoom)
	yi = p[1] / G.zoom - int(G.imgpos[1] / G.zoom)
	if xi < 0: xi = 0
	if yi < 0: yi = 0
	if G.image is not None:
		if xi > G.iw: xi = G.iw
		if yi > G.ih: yi = G.ih
	return (xi, yi)


def setZoom(p, zoom):

	"""
	Sets the zoom factor, G.zoom.  This method maintains the
	centre of zoom (p) at the same location before and after
	the zoom.
	
	@param p:    Centre of the zoom (xw,yw), in window coordinates.
	@param zoom: New zoom factor.
	"""
	
	# bounds-check the zoom factor
	if zoom < ZOOM_MIN: zoom = ZOOM_MIN
	if zoom > ZOOM_MAX: zoom = ZOOM_MAX
	
	# if we have no image then just set the zoom
	global G
	if G.image is None:
		G.zoom = zoom
		Draw.Redraw(1)
		return
	
	# find centre of zoom in image coordinates
	xi = (p[0] - G.imgpos[0]) / G.zoom
	yi = (p[1] - G.imgpos[1]) / G.zoom
	
	# set the new zoom
	G.zoom = zoom

	# find the new centre of zoom
	xinew = (p[0] - G.imgpos[0]) / G.zoom
	yinew = (p[1] - G.imgpos[1]) / G.zoom
	
	# find the offset for the origin
	xofs, yofs = int(xinew - xi), int(yinew - yi)
	
	# change the origin
	G.imgpos[0] += xofs * G.zoom
	G.imgpos[1] += yofs * G.zoom
	
	# queue a redraw
	Draw.Redraw(1)


def calibrate(camera = None):

	"""
	Performs calibration.
	"""
	
	# TODO: Add warnings about insufficient numbers of points

	if camera == None:
		camera = G.selection[0]
	
	if G.coplanar: target_type = 'coplanar'
	else:          target_type = 'noncoplanar'
	
	if G.fullopt:  optimization_type = 'full'
	else:          optimization_type = 'three-param'
	
	origin_offset = (G.ofsx, G.ofsy, G.ofsz)
	
	calcoords = []
	for (emptyname, imgcoords) in G.coordmap.items():
		m = Object.Get(emptyname).getMatrix('worldspace')
		xs = m[3][0]
		ys = m[3][1]
		zs = m[3][2]
		xi = imgcoords[0]
		yi = G.ih - imgcoords[1]    # note the change to the y-axis there!
		#print xs, ys, zs, xi, yi
		calcoords.append((xs, ys, zs, xi, yi))
	
	cp = Tsai.CameraParameters(image_dim=(G.iw, G.ih))
	try:
		cp = Tsai.calibrate(target_type, optimization_type, calcoords, cp, origin_offset)
	except Tsai.CalibrationError, msg:
		menuError(msg)
		return
	
	cp.setBlenderCamera(camera, G.iw, G.ih)
	
	Window.RedrawAll()
	

def loadTGARAW(fileName):
	
	"""
	Loads a RAW TGA image from a file.  On success, the method
	returns: (image, imagebuffer).  On failure, the method
	returns None.
	
	@param filename: Name of the RAW TGA file to load.
	
	@return: (image, imagebuffer) 2-tuple on success, or None
	         on failure.
	"""
	
	# create the image and allocate the image buffer using
	#  Blender's own methods
	image = Image.Load(fileName)
	#ibuf = BGL.Buffer(BGL.GL_BYTE, (image.size[1], image.size[0], 4))
	
	## read in the file
	#infile = open(fileName, 'rb')
	#array = map(ord, infile.read())
	#infile.close()
	
	## check the validity of the file and offset the data
	##  array to the start of the image data
	#if array[2] is not 2:
	#	menuError('Can only load RAW TGA images.')
	#	return None
	#array = array[(18+array[0]):]

	# read the array into the image buffer
	#index = 0
	#for y in range(0, image.size[1]):
	#	for x in range(0, image.size[0]):
	#		ibuf[y][x][2] = array[index]
	#		ibuf[y][x][1] = array[index+1]
	#		ibuf[y][x][0] = array[index+2]
	#		index += 3
	
	# return the (image, imagebuffer) tuple
	return (image,)


def loadImage(fileName):

	"""
	Handles loading an image.  This method populates
	G.image, G.ibuf, G.iw, G.ih.
	
	@param fileName: Name of the file to load.
	"""
	
	Window.WaitCursor(True)
	
	# load the image and set up parameters
	global G
	try:
		iloaded = loadTGARAW(fileName)
	except Exception:
		menuError('Could not load image.')
		iloaded = None
	if iloaded is None:
		return
	(G.image,) = iloaded
	G.iw, G.ih = G.image.size
	G.imagename = fileName
	
	# centre the image in the window rectangle
	xw,yw,ww,hw = getWinRect()
	G.imgpos[0] = (ww-G.iw) / 2
	G.imgpos[1] = (hw-G.ih) / 2
	
	Window.WaitCursor(False)

	



def removeUnknownsFromCoords():

	"""
	Checks that all names of empties in G.coordmap
	actually exist in the scene.  If they don't, they're
	pruned from the map.
	"""
	
	global G
	
	for emptyname in G.coordmap.keys():
		try:
			Object.Get(emptyname)
		except AttributeError:
			del G.coordmap[emptyname]


def writeDataToRegistry():

	"""
	Writes data to the Blender registry.
	"""

	global G

	dict = {}	
	dict['coordmap'] = G.coordmap
	if G.imagename is not None:
		dict['imagename'] = G.imagename
	
	Registry.SetKey('TsaiCC', dict)
	
	
def readDataFromRegistry():

	"""
	Reads data from the Blender registry.
	"""

	global G
	
	dict = Registry.GetKey('TsaiCC')
	if dict:
		G.coordmap = dict['coordmap']
		removeUnknownsFromCoords()
		# TODO:
		try:
			loadImage(dict['imagename'])
		except KeyError:
			pass
		

def renderGUI():

	"""
	Renders the GUI for the script.
	"""
	
	global G

	# find the selection set and update some selection
	#  related flags
	haveEmpty  = False
	haveCamera = False
	emptyName  = ""
	G.selection = Object.GetSelected()
	G.curempty  = None
	if G.selection is not None and len(G.selection) > 0:
		mso = G.selection[0]
		msotype = mso.getType()
		if msotype == 'Empty':
			haveEmpty = True
			G.curempty = mso
			emptyName = G.curempty.getName()
		elif msotype == 'Camera':
			haveCamera = True
	emptyHasCoords = G.coordmap.has_key(emptyName)
	removeUnknownsFromCoords()
	
	# clear any buttons that need to have set states
	G.buttons.add = G.buttons.delete = None
	
	# clear the window
	c = COLOR_BACKGROUND
	BGL.glClearColor(c[0], c[1], c[2], c[3])
	BGL.glClear(BGL.GL_COLOR_BUFFER_BIT)

	# paint the image in the background
	if G.image is not None:
		#drawImage(G.image, G.imgpos, G.iw, G.ih, G.zoom)
		Draw.Image(G.image, G.imgpos[0], G.imgpos[1], G.zoom, G.zoom)
	# paint 2D vertices in the image
	BGL.glPushAttrib(BGL.GL_POINT_BIT | BGL.GL_CURRENT_BIT)
	BGL.glPointSize(POINT_SIZE)
	x0 = int(G.imgpos[0] / G.zoom)
	y0 = int(G.imgpos[1] / G.zoom)
	def drawvc(ec):
		emptyname, coord = ec
		if Object.Get(emptyname) in G.selection:
			c = COLOR_VERTSEL
		else:
			c = COLOR_VERTUNSEL
		BGL.glColor4f(c[0], c[1], c[2], c[3])
		BGL.glVertex2f(G.zoom*(coord[0]+x0), G.zoom*(coord[1]+y0))
	BGL.glBegin(BGL.GL_POINTS)
	map(drawvc, G.coordmap.items())
	BGL.glEnd()
	BGL.glPopAttrib()
	
	# if we're in add mode then draw some extra stuff
	if G.mode == MODE_ADD:
	
		xm,ym = map(int, getWMCoords())
		xm -= 10
		ym += 10
		(xw,yw,ww,hw) = getWinRect()
		
		# draw crosshairs
		c = COLOR_CROSSHAIRS
		BGL.glColor4f(c[0], c[1], c[2], c[3])		
		verts = [ (xm,0), (xm,hw), (0,ym), (ww,ym) ]
		BGL.glBegin(BGL.GL_LINES)
		map(lambda x: BGL.glVertex2d(x[0],x[1]), verts)
		BGL.glEnd()

		#############################################
		# UNCOMMENT THIS SECTION FOR A COOL MINIMAP
		# EFFECT - NOT VERY USEFUL THOUGH
		#
		## draw "minimap" background
		#c = COLOR_MINIMAP
		#BGL.glColor4f(c[0], c[1], c[2], c[3])
		#verts = [ (119,10), (221,10), (221,111), (119,111) ]
		#BGL.glBegin(BGL.GL_QUADS)
		#map(lambda x: BGL.glVertex2i(x[0],x[1]), verts)
		#BGL.glEnd()
		#
		## paint the image into the minimap
		#ix,iy = wc2ic((xm,ym))
		#ix,iy = map(int, [ix, iy])
		#drawImage(G.ibuf, (120,10), G.iw, G.ih, 10.0, (ix-5,iy-5,10,10))
		#
		# END OF MINIMAP SECTION
		#############################################
	
	# paint the current empty name
	if haveEmpty:
		c = COLOR_TEXT
		BGL.glColor4d(c[0], c[1], c[2], c[3])
		BGL.glRasterPos2i(220, 10)
		Draw.Text(emptyName, 'small')
	
	# paint the normal GUI buttons
	G.buttons.quit = Draw.PushButton('Quit', BUTTON_QUIT, 5, 5, 100, 20, 'Exits the script.')
	G.buttons.load = Draw.PushButton('Load Image', BUTTON_LOAD, 5, 25, 100, 20, 'Loads an image.')
	G.buttons.zoom = Draw.Number('Zoom', BUTTON_ZOOM, 5, 45, 100, 20, G.zoom, ZOOM_MIN, ZOOM_MAX, 'Adjusts image zoom.')

	# paint camera-specific stuff
	if haveCamera:
		G.buttons.fullopt   = Draw.Toggle('Full Optimization', BUTTON_FULLOPT, 110, 5, 120, 20, G.fullopt, 'Full or partial optimization.')
		G.buttons.coplanar  = Draw.Toggle('Coplanar', BUTTON_COPLANAR, 110, 25, 120, 20, G.coplanar, 'Coplanar or non-coplanar target.')
		# Origin offset is not currently working in the Tsai module.
		#  It should be brought back here when it is.
		#G.buttons.ofsz      = Draw.Number('OfsZ', BUTTON_OFSZ, 110, 50, 100, 20, G.ofsz, OFS_MIN, OFS_MAX, 'Z origin offset.')
		#G.buttons.ofsy      = Draw.Number('OfsY', BUTTON_OFSY, 110, 70, 100, 20, G.ofsy, OFS_MIN, OFS_MAX, 'Y origin offset.')
		#G.buttons.ofsx      = Draw.Number('OfsX', BUTTON_OFSX, 110, 90, 100, 20, G.ofsx, OFS_MIN, OFS_MAX, 'X origin offset.')
		G.buttons.calibrate = Draw.PushButton('Calibrate', BUTTON_CALIBRATE, 235, 5, 100, 20, 'Calibrates the selected camera.')
	
	# paint empty-specific stuff
	elif haveEmpty and (G.mode == MODE_NORMAL):
		if emptyHasCoords:
			G.buttons.delete = Draw.PushButton('Delete', BUTTON_DELETE, 110, 5, 100, 20, 'Adds an image calibration coordinate.')
		else:
			G.buttons.add = Draw.PushButton('Add', BUTTON_ADD, 110, 5, 100, 20, 'Removes an image calibration coordinate.')


def eventHandler(event, value):

	"""
	General GUI event handler.

	@param event: Event type.
	@param value: Value of the event.
	"""
	
	global G
	
	if event == Draw.WHEELDOWNMOUSE:
		setZoom(getWMCoords(), G.zoom*(1.0-ZOOM_DELTA))

	elif event == Draw.WHEELUPMOUSE:
		setZoom(getWMCoords(), G.zoom*(1.0+ZOOM_DELTA))
	
	elif event == Draw.MIDDLEMOUSE:
		if value: G.mousepos = Window.GetMouseCoords()
		else:     G.mousepos = None

	elif event == Draw.MOUSEX or event == Draw.MOUSEY:
	
		mouseButs = Window.GetMouseButtons()
		
		if (mouseButs & Draw.MIDDLEMOUSE) and (G.mousepos is not None):
			nx,ny = Window.GetMouseCoords()
			dx,dy = nx-G.mousepos[0], ny-G.mousepos[1]
			G.mousepos = (nx,ny)
			G.imgpos = [int(G.imgpos[0]+dx), int(G.imgpos[1]+dy)]
			Draw.Redraw(1)

		elif G.mode == MODE_ADD:
			Draw.Redraw(1)
			
		elif G.mode == MODE_GRAB:
			G.havebupclik = True
			nx,ny = Window.GetMouseCoords()
			dx,dy = (nx-G.grabpos[0])/G.zoom, (ny-G.grabpos[1])/G.zoom
			G.grabpos = [nx,ny]
			def translate(x):
				name = x.getName()
				if G.coordmap.has_key(name):
					c = G.coordmap[name]
					G.coordmap[name] = (c[0]+dx, c[1]+dy)
			map(translate, G.selection)
			## autocalibration.. gets stuck some times..
			#if G.last_camera:
			#	calibrate(G.last_camera)
			Draw.Redraw(1)

	elif (event == Draw.LEFTMOUSE) or (event == Draw.RETKEY):
	
		if (G.mode == MODE_ADD) and (value == 1):
			x,y = map(int, getWMCoords())
			x -= 10
			y += 10
			G.coordmap[G.curempty.getName()] = wc2ic((x,y))
			G.mode = MODE_NORMAL
			Draw.Redraw(1)
		
		elif (G.mode == MODE_GRAB) and (value == 1):
			G.mode = MODE_NORMAL
			Draw.Redraw(1)
	
	elif (event == Draw.RIGHTMOUSE) and (value == 1):
	
		if G.mode == MODE_NORMAL:
			xi,yi = wc2ic(getWMCoords())
			closest = None
			for (emptyname, coord) in G.coordmap.items():
				dist = math.sqrt((coord[0]-xi)**2 + (coord[1]-yi)**2) * G.zoom
				if (closest == None) or (dist < closest[0]):
					closest = (dist, emptyname)
			if closest[0] < MAX_SELECT_DIST:
				obj = Object.Get(closest[1])
				kq = Window.GetKeyQualifiers()
				if (kq & Window.Qual.LSHIFT) or (kq & Window.Qual.RSHIFT):
					obj.select(True)
				else:
					map(lambda x: x.select(False), G.selection)
					obj.select(True)
				Window.RedrawAll()
	
	elif (event == Draw.AKEY) and (value == 1):
	
		if G.mode == MODE_NORMAL:
			someSelected = False
			for (emptyname, coord) in G.coordmap.items():
				if Object.Get(emptyname).isSelected():
					someSelected = True
					break
			newselect = (someSelected == False)
			map(lambda x: Object.Get(x[0]).select(newselect), G.coordmap.items())
			Window.RedrawAll()
	
	elif (event == Draw.GKEY) and (value == 1):

		if G.mode == MODE_NORMAL:
			G.mode = MODE_GRAB
			G.grabpos = Window.GetMouseCoords()
			Draw.Redraw(1)
	
	elif event == Draw.UPARROWKEY and value == 1:
	
		p = Window.GetMouseCoords()
		Window.SetMouseCoords(p[0], p[1]+1)
		
	elif event == Draw.DOWNARROWKEY and value == 1:

		p = Window.GetMouseCoords()
		Window.SetMouseCoords(p[0], p[1]-1)
		
	elif event == Draw.LEFTARROWKEY and value == 1:
	
		p = Window.GetMouseCoords()
		Window.SetMouseCoords(p[0]-1, p[1])
	
	elif event == Draw.RIGHTARROWKEY and value == 1:
	
		p = Window.GetMouseCoords()
		Window.SetMouseCoords(p[0]+1, p[1])
		
	elif event == Draw.XKEY and value == 1:
	
		if len(G.selection) > 0:
			result = Draw.PupMenu('OK?%t|Erase selected')
			if result == 1:
				buttonEventHandler(BUTTON_DELETE)
	
	elif event == Draw.SPACEKEY and value == 1:
	
		if (G.curempty is not None) and not (G.coordmap.has_key(G.curempty.getName())):
			buttonEventHandler(BUTTON_ADD)
		
	elif event == Draw.ZKEY and value == 1:
	
		x,y,w,h = getWinRect()
		setZoom((w/2, h/2), 1.0)

	elif event == Draw.RKEY and value == 1:
	
		Draw.Redraw(1)


def buttonEventHandler(button):

	"""
	Event handler for button presses.
	
	@param button: Button ID.
	"""

	global G
	G.havebupclik = False
	
	if button == BUTTON_QUIT:
		removeCustomScriptLink()
		writeDataToRegistry()
		Draw.Exit()
		Window.RedrawAll()
	
	elif button == BUTTON_LOAD:
		G.havebupclik = True
		Window.ImageSelector(loadImage)

	elif button == BUTTON_ZOOM:
		x,y,w,h = getWinRect()
		setZoom((w/2, h/2), G.buttons.zoom.val)
	
	elif button == BUTTON_FULLOPT:
		G.fullopt = G.buttons.fullopt.val
		
	elif button == BUTTON_COPLANAR:
		G.coplanar = G.buttons.coplanar.val
		
	elif button == BUTTON_OFSX:
		G.ofsx = G.buttons.ofsx.val
		
	elif button == BUTTON_OFSY:
		G.ofsy = G.buttons.ofsy.val
		
	elif button == BUTTON_OFSZ:
		G.ofsz = G.buttons.ofsz.val
	
	elif button == BUTTON_ADD:
		G.mode = MODE_ADD
		Draw.Redraw(1)
	
	elif button == BUTTON_DELETE:
		def delmap(x): 
			del G.coordmap[x.getName()]
			x.select(False)
		map(delmap, G.selection)
		Window.RedrawAll()
	
	elif button == BUTTON_CALIBRATE:
		Window.WaitCursor(True)
		G.last_camera = G.selection[0]
		calibrate()


###############
# ENTRY POINT #
###############

# remove any spurious custom script link that may be
#  hanging around from a previous invokation
removeCustomScriptLink()

# add a new script link
addCustomScriptLink()

# read any previous coordinate data that has been
#  stored
readDataFromRegistry()

# register event handlers
Draw.Register(renderGUI, eventHandler, buttonEventHandler)