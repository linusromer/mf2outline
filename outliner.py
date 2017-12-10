#!/usr/bin/env python

import math,fontforge


# make a homogeneous transformation m*(x,y,1) or
# truncated matrix m*(x,y)
# and return the transformed x and y
def homogeneous(m,x,y):
	if len(m)==4: # truncated matrix (dtransform)
		return [m[0]*x+m[2]*y,m[1]*x+m[3]*y]
	elif len(m)==6: # general case
		return [m[0]*x+m[2]*y+m[4],m[1]*x+m[3]*y+m[5]]
	else:
		return [x,y] # wrong matrix dimensions, no change

# inverts the transformationmatrix
def invertmatrix(m):
	if len(m)>3:
		c = 1.0/(m[0]*m[3]-m[1]*m[2]) # determinant coefficient
	if len(m)==4: # truncated matrix (dinvert)
		return [c*m[3],-c*m[1],-c*m[2],c*m[0]]
	elif len(m)==6: # general case
		return [c*m[3],-c*m[1],-c*m[2],c*m[0],-m[4],-m[5]]
	else:
		return m # wrong matrix dimensions, no change

# For the own postscript interpreter we use the following 
# convention for cubic bezier paths:
# (0,1)--(3,4)..controls (1,2) and (-2,7)..(8,5)--(2,9)
# maps bijective to
# [[(0,1)],[(3,4)],[(1,2),(-2,7),(8,5)],[2,9]]
# We will NOT store, if the path is cyclic, because for
# elliptical pens this is equivalent to make the first 
# and the last point the same

def vecadd(a,b):
	return (a[0]+b[0],a[1]+b[1])

def vecdiff(a,b):
	return (a[0]-b[0],a[1]-b[1])

def veclen(a):
	return (a[0]**2+a[1]**2)**.5
	
def vecnorm(a):
	return (a[0]/veclen(a),a[1]/veclen(a))
	
	
# scale a vector a such that it has length l
# (if l is negative, it will be in opposite direction)
def vecscaleto(a,l):
	return (a[0]/veclen(a)*l,a[1]/veclen(a)*l)

# directed angle from vector a to b (inbetween -180 and 180)	
def vecangle(a,b):
	return math.degrees(math.atan2(a[0]*b[1]-a[1]*b[0],
	a[0]*b[0]+a[1]*b[1]))

# returns a points right of point p where right
# means rectangular to direction d in distance r
# (choose r negative for a left point)
def pointright(p,d,r):
	dscaled = vecscaleto(d,r)
	return (p[0]+dscaled[1],p[1]-dscaled[0])

# bezierinterpolate returns a bezier path as list (see the convention
# as above) which starts in point za, heading in direction dira, 
# passing zb at time 0.5 and ending in time 1
# at zc heading in direction dirc
def bezierinterpolate(za,dira,zb,zc,dirc):
	# (xp,yp) is the first control point
	# (xq,yq) is the second control point
	# solve([(yp-ya)*xdira=(xp-xa)*ydira,
	#(yq-yc)*xdirc=(xq-xc)*ydirc,
	# xb=0.125*(xa+3*xp+3*xq+xc),
	# yb=0.125*(ya+3*yp+3*yq+yc)],[xp,yp,xq,yq])
	# yields to
	xa = float(za[0])
	ya = float(za[1])
	xdira = float(dira[0])
	ydira = float(dira[1])
	xb = float(zb[0])
	yb = float(zb[1])
	xc = float(zc[0])
	yc = float(zc[1])
	xdirc = float(dirc[0])
	ydirc = float(dirc[1])
	xp=(-((-xdira*ydirc-3*ydira*xdirc)*xa+4*xdira*xdirc*ya+xdira
	*(4*xdirc*yc-8*xdirc*yb-4*ydirc*xc+8*ydirc*xb))
	/(3*ydira*xdirc-3*xdira*ydirc))
	yp=(-(-4*ydira*ydirc*xa+(3*xdira*ydirc+ydira*xdirc)*ya+ydira
	*(4*xdirc*yc-8*xdirc*yb-4*ydirc*xc+8*ydirc*xb))
	/(3*ydira*xdirc-3*xdira*ydirc))
	xq=((-4*ydira*xdirc*xa+ydira*(8*xdirc*xb-xdirc*xc)
	+4*xdira*xdirc*ya+xdira*(4*xdirc*yc-8*xdirc*yb-3*ydirc*xc))
	/(3*ydira*xdirc-3*xdira*ydirc))
	yq=((-4*ydira*ydirc*xa+4*xdira*ydirc*ya
	+ydira*(3*xdirc*yc-4*ydirc*xc+8*ydirc*xb)
	+xdira*(ydirc*yc-8*ydirc*yb))/(3*ydira*xdirc-3*xdira*ydirc))
	return [[za],[(xp,yp),(xq,yq),zc]]
	
# parallel of a cubic bezier path p (s,cs,ce,e)
# where s = start of p, cs = starting control of p
# ce = ending control of p, e = end of p 
# in distance r at the right iff r>0 (else left)
def beziersidesegment(s,cs,ce,e,r):
	# compute the midpoint (t = 0.5) of the original path
	mid = (0.125*(s[0]+3*cs[0]+3*ce[0]+e[0]),
	0.125*(s[1]+3*cs[1]+3*ce[1]+e[1]))
	middir = (-s[0]-cs[0]+ce[0]+e[0],
	-s[1]-cs[1]+ce[1]+e[1])
	midside = pointright(mid,middir,r)
	if (cs[0]-s[0] == 0) and (cs[1]-s[1] == 0) \
	and (ce[0]-e[0] == 0) and (ce[1]-e[1] == 0): 
		# this means both control points lie on the start and end point
		# hence the path must be a straight line
		xstartdir = e[0]-s[0]
		ystartdir = e[1]-s[1]
		xenddir = e[0]-s[0]
		yenddir = e[1]-s[1]
	elif (cs[0]-s[0] == 0) and (cs[1]-s[1] == 0):
		# the first control point lies on the start point
		xstartdir = s[0]-2*cs[0]+ce[0] # proportional to second derivative
		ystartdir = s[1]-2*cs[1]+ce[1]
		xenddir = e[0]-ce[0]
		yenddir = e[1]-ce[1]
	elif (ce[0]-e[0] == 0) and (ce[1]-e[1] == 0):
		# the second control point lies on the end point
		xstartdir = cs[0]-s[0]
		ystartdir = cs[1]-s[1]
		xenddir = e[0]-2*ce[0]+cs[0]
		yenddir = e[1]-2*ce[1]+cs[1]
	else:
		xstartdir = cs[0]-s[0]
		ystartdir = cs[1]-s[1]
		xenddir = e[0]-ce[0]
		yenddir = e[1]-ce[1]
	start = ( s[0]+ystartdir
	/(xstartdir**2+ystartdir**2)**.5*r ,
	s[1]-xstartdir
	/(xstartdir**2+ystartdir**2)**.5*r )
	end = ( e[0]+yenddir
	/(xenddir**2+yenddir**2)**.5*r ,
	e[1]-xenddir
	/(xenddir**2+yenddir**2)**.5*r )
	return bezierinterpolate(start,(xstartdir,ystartdir),midside,end,(xenddir,yenddir))
	
# bezierjoin returns the join of two bezierpaths
def bezierjoin(first,second):
	if first[-1][-1] == second[0][0]: 
		# the last point of the first path equals the first
		# point of the second path
		return first + second[1:]
	else:
		return first + second

# bezierarc returns a bezierpath which is a part of a circle
# around center c with the radius r starting in direction s
# and ending in direction e
def bezierarc(c,s,e,r):
	if (vecangle(s,e) % 360) <= 90: # not more than a quarter circle
		mid = vecadd(c,vecscaleto(vecadd(vecnorm(s),vecscaleto(e,-1)),r))
		start = pointright(c,s,r)
		end = pointright(c,e,r)
		return bezierinterpolate(start,s,mid,end,e)
	else:
		middir = pointright((0,0),vecadd(vecnorm(s),vecscaleto(e,-1)),-1) 
		return bezierjoin(bezierarc(c,s,middir,r),
		bezierarc(c,middir,e,r))
		
# bezierreverse returns a reversed copy of the path p
def bezierreverse(p):
	reverse = [] 
	overlap = [] # overlapping list items
	l = len(p)
	for i in range (0,l):
		reverse.append([])
		for j in range (0,len(overlap)):
			reverse[i].append(overlap[j])
		reverse[i].append(p[l-1-i][-1])
		overlap = p[l-1-i][-2::-1]
	return reverse
	
# bezierfontforge takes a (list) path p
# and returns a fontforge contour
def bezierfontforge(p):
	c = fontforge.contour()
	c.moveTo(p[0][0][0],p[0][0][1])
	for i in range (1,len(p)):
		if len(p[i]) == 1: # this means a straight line
			c.lineTo(p[i][0][0],p[i][0][1])
		elif len(p[i]) == 3: # this means a cubic bezier
			c.cubicTo(p[i][0],p[i][1],p[i][2])
	return c

# bezierrightpath returns the right part of an outline of a 
# (list) path p when drawn with a circular pen of radius r
def bezierrightpath(p,r):
	outline = [] 
	for i in range(1,len(p)):
		if (len(p[i]) == 1) or (len(p[i]) == 3):
			if len(p[i]) == 3: # cubic bezier path
				if i == 1: # include initial point
					outline += beziersidesegment(p[i-1][-1],p[i][0],p[i][1],p[i][2],r)
				else:
					outline += beziersidesegment(p[i-1][-1],p[i][0],p[i][1],p[i][2],r)[1:]
				startdir = vecdiff(p[i][2],p[i][1]) # for the following arc
			elif len(p[i]) == 1: # straight line
				if i == 1: # include initial point
					outline += [[pointright(p[i-1][0],vecdiff(p[i][0],p[i-1][-1]),r)]]
				outline += [[pointright(p[i][0],vecdiff(p[i][0],p[i-1][-1]),r)]]
				startdir = vecdiff(p[i][0],p[i-1][-1]) # for the following arc
			# append a round joint:
			if i == len(p)-1: # last point of the path, so enddir is the reverse direction
				if len(p[i]) == 3:
					enddir = vecdiff(p[i][1],p[i][2]) # 
				else: # assume a straight line
					enddir = vecdiff(p[i-1][-1],p[i][0])
			elif len(p[i+1]) == 3: # cubic bezier path follows
				enddir = vecdiff(p[i+1][0],p[i][-1])
			else: # assume that a straight line follows
				enddir = vecdiff(p[i+1][0],p[i][-1])
			if not (abs(startdir[0]*enddir[1]-enddir[0]*startdir[1])<0.00001 \
			and (math.copysign(1,enddir[0]) == math.copysign(1,startdir[0])) \
			and (math.copysign(1,enddir[1]) == math.copysign(1,startdir[1]))): 
				# if not nearly same direction (reverse direction is okay)
				outline += bezierarc(p[i][-1],startdir,enddir,r)
	return outline
	
# beziercircularoutline returns the outline of a 
# (list) path p when drawn with a circular pen of diameter d
# the list has to be reversed, as fontforge determines outer = clockwise
def beziercircularoutline(p,d):
	return bezierreverse(bezierjoin(bezierrightpath(p,0.5*d),
	bezierrightpath(bezierreverse(p),0.5*d)))

# bezierhomogeneous transforms the bezier path p homogeneously with
# the transformation matrix m
def bezierhomogeneous(p,m):
	transformed = []
	for i in range(0,len(p)):
		transformed.append([])
		for j in range(0,len(p[i])):
			transformedpoint = homogeneous(m,p[i][j][0],p[i][j][1])
			transformed[i].append((transformedpoint[0],transformedpoint[1]))
	return transformed
	
# bezieroutline returns the outline of a (list) path p when drawn with a
# elliptical pen of a horizontal diameter dx and a vertical diameter dy
# rotated by the angle alpha
def bezieroutline(p,dx,dy,alpha):
	ca = math.cos(math.radians(alpha))
	sa = math.sin(math.radians(alpha))
	m = [ca,sa,-sa*dy/dx,ca*dy/dx] # transformation matrix (squeeze and rotate)
	return bezierhomogeneous(beziercircularoutline(
	bezierhomogeneous(p,invertmatrix(m)),dx),m)

closedpath = [[(0,0)],[(100,0),(200,100),(200,200)],[(0,200)],[(0,0)]]
openpath = [[(0,0)],[(100,0),(200,100),(200,200)],[(0,200)],[(0,150)],[(0,100)]]

font = fontforge.font()
font.createChar(65)
font.createChar(66)
font["A"].foreground += bezierfontforge(beziercircularoutline(openpath,30))
font["B"].foreground += bezierfontforge(bezieroutline(closedpath,40,20,30))

font.save("outliner.sfd")
