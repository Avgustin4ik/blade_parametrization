import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import minimize
import scipy
from vertex import *

x = [i for i in range(0,11,1)]
R = 10
y = [sqrt(R*R - i*i) for i in range(0,11)]
#default spline
tck = interpolate.splrep(x,y)
spline = interpolate.BSpline(tck[0],tck[1],tck[2])
sx = np.arange(0,10.01,0.01)
sy = spline(sx)
dsx,dsy = sx,spline.derivative()(sx)
ddsx,ddsy = sx,spline.derivative(nu=2)(sx)
splCurvature = abs(ddsy)/pow(1+pow(dsy,2),3/2)
#parametric spline
tck2,u = interpolate.splprep([x,y])
px,py = interpolate.splev(sx,tck2)
dpx,dpy = interpolate.splev(sx,tck2,der=1)
ddpx,ddpy = interpolate.splev(sx,tck2,der=2)
psplCurvature = abs(ddpy)/pow(1+pow(dpy,2),3/2)
#from paper
n=len(x)
plotpoints = 100
k=3
knotspace = range(n)
knots = interpolate.InterpolatedUnivariateSpline(knotspace,knotspace,k=k).get_knots()
knots_full = np.concatenate(([knots[0]]*k,knots,[knots[-1]]*k))
tckX = knots_full,x,k
tckY = knots_full,y,k
splineX = interpolate.UnivariateSpline._from_tck(tckX)
splineY = interpolate.UnivariateSpline._from_tck(tckY)
tP = np.linspace(knotspace[0],knotspace[-1],plotpoints)
xP = splineX(tP)
yP = splineY(tP)
f = interpolate.interp1d(x,y,kind='cubic')
xnew = np.arange(0,10.01,0.1)
a = f._call_spline(xnew)
plt.figure()
plt.subplot(231)
plt.axis('equal')
plt.scatter(x,y)
plt.plot(sx,sy)
plt.subplot(232)
plt.axis('equal')
plt.plot(dsx,dsy)
plt.subplot(233)
plt.axis('equal')
plt.plot(ddsx,ddsy)
plt.plot(ddsx,splCurvature)
####    PARAMETRIC  ####
plt.subplot(234)
plt.axis('equal')
plt.scatter(x,y)
# plt.plot(px,py)
# plt.plot(xP,yP)
plt.plot(xnew,f(xnew))
plt.subplot(235)
plt.axis('equal')
plt.plot(dpx,dpy)
plt.subplot(236)
plt.axis('equal')
plt.plot(ddpx,ddpy)
plt.plot(ddpx,psplCurvature)
# plt.plot(u,[-1/i for i in ddpy])
plt.show()