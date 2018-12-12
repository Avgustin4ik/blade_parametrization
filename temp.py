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
#parametric spline
tck2,u = interpolate.splprep([x,y])
px,py = interpolate.splev(u,tck2)
dpx,dpy = interpolate.splev(u,tck2,der=1)
ddpx,ddpy = interpolate.splev(u,tck2,der=2)

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
# plt.plot(ddsx,[-1/i  for i in ddsy])
plt.subplot(234)
plt.axis('equal')
plt.scatter(x,y)
plt.plot(px,py)
plt.subplot(235)
plt.axis('equal')
plt.plot(u,dpy)
plt.subplot(236)
plt.axis('equal')
plt.plot(u,ddpy)
# plt.plot(u,[-1/i for i in ddpy])
plt.show()