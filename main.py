import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import minimize
import scipy
from vertex import *
def Reading(fileName):
    with open(fileName) as file:
        data = file.readlines()
        x = []
        y = []
        for i in data:
            a = i.split('\t')
            x.append(float(i.split('\t')[0]))
            y.append(float(i.split('\t')[1]))
    return [x,y]
def sorting(matrix,raw=0):
    if raw == 0:
        j = 1
    else:
            j = 1
    n = 1 
    while n < len(matrix[raw]):
        for i in range(len(matrix[raw])-n):
            if matrix[raw][i] > matrix[raw][i+1]:
                matrix[raw][i],matrix[raw][i+1] = matrix[raw][i+1],matrix[raw][i]
                matrix[j][i],matrix[j][i+1] = matrix[j][i+1],matrix[j][i]
        n += 1
# Reading files
# fileName = 'NACA4415.txt'
# spline1 = 'spline1.txt'
# spline2 = 'spline2.txt'
# data = Reading(fileName)
spline1 = 'c14bsp1.txt'
spline2 = 'c14bsp2.txt'
data1 = Reading(spline2)
data2 = Reading(spline1)
x1 = np.array(data1[0])
y1 = np.array(data1[1])
x2 = np.array(data2[0])
y2 = np.array(data2[1])
xy1 = np.array([x1,y1])
xy2 = np.array([x2,y2])
def scaling(data,b=1.0):
    scaleFactor = b/max(data[0])
    return [i*scaleFactor for i in data]
xy1 = scaling(xy1)
xy2 = scaling(xy2)
sorting(xy1,0)
x1 = np.array(xy1[0])
y1 = np.array(xy1[1])
tck1 = interpolate.splrep(x1,y1)
x1n = np.arange(0, 1, 0.001)
y1n = interpolate.splev(x1n, tck1, der=0)
sorting(xy2,0)
x2 = np.array(xy2[0])
y2 = np.array(xy2[1])
tck2 = interpolate.splrep(x2,y2)
x2n = np.arange(0,1,0.001)
y2n = interpolate.splev(x2n,tck2,der=0)
# get Splines and Derivative
suctionSpline = interpolate.BSpline(tck1[0],tck1[1],tck1[2])
pressureSpline = interpolate.BSpline(tck2[0],tck2[1],tck2[2])
dtSuctionSpline = suctionSpline.derivative()
dtPressureSpline = pressureSpline.derivative()
# Custom classes and functions
def getLineKB(vertex1, vertex2):
    k = (vertex2.y - vertex1.y)/(vertex2.x - vertex1.x)
    b = vertex1.y - vertex1.x * k
    return k,b
def plotline(k,b,x1,x2):
    x = []
    y = []
    x.append(x1)
    x.append(x2)
    y.append(x1*k + b)
    y.append(x2*k + b)
    return x,y
def getNormal(spline,x, isUp = False):
    y = spline(x)
    k = spline.derivative()(x)
    # if k>0 and isUp == False:
    #     k*=-1
    # if k<0 and isUp == True:
    #     k*=-1
    b = y - k*x
    dx = 0.01
    point1 = Vertex(x-dx,k*(x-dx)+b)
    point2 = Vertex(x+dx,k*(x+dx)+b)
    vec = getVectorFromPoints(point1,point2).normal2vector()
    if isUp == True and vec.y < 0:
        vec.reverse()
    if isUp == False and vec.y > 0:
        vec.reverse()
    tempPoint = Vertex(x,y).move(vec)
    return getLineKB(Vertex(x,y),tempPoint)
def half_divide_method(a, b, f):
    e = 1e-6
    x = (a + b) / 2
    while math.fabs(f(x)) >= e:
        x = (a + b) / 2
        a, b = (a, x) if f(a) * f(x) < 0 else (x, b)
    return (a + b) / 2    
class objectiveFunction(object):
    
    res = []
    def __init__(self,_suctionSpline,_pressureSpline,_xInit):
        self.suctionSpline = _suctionSpline
        self.pressureSpline = _pressureSpline
        self.xInit = _xInit
        self.initPoint = Vertex(_xInit, suctionSpline(_xInit))
    def value(self,x):
        k1,b1 = getNormal(self.suctionSpline,self.xInit)
        k2,b2 = getNormal(self.pressureSpline, x)
        point = Vertex(x,self.pressureSpline(x))
        xr = (b2-b1)/(k1-k2)
        yr = k2*xr+b2
        centr = Vertex(xr,yr)
        r1 = centr.length(self.initPoint)
        r2 = centr.length(point)
        self.centr = centr
        self.res.append(math.fabs(r1-r2))
        self.lastX = x
        return math.fabs(r2-r1),r1,r2,centr,self.lastX
    def justValue(self,x):
        return self.value(x)[0]
        
    def plot(self,x):
        #plotig
        plt.clf()
        plt.axis('equal')
        plt.grid()
        plt.plot(x1n,y1n,'r',x2n,y2n,'b')
        k1,b1 = getNormal(self.suctionSpline,self.xInit)
        k2,b2 = getNormal(self.pressureSpline, x)
        dxData = 0.01
        xx1,yy1 = plotline(k1,b1,self.xInit-dxData,self.xInit+dxData)
        xx2,yy2 = plotline(k2,b2,x-xData*0.1,x+xData*0.1)
        plt.plot(xx1,yy1,xx2,yy2)
        plt.scatter(self.centr.x,self.centr.y)
        plt.pause(0.001)
        # plt.show()
    def setInialX(self,x):
        self.xInit = x
def dihotomia(a,b,f,dx):
    while math.fabs(a-b)>1e-6:
        x = (a+b)/2
        if f.value(x-dx)[0]<f.value(x+dx)[0]:
            b = x
        else:
            a = x
        # f.plot(x)
    x = (a+b)/2
    fm = f.value(x)
    return fm,x
    # main
def der(f,x):
    dx = 0.01
    return (f(x+dx)-f(x))/dx
def newtons_method(x0, f, e):
    #f1 - производная
    dx = 0.001
    def f1(x):
        return (f(x+dx)-f(x))/(dx)
    x0 = float(x0)
    while True:
        x1 = x0 - (f(x0) / f1(x0))
        if abs(x1 - x0) < e:
            return x1
        x0 = x1    
xData = 0.1
dxData = 0.01
k1,b1 = getNormal(suctionSpline,xData)
xx,yy = plotline(k1,b1,xData-dxData,xData+dxData)
f = objectiveFunction(suctionSpline,pressureSpline,xData)
a = 0.0
b = 1.0
# result,x = dihotomia(0.0,1.0,f,0.001)

# result = half_divide_method(a,b,f.justValue)
dataXarray = [i/10 for i in range(1,9,1)]
points = []
for i in dataXarray:
    x0 = i
    f = objectiveFunction(suctionSpline,pressureSpline,x0)
    
    res = minimize(f.justValue,x0, method='nelder-mead',options={'xtol':1e-5,'disp':True})
    # result = newtons_method(x0,f.justValue,1e-8)
    points.append(f.centr)
# def removeError(data):
#     for i in range(len(data)):
#         if data
# x = f.lastX
# k2,b2 = getNormal(pressureSpline, x)
# xxp,yyp = plotline(k2,b2,x-dxData,x+dxData)
plt.figure()
plt.axis('equal')
plt.grid()
# point = f.centr
for i in points:
    plt.scatter(i.x,i.y, c="green", marker='x')
# plt.scatter(point.x,point.y)
# plt.scatter(x,pressureSpline(x))
plt.plot(x1n,y1n,'r')
plt.plot(x2n,y2n,'b')
plt.plot(xx,yy)
# plt.plot(xxp,yyp)
plt.show()
