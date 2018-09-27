import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from vertex import *
def Reading(fileName):
    with open(fileName) as file:
        data = file.readlines()
        x = []
        y = []
        for i in data:
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
fileName = 'NACA4415.txt'
spline1 = 'spline1.txt'
spline2 = 'spline2.txt'
data = Reading(fileName)
data1 = Reading(spline2)
data2 = Reading(spline1)
x1 = np.array(data1[0])
y1 = np.array(data1[1])
x2 = np.array(data2[0])
y2 = np.array(data2[1])
xy1 = np.array([x1,y1])
xy2 = np.array([x2,y2])
sorting(xy1,0)
x1 = np.array(xy1[0])
y1 = np.array(xy1[1])
tck1 = interpolate.splrep(x1,y1)
x1n = np.arange(0, 1, 0.00001)
y1n = interpolate.splev(x1n, tck1, der=0)
sorting(xy2,0)
x2 = np.array(xy2[0])
y2 = np.array(xy2[1])
tck2 = interpolate.splrep(x2,y2)
x2n = np.arange(0,1,0.00001)
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
    tempPoint = Vertex(x,y).move(vec)
    return getLineKB(Vertex(x,y),tempPoint)
def half_divide_method(a, b, f):
    e = 1e-6
    x = (a + b) / 2
    while math.fabs(f(x)[0]) >= e:
        x = (a + b) / 2
        a, b = (a, x) if f(a)[0] * f(x)[0] < 0 else (x, b)
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
        return math.fabs(r2-r1),r1,r2,centr
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
        # plt.pause(0.5)
        plt.show()
    def setInialX(self,x):
        self.xInit = x
def dihotomia(a,b,f,dx):
    while math.fabs(a-b)>1e-6:
        x = (a+b)/2
        if f(x-dx)[0]<f(x+dx)[0]:
            b = x
        else:
            a = x
    x = (a+b)/2
    fm = f(x)
    return fm
    # main
xData = 0.6
dxData = 0.01
k1,b1 = getNormal(suctionSpline,xData)
xx,yy = plotline(k1,b1,xData-dxData,xData+dxData)
f = objectiveFunction(suctionSpline,pressureSpline,xData)
a = 0.0
b = 1.0
result = dihotomia(a,b,f.value,1e-3)

    # get array of camber line coords
# numberOfCircles = 10
# dataXarray = [i/10 for i in range(1,9,1)]
# result = []
# for i in dataXarray:
#     f.setInialX(float(i))
#     result.append(dihotomia(a,b,f.value,1e-5))
plt.figure()
plt.axis('equal')
plt.grid()
point = result[3]
# for i in result:
#     point = i[3]
plt.scatter(point.x,point.y)
plt.plot(x1n,y1n,'r')
plt.plot(x2n,y2n,'b')
plt.plot(xx,yy)
plt.show()
