import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from vertex import *
# import mathscipy
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


def getLineKB(vertex1, vertex2):
    k = (vertex2.y - vertex1.y)/(vertex2.x - vertex1.x)
    b = vertex1.y - vertex1.x * k
    return k,b
fileName = 'NACA4415.txt'
spline1 = 'spline1.txt'
spline2 = 'spline2.txt'
data = Reading(fileName)
data1 = Reading(spline1)
data2 = Reading(spline2)
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
BSpline = interpolate.BSpline(tck1[0],tck1[1],tck1[2])
a = BSpline.__call__([0.5])
a = BSpline(0.5)
dtSpline = BSpline.derivative(1)
dt = dtSpline(0.5)
temp = interpolate.splev([0.5],tck1)

x1n = np.arange(0, 1, 0.00001)
y1n = interpolate.splev(x1n, tck1, der=0)
sorting(xy2,0)
x2 = np.array(xy2[0])
y2 = np.array(xy2[1])
tck2 = interpolate.splrep(x2,y2)
x2n = np.arange(0,1,0.00001)
y2n = interpolate.splev(x2n,tck2,der=0)
def getMidPoint(v1,v2):
    vec = getVectorFromPoints(v1,v2)
    vec.setLength(v1.length(v2)/2)
    return v1.move(vec)
t = 0.5
vL1 = Vertex(x1n[49],y1n[49])
vR1 = Vertex(x1n[50],y1n[50])
topPoint = getMidPoint(vL1,vR1)
vector1 = getVectorFromPoints(vL1,vR1).normal2vector()
temp1 = topPoint.move(vector1)
kr1,br1 = getLineKB(topPoint,temp1)

def f(x):
    vL2 = Vertex(x2n[x],y2n[x])            
    vR2 = Vertex(x2n[x+1],y2n[x+1])
    bottomPoint = getMidPoint(vL2,vR2)
    # topPoint = getMidPoint(Vertex(x2n[2],y2n[2]),Vertex(y2n[3],y1n[3]))
    vector2 = getVectorFromPoints(vL2,vR2).normal2vector().reverse()
    temp2 = bottomPoint.move(vector2)
    kr2,br2 = getLineKB(bottomPoint,temp2)
    # kr2,br2 = getLineKB(vL2,vR2)
    xr = (br2-br1)/(kr1-kr2)
    yr = kr2*xr+br2
    r2 = Vertex(xr,yr).length(bottomPoint)
    r1 = Vertex(xr,yr).length(topPoint)
    return abs(r2-r1),r1,r2,Vertex(xr,yr)

    
def findCentr():
    eps = 1e-3
    a = 0
    b = len(x2n)-1
    x = math.floor((a+b)/2)
    if f(x)[0]<eps:
        return f(x)
    if abs(f(x)[0]*f(a)[0] < 0):
        b = x
    else:
        a = x
    xn = math.floor((a+b)/2)
    while abs(f(x)[0] - f(xn)[0])>eps:
        x = xn
        if f(x)[0]<eps:
            return f(x)
        if abs(f(x)[0]*f(a)[0] < 0):
            b = x
        else:
            a = x
        xn = math.floor((a+b)/2)
    return f(xn) 
result = findCentr()
# result = f(50)

plt.figure()
plt.axis('equal')
def plotline(k,b,x1,x2):
    x = []
    y = []
    x.append(x1)
    x.append(x2)
    y.append(x1*k + b)
    y.append(x2*k + b)
    return x,y
xx1,yy1 = plotline(kr1,br1,0,1)
plt.scatter(topPoint.x,topPoint.y)

# xx2,yy2 = plotline(kr2,br2,0,1)
# plt.plot(xx1,yy1)
# plt.scatter(x1,y1)
# plt.scatter(x2,y2)
# plt.scatter(result[3].x,result[3].y)
# plt.plot(x1n, y1n)
# plt.plot(x2n, y2n)
# plt.show()

z = 0