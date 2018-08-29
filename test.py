import matplotlib.pyplot as plt
import math
from vertex import *

xTop = [0,1,2,3,4,5]
yTop = [4,4.5,5,5.5,6,6.5]
# xBottom = [-1,0,1,2,3,4]
# yBottom = [i * -1 for i in yTop]
xBottom = [i for i in range (-1,9)]
yBottom = [-0.5*i-3 for i in range (-1,9)]


class Function(object):
    topLine = []
    bottomLine = []
    initIndex = 0
    def __init__(self, topData, bottomData, initIndex):
        self.topLine = topData
        self.bottomLine = bottomData
        self.initIndex = initIndex
    def value(self, x):
        index = self.initIndex
        k1,b1 = getNormalLineKB(self.topLine[index],self.topLine[index+1],False)
        k2,b2 = getNormalLineKB(self.bottomLine[x],self.bottomLine[x+1])
        xr = (b2-b1)/(k1-k2)
        yr = k2*xr+b2
        r1 = Vertex(xr,yr).length(getMidPoint(self.topLine[index],self.topLine[index+1]))
        r2 = Vertex(xr,yr).length(getMidPoint(self.bottomLine[x],self.bottomLine[x+1]))
        return abs(r2-r1),r1,r2,Vertex(xr,yr)

def getLineKB(vertex1, vertex2):
    k = (vertex2.y - vertex1.y)/(vertex2.x - vertex1.x)
    b = vertex1.y - vertex1.x * k
    return k,b
def getNormalLineKB(vertex1, vertex2, isUp = True):
    point = getMidPoint(vertex1,vertex2)
    vec = getVectorFromPoints(vertex1,vertex2).normal2vector()
    if isUp == True and vec.y < 0:
        vec.reverse()
    if isUp == False and vec.y > 0:
        vec.reverse()
    temp = point.move(vec)
    return getLineKB(point,temp)
def getMidPoint(v1,v2):
    vec = getVectorFromPoints(v1,v2)
    vec.setLength(v1.length(v2)/2)
    return v1.move(vec)
kr1,br1 = getNormalLineKB(Vertex(xTop[2],yTop[2]),Vertex(xTop[3],yTop[3]), False)
kr2,br2 = getNormalLineKB(Vertex(xBottom[1],yBottom[1]),Vertex(xBottom[2],yBottom[2]),True)
# def f(x):
#     vL2 = Vertex(x2n[x],y2n[x])            
#     vR2 = Vertex(x2n[x+1],y2n[x+1])
#     bottomPoint = getMidPoint(vL2,vR2)
#     topPoint = getMidPoint(vL1,vR1)
#     k2,b2 = getLineKB(vL2,vR2)
#     xr = (b2-b1)/(k1-k2)
#     yr = k2*xr+b2
#     r2 = Vertex(xr,yr).length(bottomPoint)
#     r1 = Vertex(xr,yr).length(topPoint)
#     return abs(r2-r1),r1,r2,Vertex(xr,yr)

def plotline(k,b,x1,x2):
    x = []
    y = []
    x.append(x1)
    x.append(x2)
    y.append(x1*k + b)
    y.append(x2*k + b)
    return x,y

def findCentr(someFunction,x):
    f = someFunction
    eps = 1e-3
    a = 0
    b = len(xBottom)-1
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
topData = [Vertex(xTop[i],yTop[i]) for i in range(len(xTop))]
bottomData = [Vertex(xBottom[i],yBottom[i]) for i in range(len(xBottom))]
f = Function(topData,bottomData,2)
res,r1,r2,centr = findCentr(f.value,3)

plt.figure()
plt.grid()
plt.axis('equal')
plt.plot(xTop,yTop)
plt.plot(xBottom,yBottom)
plt.scatter(xTop,yTop)
plt.scatter(xBottom,yBottom)
plt.scatter(centr.x,centr.y)
n1x,n1y = plotline(kr1,br1,xTop[0],xTop[-1])
n2x,n2y = plotline(kr2,br2,xBottom[0],xBottom[-1])
plt.plot(n1x,n1y)
plt.plot(n2x,n2y)
plt.show()
