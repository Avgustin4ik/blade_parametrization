import matplotlib.pyplot as plt
import math
from vertex import *
def getLineKB(vertex1, vertex2):
    k = (vertex2.y - vertex1.y)/(vertex2.x - vertex1.x)
    b = vertex1.y - vertex1.x * k
    return k,b
def getMidPoint(v1,v2):
    vec = getVectorFromPoints(v1,v2)
    vec.setLength(v1.length(v2)/2)
    return v1.move(vec)
def f(x):
    vL2 = Vertex(x2n[x],y2n[x])            
    vR2 = Vertex(x2n[x+1],y2n[x+1])
    bottomPoint = getMidPoint(vL2,vR2)
    topPoint = getMidPoint(vL1,vR1)
    k2,b2 = getLineKB(vL2,vR2)
    xr = (b2-b1)/(k1-k2)
    yr = k2*xr+b2
    r2 = Vertex(xr,yr).length(bottomPoint)
    r1 = Vertex(xr,yr).length(topPoint)
    return abs(r2-r1),r1,r2,Vertex(xr,yr)
line2x = [-1,0,1,2,3,4]
line2y = [0,-0.5,-1,-1.5,-2,-2.5]
line1x = [0,1,2,3,4,5]
line1y = [4,4.5,5,5.5,6,6.5]
p2 = getMidPoint(Vertex(line2x[2],line2y[2]),Vertex(line2x[3],line2y[3]))
k2,b2 = getLineKB(Vertex(line2x[2],line2y[2]),Vertex(line2x[3],line2y[3]))
vectorR2 = getVectorFromPoints(Vertex(line2x[2],line2y[2]),Vertex(line2x[3],line2y[3])).normal2vector()
# vectorR2.reverse()
temp2 =  p2.move(vectorR2)
kr2,br2 = getLineKB(p2,temp2)
p1 = getMidPoint(Vertex(line1x[2],line1y[2]),Vertex(line1x[3],line1y[3]))
vectorR1 = getVectorFromPoints(Vertex(line1x[2],line1y[2]),Vertex(line1x[3],line1y[3])).normal2vector()
temp1 = p1.move(vectorR1)
kr1,br1 = getLineKB(p1,temp1)
def plotline(k,b,x1,x2):
    x = []
    y = []
    x.append(x1)
    x.append(x2)
    y.append(x1*k + b)
    y.append(x2*k + b)
    return x,y
def f(x):
    vL2 = Vertex(line2x[x],line2y[x])            
    vR2 = Vertex(line2x[x+1],line2y[x+1])
    bottomPoint = getMidPoint(vL2,vR2)
    topPoint = getMidPoint(Vertex(line1x[2],line1y[2]),Vertex(line1x[3],line1y[3]))
    # kr2,br2 = getLineKB(vL2,vR2)
    xr = (br2-br1)/(kr1-kr2)
    yr = kr2*xr+br2
    r2 = Vertex(xr,yr).length(bottomPoint)
    r1 = Vertex(xr,yr).length(topPoint)
    return abs(r2-r1),r1,r2,Vertex(xr,yr)
def findCentr():
    eps = 1e-3
    a = 0
    b = len(line2x)-1
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
centr = findCentr()
plt.figure()
plt.axis('equal')
plt.plot(line1x,line1y)
plt.plot(line2x,line2y)
plt.scatter(line1x,line1y)
plt.scatter(line2x,line2y)
plt.scatter(p2.x,p2.y)
x,y = plotline(kr2,br2,2.5,line1x[-1])
plt.plot(x,y)
res = f(2)
plt.scatter(centr[3].x,centr[3].y)
x,y = plotline(kr1,br1,2.5,line1x[-1])
plt.plot(x,y)
plt.show()
