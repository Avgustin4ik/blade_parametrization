import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import minimize
import scipy
from vertex import *
from main import pressureSpline,suctionSpline

#####   ploting circle  ####
def getCircle(centr,r):
    x = []
    y = []
    t = np.linspace(0,2*math.pi,100)
    for a in t:
        x.append(math.cos(a)*r + centr.x)
        y.append(math.sin(a)*r + centr.y)
    return x,y

    # inlet circle fitting
xInlet = np.arange(0.2,-0.01,-0.01)
yInlet = np.array(pressureSpline(xInlet))
yInlet = np.append(yInlet,np.array(suctionSpline(np.arange(0,0.21,0.01))))
xInlet = np.append(xInlet,np.arange(0,0.21,0.01))
dataInlet = np.array([xInlet,yInlet])

# sorting(dataInlet,0)
tckInlet,u = interpolate.splprep([dataInlet[0],dataInlet[1]],s=0,nest=-1)
uInlet = np.arange(0,1.01,0.01)
vInlet = interpolate.splev(uInlet,tckInlet)
dx,dy = interpolate.splev(uInlet,tckInlet,der=1)
ddx,ddy = interpolate.splev(uInlet,tckInlet,der=2)

####   TRUE CODE MASTER PIECE  ####
normals = []
for i in range(len(dx)):
    v = Vector(-dy[i],dx[i]).normalize().setLength(0.1)
    normals.append(v)
curvature = abs(dx*ddy-dy*ddx)/pow(dx*dx+dy*dy,3/2)
####    FINDING INLET CIRCLE    ####
r = []
r = [-1.0/i for i in curvature]
start = 0
end = len(uInlet)
eps = 1e-3
count_r = 3
def findDoublesIndexes(data,eps,index):
    value = data[index]
    result = []
    for i in range(0,len(data)):
        condition = abs(data[i]-value)
        if condition<eps:
            result.append(i)
    return result
result = []
for i in range(0,len(r)//2):
    a = findDoublesIndexes(r,eps,i)
    if len(a)>count_r:
        result.append(a)   
from main import normals,xy1,xy2
plt.figure()
for j in result:
    plt.clf()
    plt.subplot(121)

    plt.scatter([xy1[0][i] for i  in range(0,5)],[xy1[1][i] for i  in range(0,5)],marker='X',c='green',s=100)
    plt.scatter([xy2[0][i] for i  in range(0,10)],[xy2[1][i] for i  in range(0,10)],marker='X',c='orange',s=100)
    plt.axis('equal')
    plt.plot(vInlet[0],vInlet[1])
    plt.scatter(vInlet[0][len(r)//2],vInlet[1][len(r)//2],marker='x',c='red')
    plt.scatter(vInlet[0][j],vInlet[1][j])
    px = []
    py = []
    for i in j:
        normals[i].reverse()
        DX = [vInlet[0][i],vInlet[0][i]+normals[i].x]
        DY = [vInlet[1][i],vInlet[1][i]+normals[i].y]
        px.append(uInlet[i])
        py.append(r[i])
        plt.plot(DX,DY)
    index = j[0]
    R1 = r[index]
    DT = normals[index].setLength(R1)
    if DT.y<0:
        DT.reverse()
    Centr1 = Vertex(vInlet[0][index],vInlet[1][index]).move(DT)
    circle = getCircle(Centr1,R1)
    #### определение лучшего решения ####
    vec = normals[0].normal2vector()
    k1 = vec.y/vec.x
    b1 = xy1[1][0] - k1*xy1[0][0]
    temp = Vector(normals[0].y,normals[0].x)
    k2 = temp.y/temp.x
    b2 = Centr1.y - k2*Centr1.x
    # crossPoint = Vertex((b2-b1)/(k1-k2),k1*(b2-b1)/(k1-k2)+b1)
    # length = crossPoint.length(Centr1)

    plt.plot(circle[0],circle[1])
    plt.subplot(122)
    plt.scatter(px,py)
    plt.plot(uInlet,r)
    plt.plot([uInlet[len(r)//2],uInlet[len(r)//2]],[-1,1],'red')
    plt.draw()
    plt.waitforbuttonpress()
    

inletCircle = [Centr1,R1]
plt.figure()
plt.axis('equal')
plt.plot(uInlet,r)
plt.show()
####   TRUE CODE MASTER PIECE  ####
plt.figure()
plt.subplot()
plt.axis('equal')
plt.plot(vInlet[0],vInlet[1])
plt.figure()
plt.axis('equal')
plt.plot(vInlet[0],vInlet[1])
for i in range(len(dx)):
    x = [vInlet[0][i],vInlet[0][i]+normals[i].x]
    y = [vInlet[1][i],vInlet[1][i]+normals[i].y]
    plt.plot(x,y)
plt.figure()
plt.plot(uInlet,curvature)
plt.show()