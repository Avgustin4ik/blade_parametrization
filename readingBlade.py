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
def scaling(data,b=1.0):
    scaleFactor = b/max(data[0])
    return [i*scaleFactor for i in data]
def read(fileName):
        data = Reading(fileName)
        xy = np.array(data)
        return xy            
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
def getSplineFromPoints(pointsArray):
        data = sorting(pointsArray)
        tck = interpolate.splrep(data[0],data[1])
        spline = interpolate.BSpline(tck[0],tck[1],tck[2])
        return spline
#### MAIN ####
points_pressure = read('c14bsp1.txt')
points_suction = read('c14bsp2.txt')
spline_suction = getSplineFromPoints(points_suction)
spline_pressure = getSplineFromPoints(points_pressure)


plt.figure()
# there are need a legend
plt.axis('equal')
plt.scatter(points_suction[0],points_suction[1])
plt.scatter(points_pressure[0],points_pressure[1])
plt.show()

# spline1 = 'c14bsp1.txt'
# fileName = 'c14bsp2.txt'
# data2 = Reading(spline1)
# x1 = np.array(data1[0])
# y1 = np.array(data1[1])
# x2 = np.array(data2[0])
# y2 = np.array(data2[1])
# xy1 = np.array([x1,y1])
# xy2 = np.array([x2,y2])