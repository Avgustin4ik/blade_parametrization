from scripts import *
#### MAIN ####
# ps - pressure side
# ss - suction side
points_pressure = read('c14bsp1.txt')
points_suction = read('c14bsp2.txt')
points_pressure = scaling(points_pressure,r1=2.6,r2=0.25,B=42.47)
points_suction = scaling(points_suction,r1=2.6,r2=0.25,B=42.47)
spline_suction = getSplineFromPoints(points_suction)
spline_pressure = getSplineFromPoints(points_pressure)
result_data = FindCamberPoints(spline_suction,spline_pressure,20,eps=1e-5,leftborderX=points_pressure[0][0],rightborder=points_pressure[0][-1],border=0.005)
# inletEdgePoints = np.vstack(( np.hstack((points_pressure[0][0:2],points_suction[0][0:10])) ,np.hstack((points_pressure[1][0:2],points_suction[1][0:10])) )) 
inletEdgePoints = np.vstack((points_pressure[0][0:5],points_pressure[1][0:5]))
result_inlet = FitInletEdge(inletEdgePoints[0],inletEdgePoints[1])
points_camber = [(i[0],i[1]) for i in result_data]
points_camber.insert(0,(0,0))
points_camber.append((1,0))
spline_camber = getSplineFromPoints(np.transpose(points_camber))
r1 = 2.6
r2 = 0.25
B = 42.47
sf = 1/(B-r1-r2)
p2up,p2down = FindTrailingEdgePoints(spline_camber,r2*sf,Vertex(1,0))

plt.figure()
# there are need a legend
plt.axis('equal')
plt.plot([0,1],[0,0],'red')
for i in range(len(result_data)):
    c = Vertex(result_data[i][0],result_data[i][1])
    x2 = result_data[i][2]
    r = Vertex(x2,spline_pressure(x2)).length(c)
    t = np.arange(0,2*np.pi,.01)
    x = [math.cos(i)*r + c.x for i in t]
    y = [math.sin(i)*r + c.y for i in t]
    plt.scatter(x2,spline_pressure(x2),marker='x',c='red')
    plt.plot(x,y,'g--')
plt.scatter(points_suction[0],points_suction[1])
plt.scatter(points_pressure[0],points_pressure[1])
x = np.arange(points_pressure[0][0],points_pressure[0][-1],0.01)
plt.plot(x, spline_pressure(x),'black',x, spline_suction(x),'black')
r = result_inlet[2]
t = np.arange(0,2*np.pi,.01)
x = [math.cos(i)*r + result_inlet[0] for i in t]
y = [math.sin(i)*r + result_inlet[1] for i in t]
scaleFactor = 1/(42.47-2.6-0.25)
print(r,2.6*scaleFactor)
r=0.25*scaleFactor
x2 = [math.cos(i)*r + 1 for i in t]
y2 = [math.sin(i)*r + 0 for i in t]
x3 = [math.cos(i)*r1*sf + 0 for i in t]
y3 = [math.sin(i)*r1*sf + 0 for i in t]
plt.plot(x,y,'r--',x2,y2,'r--',x3,y3,'b--')
k = spline_pressure.derivative(nu=1)(points_pressure[0][-1])
b = points_pressure[1][-1] - k*points_pressure[0][-1]
plt.plot([points_pressure[0][-1],1.1],[points_pressure[1][-1],k*1.1+b])
k = spline_suction.derivative(nu=1)(points_suction[0][-1])
b = points_suction[1][-1] - k*points_suction[0][-1]
plt.plot([points_suction[0][-1],1.1],[points_suction[1][-1],k*1.1+b])
x = np.arange(np.transpose(points_camber)[0][0],np.transpose(points_camber)[0][-1]+0.01,0.01)
plt.plot(x,spline_camber(x),'green','solid')
plt.scatter(p2down.x,p2down.y,marker='d',c='red')
plt.scatter(p2up.x,p2up.y,marker='d',c='red')
plt.show()