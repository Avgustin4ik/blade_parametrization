from scripts import *
from tkinter import *
from tkinter import filedialog as fd
import numpy as np
import os
import matplotlib.pyplot as plt
from statistics import median

#### MAIN ####
# ps - pressure side
# ss - suction side
from iapws import IAPWS97
sat_steam=IAPWS97(P=1,x=1)                #saturated steam with known P
sat_liquid=IAPWS97(T=370, x=0)            #saturated liquid with known T
steam=IAPWS97(P=2.5, T=500)               #steam with known P and T
print(sat_steam.h, sat_liquid.h, steam.h) #calculated enthalpies
#### reading data block ####
filename = fd.askopenfilename()
# filename = r'W:\Geometry\C-9015А\С-9015Авог.txt'
points_ps = np.transpose(np.loadtxt(filename,skiprows=1,usecols=(0,1)))
filename = fd.askopenfilename()
# filename = r'W:\Geometry\C-9015А\С-9015Асп.txt'
points_ss = np.transpose(np.loadtxt(filename,skiprows=1,usecols=(0,1)))
# filename = r'W:\Geometry\C-9015А\С-9015Агеом.txt'
if median(points_ps[1])>median(points_ss[1]):
    temp = points_ps
    points_ps = points_ss
    points_ss = temp
filename = fd.askopenfilename()
CHORD,INSTALL_ANGLE,PITCH,R_INLET,R_OUTLET = np.loadtxt(filename,skiprows=1,usecols=(0,1,2,3,4))
inletEdgePoints = np.vstack((points_ps[0][0:3],points_ps[1][0:3]))
result_inlet = FitInletEdge(inletEdgePoints[0],inletEdgePoints[1])
R_INLET = result_inlet[2]
scalefactor = 1.0/(CHORD - R_INLET - R_OUTLET)
chord = CHORD * scalefactor
pitch = PITCH * scalefactor
r_inlet = R_INLET * scalefactor
r_outlet = R_OUTLET * scalefactor
#### data postprocessing ####
points_ps = scaling(points_ps,1.0,R_INLET,R_OUTLET,CHORD)
points_ss = scaling(points_ss,1.0,R_INLET,R_OUTLET,CHORD)
spline_ps = getSplineFromPoints(points_ps)
spline_ss = getSplineFromPoints(points_ss)

# detecting r_max,xr_max and camber line
res = FindCamberPoints(spline_ss,spline_ps,count=15,eps=1e-4,border=0.01,rightborder=0.9)
radius_array = [Vertex(res[i][0],res[i][1]).length(Vertex(res[i][2],float(spline_ps(res[i][2])))) for i in range(0,len(res))]
xr_array = [i[0] for i in res]
points_rDistrub = np.transpose([(xr_array[i],radius_array[i]) for i in range(0,len(radius_array))])
spline_rDisturb = getSplineFromPoints(points_rDistrub)
r_max = np.max(radius_array)
xr_max = float(np.transpose(res)[0][np.where(radius_array == r_max)])
yr_max = float(np.transpose(res)[1][np.where(radius_array == r_max)])
points_camber = [(i[0],i[1]) for i in res]
points_camber.insert(0,(0,0))
points_camber.append((1,0))
spline_camber = getSplineFromPoints(np.transpose(points_camber))

#### detecting x_max, y_max
points_camber = np.transpose([(x,float(spline_camber(x))) for x in np.arange(0,1,0.001)])
y_max = np.max(points_camber[1])
x_max = float(points_camber[0][np.where(points_camber[1] == y_max)])

#### detecting omega 1
p2up,p2down = FindTrailingEdgePoints(spline_camber,r_outlet,Vertex(1.0,0.0))
W1 = FindPoint(spline_ps,r_inlet)
angle = lambda x: math.degrees(math.atan(x)) if math.atan(x)>0 else math.degrees(math.atan(x)+math.pi/2.0)

# omega1= 2 * (angle(spline_camber.derivative(nu=1)(0.0))-angle(spline_ps.derivative(nu=1)(W1.x)))
k = spline_camber.derivative(nu=1)(0.0)
k_temp,b_temp=getLineKB(Vertex(0.0,0.0),W1)
omega1= 2 * (angle(spline_camber.derivative(nu=1)(0.0))-angle(k_temp))
print(spline_ps.derivative(nu=1)(W1.x))
# k2 = spline_ps.derivative(nu=1)(W1.x)
k2 = -1.0/k_temp
b2 = W1.y-k2*W1.x
rorate = lambda a : [[math.cos(a),math.sin(a)],[-math.sin(a),math.cos(a)]]
rotation_angle = (90.0-omega1/2.0)*2.0
rMatrix = rorate(math.radians(rotation_angle))
w = np.dot(rMatrix,[W1.x,W1.y])
W2 = Vertex(w[0],w[1])

# detecting BEND
angle_scale = math.degrees(math.atan2(R_INLET-R_OUTLET, (CHORD-R_OUTLET)-R_INLET))
angle_new = 180 - (INSTALL_ANGLE - angle_scale)
# t_opt = 0.75#RECOMMENDET
t = PITCH
if angle_new<0:
    angle_new = angle_new + 180
step_vec = Vector(math.cos(math.radians(angle_new)),math.sin(math.radians(angle_new)))
step_vec.setLength(pitch)
step_point = p2down.move(step_vec)
result_bend = FindBend(spline_ss,step_point)
bend_point = Vertex(result_bend[0],spline_ss(result_bend[0]))
dt_bend = result_bend[1]
k_bend = -1/dt_bend
b_bend = bend_point.y - k_bend*bend_point.x
f = lambda x : k_bend*x+b_bend - float(spline_camber(x))
x0 = 0.5
xr_bend = float(scipy.optimize.root(f,x0).x)
BendPoint = Vertex(xr_bend,float(spline_camber(xr_bend)))
RBend = BendPoint.length(bend_point)
dt_camber = float(spline_camber.derivative(nu=1)(BendPoint.x))
deltaBend = math.degrees(2*(abs(math.atan(dt_bend))-abs(math.atan(dt_camber))))
AngleInlet = abs(math.degrees(math.atan(float(spline_camber.derivative(nu=1)(0.0)))))
AngleOutlet = abs(math.degrees(math.atan(float(spline_camber.derivative(nu=1)(1.0)))))
omega2 = 0.0

PARAMETRS = {
    '0':x_max,
    '1':y_max,
    '2':AngleInlet,
    '3':AngleOutlet,
    '4':r_inlet,
    '5':r_outlet,
    '6':omega1,
    '7':omega2,
    '8':r_max,
    '9':xr_max,
    '10':RBend,
    '11':xr_bend,
    '12':deltaBend
}
filename = fd.asksaveasfilename()
Writing(filename,PARAMETRS)

plt.figure()
ax = plt.gca()

#### Добавление окружностей на график ####
InletEdge=plt.Circle((0.0,0.0),r_inlet,fill=False)
OutletEdge=plt.Circle((1.0,0.0),r_outlet,fill=False)
RMax = plt.Circle((xr_max,yr_max),r_max,fill=False)
BendCircle = plt.Circle((BendPoint.x,BendPoint.y),RBend,fill=False)

ax.add_artist(InletEdge)
ax.add_artist(OutletEdge)
ax.add_artist(RMax)
ax.add_artist(BendCircle)
for i in range(len(res)):
    ax.add_artist(plt.Circle((res[i][0],res[i][1]),radius_array[i],fill=False,color='green',ls='--'))


plt.scatter(points_ps[0],points_ps[1],color='red',marker='+') #### Цетр входной кромки
plt.scatter(points_ss[0],points_ss[1],color='blue',marker='+')  #### Центр выходной кромки
plt.scatter(xr_max,yr_max,color='black',marker='+') #### Центр R MAX

plt.scatter(0,0,marker='+')
plt.scatter(1,0,marker='+')
plt.scatter(x_max,y_max,color='black',marker='+',s=200)

plt.scatter(W1.x,W1.y,marker='o',s=100,color='red')
plt.scatter(W2.x,W2.y,marker='o',s=100,color='black')
plt.scatter(p2up.x,p2up.y,color='black',marker='o')
plt.scatter(p2down.x,p2down.y,color='black',marker='o')
plt.plot((-0.1,0.0),(k*-0.1,k*0.0))
plt.plot((-0.1,W1.x),(k2*-0.1+b2,W1.y))
plt.show()

print('succsses')
