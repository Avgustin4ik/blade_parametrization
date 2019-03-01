from scripts import *

k1 = 1
b1 = 3
# исходная точка
x0 =1
y0 = k1*x0+b1
# функция нижней линии
def f2(x):
      k = math.atan(-0.1)
      b = 0
      return k*x + b
kr = -1

# - 'hybr'             :ref:`(see here) <optimize.root-hybr>`
# - 'lm'               :ref:`(see here) <optimize.root-lm>`
# - 'broyden1'         :ref:`(see here) <optimize.root-broyden1>`
# - 'broyden2'         :ref:`(see here) <optimize.root-broyden2>`
# - 'anderson'         :ref:`(see here) <optimize.root-anderson>`
# - 'linearmixing'     :ref:`(see here) <optimize.root-linearmixing>`
# - 'diagbroyden'      :ref:`(see here) <optimize.root-diagbroyden>`
# - 'excitingmixing'   :ref:`(see here) <optimize.root-excitingmixing>`
# - 'krylov'           :ref:`(see here) <optimize.root-krylov>`
# - 'df-sane'          :ref:`(see here) <optimize.root-dfsane>`

br = y0 - kr*x0
def fun(x):
      # return [x[1]-kr*x[0]-br,
      # f2(x[2]) - math.atan(-0.1)*x[2] - 0,
      # sqrt(pow(x0 - x[0],2) + pow(y0 - x[1],2)) - sqrt(pow(x[2] - x[0],2) + pow(f2(x[2]) - x[1],2))]
      p1 = Vertex(x0,y0)
      p2 = Vertex(x[2],f2(x[2]))
      pr = Vertex(x[0],x[1])
      kr2 = -1/math.atan(-0.1)
      br2 = f2(x[2])-kr2*x[2]
      return [
            x[1]-kr*x[0]-br,
            x[1]-kr2*x[0]-br2,
            p1.length(pr)-p2.length(pr)
      ]
sol = scipy.optimize.root(fun,[0,1,1],method='hybr')
# xr,yr,x2 = scipy.optimize.fsolve(fun,(5,5,5))
print(sol.x)
xr = sol.x[0]
yr = sol.x[1]
x2 = sol.x[2]
r = sqrt(pow((x0-xr),2)+pow((y0-yr),2))
c = Vertex(xr,yr)

plt.figure()
plt.grid()
plt.axis('equal')
t = np.arange(0,2*np.pi,.01)
x = [math.cos(i)*r + c.x for i in t]
y = [math.sin(i)*r + c.y for i in t]
plt.plot(x,y)
x = np.arange(0,10,1)
y = k1*x + b1
YR = kr*x + br
plt.scatter(x,y)
plt.scatter(xr,yr)
plt.plot(x,y,x,f2(x),x,YR)
plt.show()



from readingBlade import inletEdgePoints
from scipy.spatial.distance import cdist
from scipy.optimize import fmin

def FitInletEdge(X,Y):
      X = inletEdgePoints[0]
      Y = inletEdgePoints[1]

      # Choose the inital center of fit circle as the CM
      xm = X.mean()
      ym = Y.mean()

      # Choose the inital radius as the average distance to the CM
      cm = np.array([xm,ym]).reshape(1,2)
      rm = cdist(cm, np.array([X,Y]).T).mean()

      # Best fit a circle to these points
      def err(vec):
            w = vec[0]
            v = vec[1]
            r = vec[2]
            pts = [np.linalg.norm([x-w,y-v])-r for x,y in zip(X,Y)]
            return (np.array(pts)**2).sum()

      return xf,yf,rf = scipy.optimize.fmin(err,[xm,ym,rm])  


      # write PARAMETRS
# RInlet = r1*scaleFactor
# ROutlet = r2*scaleFactor
# AngleInlet = math.degrees(math.atan(float(spline_camber.derivative(nu=1)(0.0))))
# AngleOutlet = abs(math.degrees(math.atan(float(spline_camber.derivative(nu=1)(1.0)))))
# omega2 = 0.0
# PARAMETRS = {   'Angle inlet':AngleInlet,
#                 'Angle outlet':AngleOutlet,
#                 'R inlet':r1,
#                 'R outlet':r2,
#                 'X Max':x_max,
#                 'Y max':y_max,
#                 'R Max': r_max,
#                 'XR Max': xr_max,
#                 'omega 1': omega1,
#                 'omega 2': omega2,
#                 'Angle Bend':deltaBend,
#                 'R Bend':RBend,
#                 'XR Bend':xr_bend}
# Writing(fileName,PARAMETRS)