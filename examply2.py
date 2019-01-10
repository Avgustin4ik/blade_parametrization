from numpy import *
from matplotlib.pyplot import *
from scipy import interpolate

# first create a list of x,y points based on an airfoil design from
# https://en.wikipedia.org/wiki/NACA_airfoil
# created so that points are not uniformly distributed
c=1.0
th=0.15
x=concatenate( (linspace(1.0,0.1,40),linspace(0.95,0.05,20)**3.0 *0.1 ) )
y=5.0*th*(0.2969*np.sqrt(x/c) - 0.1260*(x/c) - 0.3516*(x/c)**2 +0.2843*(x/c)**3 - 0.1015*(x/c)**4)

x=concatenate( (x,flipud(x)))
y=concatenate( (y,-flipud(y)))


# create parameters. t is uniformly spaced on the interval 0 to 1. 
# u is spaced proportional to edge lengths
t=arange( len(x),dtype=float) / (len(x)-1)
u=zeros( (len(x)),dtype=float)
L=0.0
u[0]=0.0
for ii in range(1,len(x)):
  L += sqrt( (x[ii]-x[ii-1])**2 + (y[ii]-y[ii-1])**2 )
  u[ii]=L
u /= L

# create two cubic splines for each x and y
csxt=interpolate.splrep(t,x,k=3)
csxu=interpolate.splrep(u,x,k=3)
csyt=interpolate.splrep(t,y,k=3)
csyu=interpolate.splrep(u,y,k=3)

# tt and uu are fine samples used for plotting
tt=linspace(0.0,1.0,10000)
uu=linspace(0.0,1.0,10000)


# evaluate x(t),y(t),x(u),y(u) and their derivatives
xt0=interpolate.splev(tt,csxt)
xt1=interpolate.splev(tt,csxt,der=1)
xt2=interpolate.splev(tt,csxt,der=2)
yt0=interpolate.splev(tt,csyt)
yt1=interpolate.splev(tt,csyt,der=1)
yt2=interpolate.splev(tt,csyt,der=2)

xu0=interpolate.splev(uu,csxu)
xu1=interpolate.splev(uu,csxu,der=1)
xu2=interpolate.splev(uu,csxu,der=2)
yu0=interpolate.splev(uu,csyu)
yu1=interpolate.splev(uu,csyu,der=1)
yu2=interpolate.splev(uu,csyu,der=2)

# calculate curvature by the formula
# |x' y'' - y' x'' | / |x'^2 + y'^2|^3/2
Kt=abs( xt1*yt2 - yt1*xt2) / sqrt(xt1*xt1 + yt1*yt1)**3
Ku=abs( xu1*yu2 - yu1*xu2) / sqrt(xu1*xu1 + yu1*yu1)**3


# interpolate between t and u, so that we can plot
# apples to apples
cstu=interpolate.splrep(t,u,k=1)
uu2=interpolate.splev(tt,cstu)

# plots
ff=figure(1)
ff.clf()
plot(x,y,'.',markersize=6)
gca().set_aspect('equal')
grid('on')
xlabel('x')
ylabel('y')
axis([-0.05,1.05,-0.2,0.2])

ff=figure(2)
ff.clf()
subplot(211)
plot(tt,xt0,'b-',label=r"$x_t(t)$")
plot(tt,yt0,'r-',label=r"$y_t(t)$")
grid('on')
xlabel('t')
legend(loc='upper right',fontsize=20)
subplot(212)
plot(uu,xu0,'b-',label=r"$x_u(u)$")
plot(uu,yu0,'r-',label=r"$y_u(u)$")
grid('on')
xlabel('u')
legend(loc='upper right',fontsize=20)



ff=figure(3)
ff.clf()
subplot(211)
plot(tt,xt2,'b-',label=r"$x_t''(t)$")
plot(tt,yt2,'r-',label=r"$y_t''(t)$")
grid('on')
xlabel('t')
legend(loc='upper right',fontsize=20)
subplot(212)
plot(uu,xu2,'b-',label=r"$x_u''(u)$")
plot(uu,yu2,'r-',label=r"$y_u''(u)$")
grid('on')
xlabel('u')
legend(loc='upper right',fontsize=20)

ff=figure(4)
ff.clf()
plot(uu2,Kt,'b-',label=r"$\kappa_t(u)$")
plot(uu,Ku,'r-',label=r"$\kappa_u(u)$")
grid('on')
xlabel('u')
gca().set_xlim([0.4,0.6])
legend(loc='upper right',fontsize=20)