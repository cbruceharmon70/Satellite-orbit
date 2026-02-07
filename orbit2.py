# orbit2.py, B. Harmon 10/25/2019
# Numerical integration of satellite around Earth

import numpy as np
import matplotlib.pyplot as plt

# units are mks
G =  6.674e-11       # Univeral gravity const (N-m^2/kg^2)
Me = 5.9722e24       # Mass of Earth (kg)
re = 6.371e6         # Radius of Earth (m)

NPTS=119             # One orbit
dt=100
t=np.arange(NPTS)*dt
x=np.zeros(NPTS)
y=np.zeros(NPTS)
xdot=np.zeros(NPTS)
ydot=np.zeros(NPTS)
i=0
for _ in t:
    if _==0:
        x[i]=-re-1000000    # start 1000 km above surface
        xdot[i]=0
        y[i]=0
        ydot[i]=-8500       # counter-clockwise orbit
    else:
        r3 = (x[i-1]**2 + y[i-1]**2)**1.5    # r^3
        dxdot = dt * (-G * Me * x[i-1] / r3) # the ODE in x
        xdot[i] = xdot[i-1] + dxdot
        x[i] = x[i-1] + xdot[i]*dt
        dydot = dt * (-G * Me * y[i-1] / r3) # the ODE in y
        ydot[i] = ydot[i-1] + dydot
        y[i] = y[i-1] + ydot[i]*dt
    i += 1

# model the earth surface
thet=2*np.pi*np.arange(101)/100
xe=re*np.cos(thet)
ye=re*np.sin(thet)
# Center of the earth
xc=[0]; yc=[0]

fig=plt.figure(figsize=(10, 6))
ax=fig.add_axes([0.1,0.1,0.85,0.85])
ax.axis('equal')
ax.set(xlim=(-0.8e7,1.6e7),ylim=(-1.1e7,1.2e7))
ax.plot(x,y, color='blue')
ax.set_title('Satellite orbiting earth')
ax.set_xlabel('Distance in meters')
ax.set_ylabel('Distance in meters')
ax.plot(xe,ye, color='green')
ax.scatter(xc,yc, color='red')
ax.grid()
plt.show()
