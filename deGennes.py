import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages

from scipy import interpolate
from scipy.integrate import odeint

import time
import sys
import os

vc = 0.5
e = vc

def phi(x,y):
        return a*np.power(x,c)+b*np.power(y,c)

def vectorfield(w, t):
	x, y = w
	f = [y, -6.0*e*x+6.0*vc*x*x*x - (2.0/t)*y] #[x', y']
	return f

x0 = 1.5
y0 = 1.5

abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 15.0
dt=0.001
numpoints = int(stoptime/dt)
t = np.linspace(0.1,stoptime,numpoints)

def solveTraj(w0, t):
    wsol = odeint(vectorfield, w0, t, atol=abserr, rtol=relerr)

    xs = []
    ys = []

    for s in wsol:
        xs.append(s[0])
        ys.append(s[1])
    xax = np.array(xs)
    yax = np.array(ys)
    return xax, yax

xiv = np.linspace(0.0, 1.0, 50)
yiv = np.linspace(-10.0, 0.0, 50)

for i in xiv:
    for j in yiv:
        xax, yax = solveTraj([i,j], t)
        vx, vy = vectorfield([i,j],1)
        plt.plot( xax, yax)

for i in xiv:
    for j in yiv:
        vx, vy = vectorfield([i,j],0.1)
        plt.quiver(i,j,vx,vy, headwidth=4, headlength=6,width=0.001, pivot='middle')
plt.title("c=0.5, b>a")
plt.show()


