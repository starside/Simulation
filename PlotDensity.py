import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys

dl = []

N = [299.0,328.0]

for filename in sys.argv[1:]:
	with open(filename,"r") as fp:
		dl.append( pickle.load(fp) )

for i,data in enumerate(dl):
	xax = np.array(data[0])
	yax = np.array(data[1])
	dx=xax[1]-xax[0]
	jax = 4.0*np.pi*xax*xax*dx
	s = 2
	nc = np.trapz(yax[1:]/jax[1:],xax[1:])
	print nc
	plt.plot(xax[s:],yax[s:]/jax[s:]/nc)

plt.show()
