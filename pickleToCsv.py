import pickle
import numpy as np
import matplotlib.pyplot as plt
import sys

pfile = sys.argv[1]

with open(pfile,"r") as fd:
	data = pickle.load(fd)

x = np.array(data[0])
y = np.array(data[1])
dr = x[1] - x[0]

div = 4.0*np.pi*x*x*dr
div[0] = 1.0

div = (np.power(x[1:],3.0) - np.power(x[0:-1],3.0))*(4.0/3.0)*np.pi

print div

plt.plot(x[1:],y[1:]/div)
plt.show()

