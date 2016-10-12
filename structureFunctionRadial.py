import math
import numpy as np
import matplotlib.pyplot as plt


Ns= 100
Nf = 10
xa = np.linspace(0,math.pi*10,10000)

sq = np.zeros(np.size(xa))
for z in range(Nf):
	data = (np.random.random( (Ns,3) ) - 0.5)
	for i,u in enumerate(data):
		for j in range(i+1,len(data)):
			rij = np.sqrt( np.dot(data[i]-data[j],data[i]-data[j]) )
			sq += 2.0*np.sinc(rij*xa/math.pi)

plt.plot(xa,sq/np.size(data)/Nf)
plt.show()

