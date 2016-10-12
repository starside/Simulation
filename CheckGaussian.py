import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("build/foo.csv",delimiter=",")
norm = np.trapz(data[:,1], data[:,0])

y = data[:,1]/norm
x = data[:,0]

histmon = 2
a = 1.0135
#ac = a + 0.4
sigma = a/np.sqrt(3.0)

plt.plot(x,y,"-x")
#plt.plot(1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2))
plt.plot(x,1/(sigma*np.sqrt(2*np.pi)) * np.exp(-x*x/(2*sigma**2)) )
#sigma = ac/np.sqrt(3)
#plt.plot(x,1/(sigma*np.sqrt(2*np.pi)) * np.exp(-x*x/(2*sigma**2)) ,'--')
plt.title("Monomer "+str(histmon)+", a = "+str(a))
plt.show()