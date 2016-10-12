import csv
import numpy as np
import matplotlib.pyplot as plt

fname = "corr94.csv"

vecs = []
vmag = []
dsq = []
rg2 = []

with open(fname,'r') as cf:
    dst = csv.reader(cf)
    for row in dst:
        dsq.append(float(row[0]))
        rg2.append(float(row[1]))
        tempvec = np.array( [float(row[2]), float(row[3]), float(row[4])] )
        tmag = np.sqrt(tempvec.dot(tempvec))
        vmag.append (tmag)
        vecs.append( tempvec/tmag )

plt.plot(dsq,'x' )
dsqa = np.mean(dsq)
ste = np.std(dsq)

rg2a = np.mean(rg2)
star = np.std(rg2)

print dsqa, ste
print rg2a, star
plt.show()

vcorr = np.zeros( len(vecs) ) #array to hold correlations

for step in range( int(len(vecs)/10) ): #len(vecs)
    #print "For step size" + str(step)+" go to "+ str( len(vecs) - step )
    numav = 0.0
    for i in range(len(vecs) - step):
        vcorr[step] += np.dot( vecs[i], vecs[i+step] )
        numav += 1.0
    vcorr[step] = vcorr[step]/numav

plt.plot(vcorr[0:200],'-x')
plt.show()
