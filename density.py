import numpy as np
import math

import csv
import numpy as np
import matplotlib.pyplot as plt

fname = "logfile_dend94.csv"

def findRg(ml):
	rg2 = 0
	nm = len(ml)/3
	for i in range(nm):
		for j in range(i+1,nm):
			mi = ml[i*3:i*3+3]
			mj = ml[j*3:j*3+3]
			rg2 += np.sum( (mi-mj)*(mi-mj) )
	return rg2/(nm*nm)

with open(fname,'r') as cf:
    dst = csv.reader(cf)
    totaldist = np.zeros(94)
    trg2 = 0

    dc = 0
    for row in dst:
        if row[0] == "PARAMS":  #read in a parameter line
        	if dc > 0:
        		print trg2/dc
    			plt.hist(totaldist/dc,10)
    			plt.show()
    			dc = 0
    			trg2 = 0
    			totaldist = np.zeros(np.size(totaldist))
        	continue
        if row[0] == "C:":  #read in a calcluated properties lines
        	continue
        if row[0] == "D:":  #read in a particle data line
        	numparts = (len(row) - 1)/4
        	pos = np.zeros(numparts*3)
        	dist = np.zeros(numparts)
        	for i in range( numparts ): #load data
        		mn = int(row[i*4 + 1])
        		for j in range(3):
        			pos[mn*3 +j] = float(row[i*4+2  +j])
        	for i in range( numparts ): #calculate distances
        		x1 = pos[i*3:i*3+3]
        		dist[i] = np.sqrt( x1.dot(x1) )
        	totaldist = totaldist + dist
        	trg2 += findRg(pos)
        	dc += 1
