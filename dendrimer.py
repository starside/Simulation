import math
import random
import numpy as np

def _randomVector():
	phi = random.random()*2*math.pi
	cod = random.uniform(-1.0,1.0)
	theta = math.acos(cod)
	return np.array( [math.sin(theta)*math.cos(phi), math.sin(theta)*math.sin(phi), math.cos(theta)] )

def randomVector():
	return _randomVector()*1.5/2.0 + _randomVector()*1.5/2.0

def addMonomers(monList, f, g, x):
	if g == 1:
		return
	for i in range(f-1):
		newx = x + randomVector() #random vector
		monList.append(newx)
		addMonomers(monList,f,g-1, newx)

f = 3 #functionality
g = 5 #generations
samples = 100

def dendrimerWalk(f,g):
	monomers = []
	#initialize dendrimer
	monomers.append(np.array([0,0,0]))
	for i in range(f):
		monomers.append( randomVector() )
	#Generate generations
	for i in range(f):
		addMonomers(monomers,f,g,monomers[i+1])
	return monomers

def findRg(ml):
	rg2 = 0
	nm = len(ml)
	for i in range(nm):
		for j in range(i+1,nm):
			rg2 += np.sum( (ml[i]-ml[j])*(ml[i]-ml[j]) )
	return rg2/(nm*nm)

erg = 0.0
for i in range(samples):
	mons = dendrimerWalk(f,g)
	erg += findRg(mons)

print erg/samples