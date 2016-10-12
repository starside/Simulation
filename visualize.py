import numpy as np
import csv
import math
import matplotlib.pyplot as plt
from visual import *
from matplotlib import cm
import pickle

def readTopo(fname):
	adj = None
	with open(fname,'r') as csvfile:
		tr = csv.reader(csvfile, delimiter=',')
		myr = []
		for i,row in enumerate(tr):
			if adj == None:
				adj = np.zeros( (len(row),len(row) ) )
			for col in row:
				myr.append(int(col))
			for j,u in enumerate(myr):
				adj[i][j] = u
			myr = []
	return adj

def structureFunction(xa,data):
	dl = len(data)
	npoints = np.size(xa)
	sq = np.zeros(npoints)
	for i,u in enumerate(data):
		for j in range(i+1,dl):
			rij = np.sqrt( np.dot(data[i]-data[j],data[i]-data[j]) )
			sq += 2.0*np.sinc(rij*xa/math.pi)/(1.0*dl)
	return 1.0 + sq

def plot3DMol(data, sigma, adj, tenmap):
	global ball_list
	for b in ball_list:
		b.visible = False
		del b
	#draw links
	for i in range(len(data)):
		for j in range(i+1, len(data)):
			if adj[i][j] == 1:
				if i <= phantomStart:
					ptcolor = (tenmap[i][0], tenmap[i][1], tenmap[i][2])
				else:
					ptcolor = color.green
				maxis = data[j] - data[i]
				s = cylinder(pos=data[i],axis=maxis,radius=0.05, color=ptcolor )
				ball_list.append(s)
	#Draw monomers
	for i,p in enumerate(data):
		global phantomStart
		if i <= phantomStart:
			try:
				ptcolor = (tenmap[i][0], tenmap[i][1], tenmap[i][2])
			except KeyError:
				ptcolor = color.green
			ball = sphere(pos=vector(p[0],p[1],p[2]), radius=sigma/2,color=ptcolor, opacity=0.02)
			ball_list.append(ball)

def findRg(ml):
	rg2 = 0
	nm = len(ml)
	for i in range(nm):
		for j in range(i+1,nm):
			rg2 += np.sum( (ml[i]-ml[j])*(ml[i]-ml[j]) )
	return rg2/(nm*nm)

def readHeader(fh):
	head = np.fromfile(fh, dtype=np.uint32, count=1)
	seed = np.fromfile(fh, dtype=np.uint64, count=1)
	zCount  = np.fromfile(fh, dtype=np.uint32, count=1)
	numMonomers = np.fromfile(fh, dtype=np.int32,count=1)
	if head[0] != 0xdeadbeef:
		print "Bad header!"
		exit(0)
	return head, seed,zCount,numMonomers

def readLine(fh, numMonomers):
	head = np.fromfile(fh,dtype=np.int8, count=1)
	if chr(head[0]) == 'P':
		spacing = np.fromfile(fh,dtype=np.float64,count=1)
		sigma = np.fromfile(fh,dtype=np.float64,count=1)
		epsilon = np.fromfile(fh,dtype=np.float64,count=1)
		return "P", spacing, sigma, epsilon
	elif chr(head[0]) == 'D':
		ml = []
		zCount  = np.fromfile(fh, dtype=np.uint32, count=1) #read in the number of random calls at this frame
		monomers = np.fromfile(fh, dtype=np.float64, count = 3*numMonomers)
		for i in range(numMonomers):
			ml.append(monomers[i*3:i*3+3])
		return "D", zCount, ml
	else:
		print "Screwed up file format!" + chr(head[0])
		exit(0)

num_frames = 15#00  #number of frames in the file
maximum_q = 10.0
q_resolution = 10000

fp = open("tmp/QBS_299_7_8146527361320424505.tf","rb")

head, seed, zCount, numMonomers = readHeader(fp)
print numMonomers

dataCount = 0
paramCount = 0
rg2a = 0.0
qax = np.linspace(0,maximum_q,q_resolution)
sqa = np.zeros(q_resolution)
my_sigma = 0.0
phantomStart = 327

ball_list = []

topmat = readTopo("QBeta/Len327/topo.csv")

tenmap = {}
with open("tensioncm_0.9200.pickle","r") as fh:
	tensionmap = pickle.load(fh)

for ic in tensionmap:
	scolor = ic[2]
	tenmap[ic[0]] = scolor
	tenmap[ic[1]] = scolor

scene2 = display(title='Examples of Tetrahedrons',
     x=0, y=0, width=800, height=600,
     center=(5,0,0), background=(1,1,1))

while True:
	try:
		res = readLine(fp, numMonomers[0])
		print res[0]
		if res[0] == 'D':
			dataCount += 1
			rg2a += findRg(res[2])
			#sqa += structureFunction(qax,res[2])/(1.0*5)
			print "Completed "+str(dataCount)+"/"+str(num_frames)
			plot3DMol(res[2], my_sigma, topmat, tenmap)
			raw_input("Press Enter to continue")
		elif res[0] == 'P':
			paramCount += 1
			my_sigma = res[2][0]
		else:
			print "Unknown line type "+str(res[0])
			exit(0)
		if dataCount == 5:
			break
	except IndexError:  #Sloppy but it works
		break

if num_frames != dataCount:
	print "Error:  Did not read expected number of frames"
print "Read "+str(dataCount)+" configurations"
print "Read "+str(paramCount)+ " different parameter sets"
rg2a = rg2a / (1.0*dataCount)
print rg2a

plt.plot(qax, sqa/(1.0*numMonomers) )

qas = np.linspace(0,1.0/rg2a,100)
plt.plot(qas, 1.0 - qas*qas*(rg2a/3.0), 'x' )
plt.show()

fp.close()

