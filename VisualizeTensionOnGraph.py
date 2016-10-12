import numpy as np
import csv
import math
import random
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from graph_tool.all import *
import colormaps

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

def find2Neighbors(v, adj):
	nls = []
	nn = 0
	for j in range( len(adj) ):
		if int(adj[v][j]) == 1:
			nn += 1
			nls.append(j)
	if nn == 2:
		return nls[0], nls[1]
	else:
		print "Has more than 2 neighbors"
		exit(0)

def buildGraphFromAdj(adj):
	el = [] #edge list
	g = Graph(directed=False)
	verts = [ g.add_vertex() for x in range(len(adj))]
	for i in range(len(adj)):
		for j in range( i+1,len(adj[i]) ):
			if int(adj[i][j]) == 1:
				te = g.add_edge(verts[i],verts[j])
				el.append(te)
	return g, el

def buildGraphFromAdjPhantoms(adj,r,dm):
	el = [] #edge list
	g = Graph(directed=False)
	verts = [ g.add_vertex() for x in range(r)]
	for i in range(r):
		for j in range( i+1, r ):
			if  dm[i][j] == 2:
				te = g.add_edge(verts[i],verts[j])
				el.append(te)
	return g, el

def setColors(edgeProp, edgeList, data):
	norm = colors.Normalize()
	ndata = norm(data) #normalize data
	mcolor = colormaps.viridis(ndata)
	for i,e in enumerate(edgeList):
		edgeProp[e] = mcolor[i]


#returns graph and edge color.  Takes adjacency matrix from topo.csv
#a True/False flag about to draw phantoms, and scalar data to color edges
def graphFromTopo(topmat, hasPhantoms, cdat):
	g1, edgeList = buildGraphFromAdj(topmat)
	dmat = graph_tool.topology.shortest_distance(g1)
	if hasPhantoms:
		N = len(topmat)
		r = (N+1)/2
		g, edgeList = buildGraphFromAdjPhantoms(topmat,r,dmat)
		if (N+1)%2 != 0:
			print "This does not have phantoms!"
			exit(0)
	else:
		g = g1
	ec = g.new_edge_property("vector<double>")
	setColors(ec, edgeList, cdat)
	return g, ec

hasPhantoms = True

topmat = readTopo("Dendrimer5g_Run1/topo.csv")

cdat = [random.random() for x in range(94)]

g,ec = graphFromTopo(topmat, hasPhantoms, cdat)

pos = sfdp_layout(g, cooling_step=0.99)	
graph_draw(g,pos,vertex_size=5, edge_pen_width=1.2,edge_color=ec, output_size=(1000, 500))



