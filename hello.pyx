from libc.stdio cimport fread,fopen,fclose,ferror,feof, printf
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os
from glob import glob
import scipy.interpolate as interp
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.colorbar
from matplotlib.mlab import bivariate_normal
from graph_tool.all import *
import colormaps
import csv
import Tkinter
import tkMessageBox
import matplotlib.image as mpimg
import pickle
#import FloydWarshall as fw
cimport numpy as np

cdef extern from "dtypes.h":
	struct log_line:
		int N
		double rog2_sum 
		double de_sum
		double rmax
		double density_sum[1000]  #This is defined in definitions.h.  If you change it there, must change by hand here!
		double spacing
		double sigma
		double epsilon
		unsigned long int rg2_n
		double rg2_mean
		double rg2_m2
		unsigned long int de_n 
		double de_mean
		double de_m2
		int edgeCount
		int segCount
		unsigned long int ulongcandle

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

def genColorMapNorm(tenmin,tenmax,colorScale="Log"):
	"""if colorScale == "Log":
		norm = colors.SymLogNorm(0.03,linscale=0.2, vmin=tenmin, vmax=tenmax)
	else:"""
	norm = colors.Normalize()
	return norm

def setColors(edgeProp, edgeList, data, tenmin, tenmax, colorScale="Log"):
	norm = genColorMapNorm(tenmin,tenmax)
	ndata = norm(data) #normalize data
	mcolor = colormaps.magma(ndata)
	for i,e in enumerate(edgeList):
		edgeProp[e] = mcolor[i]


#returns graph and edge color.  Takes adjacency matrix from topo.csv
#a True/False flag about to draw phantoms, and scalar data to color edges
def graphFromTopo(topmat, hasPhantoms, cdat, tenmin, tenmax, colorScale="Log"):
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
	setColors(ec, edgeList, cdat, tenmin, tenmax, colorScale=colorScale)
	return g, ec, edgeList


def colorFromKramer(g, edgeList):
	ec = g.new_edge_property("vector<double>")
	kmm = fw._findgkm(g)
	_t1,_t2,edgeListKramer,_t3 = fw.buildAdjacency(g)  #I only need edgelist
	kc = []
	for e in edgeList:
		ek = edgeListKramer[e.source(), e.target() ] #find Kramer Matrix edge
		kc.append(kmm[ek, ek]) #append diagonal element
	setColors(ec, edgeList, kc)
	return ec


#Number of points, mean, and sum squared x2 
def OLVCombine(a, b):
	n1 = a[0]
	m1 = a[1]
	m21 = a[2]

	n2 = b[0]
	m2 = b[1]
	m22 = b[2]

	delta = m2 - m1
	xb = (n1*m1 + n2*m2)/(n1+n2)
	M2X = m21 + m22 + delta*delta*(n1*n2/(n1+n2))
	#return M2X/(n1+n2 - 1)  #This is how you get variance from an OLV
	return [n1+n2, xb, M2X]

def readFile(fname, expectFrames):
	cdef log_line b

	cdef np.ndarray edgeMean #Edge Length
	cdef np.ndarray edgeN
	cdef np.ndarray edgeM2
	cdef np.ndarray edgeRMean #Edge distance
	cdef np.ndarray edgeRN
	cdef np.ndarray edgeRM2

	cdef int numMonomers
	cdef int numPhantoms
	cdef double denbuff[1000]

	cdef unsigned long int ulb  #variables for fread buffering
	cdef double dbb
	cdef int ib

	a = fopen(fname,"rb")
	Nt = 0
	ro2t = 0
	desumt = 0
	dent = np.zeros(1000) 
	dpc = 0.0
	c_spacing = 0
	c_epsilon = 0
	c_sigma = 0

	if sizeof(b.ulongcandle) == 8:
		ulongint = np.int64
	else:
		ulongint = np.int32

	#online variables
	olv_rg2 = [0,0,0]
	olv_des = [0,0,0]

	olv_edgelen = []
	olv_edgedist = []
	olv_seglen = []

	while feof(a) == 0:  #this reads multiple entries
		nr = 0
		nrd = 0
		nr += fread(&numMonomers, sizeof(numMonomers), 1, a)
		nr += fread(&numPhantoms, sizeof(numPhantoms), 1, a)
		nr += fread(&b.N, sizeof(b.N), 1, a)
		nr += fread(&b.rog2_sum, sizeof(b.rog2_sum), 1, a)
		nr += fread(&b.de_sum, sizeof(b.de_sum), 1, a)
		nr += fread(&b.rmax, sizeof(b.rmax), 1, a)
		nr += fread(&b.density_sum[0], sizeof(b.density_sum[0]), 1000, a)
		#Now load densities
		#try:
		#	allDensities
		#except:
		#	tallDensities = np.zeros( (numMonomers - numPhantoms, 1000), dtype=np.float64 )
		#	allDensities  = np.zeros( (numMonomers - numPhantoms, 1000), dtype=np.float64 )
		for idi in range(numMonomers - numPhantoms):
			nrd += fread(&denbuff, sizeof(denbuff[0]), 1000, a)			
			#for idj in range(1000):
			#	allDensities[idi,idj] = denbuff[idj]

		dsum = np.zeros(1000, dtype=np.float64)
		for ii in range(1000):
			dsum[ii] = b.density_sum[ii]

		#print np.sum(np.sum(allDensities,axis=0) - dsum)

		nr += fread(&b.spacing, sizeof(b.spacing), 1, a)
		nr += fread(&b.sigma, sizeof(b.sigma), 1, a)
		nr += fread(&b.epsilon, sizeof(b.epsilon), 1, a)

		nr += fread(&b.rg2_n, sizeof(b.rg2_n), 1, a)
		nr += fread(&b.rg2_mean, sizeof(b.rg2_mean), 1, a)
		nr += fread(&b.rg2_m2, sizeof(b.rg2_m2), 1, a)

		nr += fread(&b.de_n, sizeof(b.de_n), 1, a)
		nr += fread(&b.de_mean, sizeof(b.de_mean), 1, a)
		nr += fread(&b.de_m2, sizeof(b.de_m2), 1, a)

		nr += fread(&b.edgeCount, sizeof(b.edgeCount), 1, a) #read in number of edges
		nr += fread(&b.segCount, sizeof(b.segCount), 1, a) #read in number of edges

		if nr != 1015 + 2: #Read EOF //This number is the number of objects read from a
			break
		else:
			print b.edgeCount, nr

		edgeMean =  np.zeros(b.edgeCount, dtype=np.float64) #allocate edge data
		edgeM2 =  np.zeros(b.edgeCount, dtype=np.float64) #allocate edge data
		edgeN =  np.zeros(b.edgeCount, dtype=ulongint) #allocate edge data

		edgeRMean =  np.zeros(b.edgeCount, dtype=np.float64) #allocate edge data
		edgeRM2 =  np.zeros(b.edgeCount, dtype=np.float64) #allocate edge data
		edgeRN =  np.zeros(b.edgeCount, dtype=ulongint) #allocate edge data]

		segLMean =  np.zeros(b.segCount, dtype=np.float64) #allocate edge data
		segLM2 =  np.zeros(b.segCount, dtype=np.float64) #allocate edge data
		segLN =  np.zeros(b.segCount, dtype=ulongint) #allocate edge data

		label_trans = np.zeros(b.segCount, dtype=np.int32) #allocate edge data

		for i in range(b.edgeCount):
			nrd += fread(&ulb, sizeof(ulb), 1, a)
			edgeN[i] = ulb
			nrd += fread(&dbb, sizeof(dbb), 1, a)
			edgeMean[i] = dbb
			nrd += fread(&dbb, sizeof(dbb), 1, a)
			edgeM2[i] = dbb

		for i in range(b.edgeCount):
			nrd += fread(&ulb, sizeof(ulb), 1, a)
			edgeRN[i] = ulb
			nrd += fread(&dbb, sizeof(dbb), 1, a)
			edgeRMean[i] = dbb
			nrd += fread(&dbb, sizeof(dbb), 1, a)
			edgeRM2[i] = dbb

		for i in range(b.segCount):
			nrd += fread(&ulb, sizeof(ulb), 1, a)
			segLN[i] = ulb
			nrd += fread(&dbb, sizeof(dbb), 1, a)
			segLMean[i] = dbb
			nrd += fread(&dbb, sizeof(dbb), 1, a)
			segLM2[i] = dbb

		for i in range(b.segCount):
			nrd += fread(&ib, sizeof(ib), 1, a)
			label_trans[i] = ib

		if nrd != (b.edgeCount*6) + b.segCount*4 + (numMonomers - numPhantoms) * 1000: #Read EOF //This number is the number of objects read from a
			break
		if dpc == 0.0:
			c_spacing = b.spacing
			c_epsilon = b.epsilon
			c_sigma	  = b.sigma
			olv_rg2 = [b.rg2_n, b.rg2_mean, b.rg2_m2]  #init OLVs
			olv_des = [b.de_n, b.de_mean, b.de_m2]
			for i in range(b.edgeCount):  #Why did I keep uncertainty on everything?
				olv_edgelen.append( (edgeN[i], edgeMean[i], edgeM2[i]) )
				olv_edgedist.append( (edgeRN[i], edgeRMean[i], edgeRM2[i]) )
			for i in range(b.segCount):
				olv_seglen.append( (segLN[i], segLMean[i], segLM2[i]) )
		elif c_spacing != b.spacing or c_epsilon != b.epsilon or c_sigma != b.sigma:
			sys.stderr.write("You CANNOT mix different parameters in an Avedat file!\n")
			exit(0)
		else:
			olv_rg2 = OLVCombine(olv_rg2, [b.rg2_n, b.rg2_mean, b.rg2_m2] )  #combine OLVs
			olv_des = OLVCombine(olv_des, [b.de_n, b.de_mean, b.de_m2] )
			for i in range(b.edgeCount):  #Why did I keep uncertainty on everything?
				olv_edgelen[i] = OLVCombine(olv_edgelen[i],  [edgeN[i], edgeMean[i], edgeM2[i]] )
				olv_edgedist[i] = OLVCombine(olv_edgedist[i],  [edgeRN[i], edgeRMean[i], edgeRM2[i]] )
			for i in range(b.segCount):
				olv_seglen[i] = OLVCombine(olv_seglen[i], [segLN[i], segLMean[i], segLM2[i]] )

		Nt += b.N
		ro2t += b.rog2_sum
		desumt += b.de_sum
		dar = np.zeros(1000)
		#Calculate normalization
		nps = 0.0
		for i in range(1000):
			nps += 1.0*b.density_sum[i]
		
		for i in range(1000):
			dar[i] = b.density_sum[i]#/nps
		dent += dar
		#tallDensities += allDensities
		#print "CHCK", np.sum(dent - np.sum(tallDensities,axis=0) )
		dpc += 1.0
	if Nt == 0 or dpc == 0:
		print fname +" is bad"
		return ()
	if Nt != expectFrames:
		print fname + " is not the correct length!  Size is " + str(Nt)
		#os.remove(fname)
		# return ()
	#print "R^2 is " + str(ro2t/Nt)
	#print "dE/dL is " + str(desumt/Nt)
	xa = np.linspace(0,b.rmax,1000)
	rv = {"NumPoints":Nt, "NumLogEntries":dpc, "Rog2":ro2t/Nt, "DeDL":desumt/Nt, "XAxis":xa, "Density":dent, "Params":(b.spacing, b.epsilon,b.sigma), "deOLV":olv_des, "Rg2OLV":olv_rg2, "EdgeLength":olv_edgelen, "EdgeDist":olv_edgedist, "SegLen":olv_seglen, "labelTrans":label_trans, "EdgeCount":b.edgeCount, "AllDensities":None}
	#plt.plot(xa,dent/dpc)
	#plt.show()
	fclose(a)
	return rv

def myKey(e):
	return e["Params"][1]

def groupAve(ydata, xdata, md):
	cat = [ [xdata[0], [ 0 ]  ] ]
	for i in range(1,len(xdata)):  #cycle through data
		foundCat = False
		for j,e in enumerate(cat): #check categories
			if (np.abs(e[0] - xdata[i]) <= md):  #add
				e[1].append(i)
				e[0] = np.mean([ xdata[k] for k in e[1]])
				foundCat = True
				break  #cannot put data point in more than one category
		if not foundCat:  #make new categroy
			cat.append( [xdata[i], [ i ]  ] )

	caty = []  #y groupings
	catx = []
	for e in cat:
		caty.append( np.mean([ydata[k] for k in e[1]]) )
		catx.append( e[0] )
	return catx, caty


def deleteme():
	global _KeepLayout
	result = tkMessageBox.askquestion("Keep this layout?", "Are You Sure?", icon='warning')
	if result == 'yes':
		print "Keeping"
		_KeepLayout = True
	else:
		_KeepLayout = False


#This is a hard coded dirty hack, to fixed next time I change the file format
#of Avedat
def findZ1(mol_spacing):
	c = 4.0
	l = mol_spacing
	pref = 16.0*np.pi*np.pi*l*l*l*l/c*(1.0 - np.exp(-1.0*c) )
	return pref

def pickleGraph(xax, yax, file):
	graph = [[],[]]
	for x in xax:
		graph[0].append(x)
	for y in yax:
		graph[1].append(y)
	with open(file,"w") as fh:
		pickle.dump(graph, fh)


def test(files, tensiondb, topofile, hasphantoms, zeroTensionLength):
	global _KeepLayout
	_KeepLayout = False
	traw = np.loadtxt(tensiondb,delimiter=",")  #Load extension to tension function (computed numerically)
	extension = traw[:,0]
	tension = traw[:,1]
	tenMin = np.amin(tension)
	tenMax = np.amax(tension)
	extToTension = interp.interp1d(extension,tension) #create interpolating function
	extMin = np.amin(extension)
	extMax = np.amax(extension)
	plt.plot(extension, tension,'-x')
	plt.xlabel(r'<$R(T)$>/<$R(0)$>')
	plt.ylabel(r'Tension (F/kT)')
	plt.grid()
	plt.savefig("TensionExtension.png",dpi=300)
	plt.close()


	bigList = []
	for fn in glob(files):
		td = readFile(fn,50000)
		bigList.append(td)
	sbl = sorted(bigList,key=myKey)
	
	h_axis = np.array([e["Params"][1] for e in sbl])
	v_axis = np.array([e["DeDL"] for e in sbl])
	v2_axis = np.array([e["Rog2"] for e in sbl])

	"""for e in sbl:
		print e["Rg2OLV"]
		print e["deOLV"]"""

	with PdfPages("my_report_debug.pdf") as pdf:
		v_ebar = np.array([ np.sqrt(e["deOLV"][2]/(e["deOLV"][0]-1.0)/e["deOLV"][0]) for e in sbl ])
		plt.errorbar(h_axis, v_axis,yerr=v_ebar,fmt="-x")
		plt.title(r'<$\frac{\partial U_{pair}}{\partial \epsilon} $>')
		plt.xlabel(r'$\epsilon$')
		plt.ylabel(r'$\frac{F}{kT}$')
		#v_ebar = np.array([ np.sqrt(e["deOLV"][2]/(e["deOLV"][0]-1.0)) for e in sbl ])
		#plt.errorbar(h_axis, v_axis,yerr=v_ebar,fmt="x")
		pdf.savefig()
		plt.savefig("dude.png",dpi=300, format="png")
		plt.close()

		v2_ebar = np.array([ np.sqrt(e["Rg2OLV"][2]/(e["Rg2OLV"][0]-1.0)/e["Rg2OLV"][0]) for e in sbl ])
		plt.errorbar(h_axis, v2_axis,yerr=v2_ebar,fmt="-x")
		plt.title("Radius of Gyration")
		plt.xlabel(r'$\epsilon$')
		plt.ylabel(r'$(R_G/a)^2$')
		#v2_ebar = np.array([ np.sqrt(e["Rg2OLV"][2]/(e["Rg2OLV"][0]-1.0)) for e in sbl ])
		#plt.errorbar(h_axis, v2_axis,yerr=v2_ebar,fmt="x")
		pdf.savefig()
		plt.savefig("Rg2.png",dpi=300, format="png")
		plt.close()

		"""rax = [i[1] for i in e["EdgeDist"]]
		tax = [i[1] for i in e["EdgeLength"]]
		plt.plot(rax,tax,'*')
		pdf.savefig()
		plt.close()"""

		#calculate ideal free energy
		nseg = sbl[0]["EdgeCount"]
		el = sbl[0]["Params"][0]
		F0 = -nseg*np.log(findZ1(el))
		fe = integrate.cumtrapz(v_axis, h_axis, initial=0)
		plt.plot(h_axis,fe + F0,"-x")
		plt.title(r'Free Energy (measured from ideal)')
		plt.xlabel(r'$\epsilon$')
		plt.ylabel(r'$\frac{F}{kT}$')
		pdf.savefig()
		plt.savefig("FreeEnergy.png",dpi=300, format="png")
		pickleGraph(h_axis, fe+F0, "FreeEnergy.pickle")
		plt.close()

		for i,e in enumerate(sbl):
			if i % 4 == 0:
				plt.plot(e["XAxis"],e["Density"],label= str.format('{0:.3f}', e["Params"][1]) )
				ad = e["AllDensities"]
				plt.plot(e["XAxis"], np.sum(ad,axis=0) )
		#plt.ylim([0,2])
		#plt.xlim([0,7])
		plt.xlabel("Distance from Origin (units = 1.0)")
		plt.ylabel("Normalized monomer count")
		plt.legend()
		pdf.savefig()
		plt.savefig("Density.png",dpi=300, format="png")
		plt.close()

		#Create Tension plots
		adjmat = readTopo(topofile)

		n=10
		color=iter(cm.summer(np.linspace(0,1,n)))
		axe = plt.subplot(111)
		for i,e in enumerate(sbl):
			if i % 4 == 0:
				try:
					c = next(color)
				except StopIteration:
					color=iter(cm.summer(np.linspace(0,1,n)))
					c= next(color)
				rax = [i[1] for i in e["EdgeDist"]]
				try:
					_elx = [i[1] for i in e["SegLen"]]
					tax = []
					for kl in _elx: #convert distance to tension
						newten = extToTension(kl/zeroTensionLength)
						tax.append(newten)
					#tax = [extToTension(i) for i in _elx]
					lt = e["labelTrans"]
					edN = len(e["SegLen"])/e["EdgeCount"]
					taxn = np.zeros(e["EdgeCount"])
					for u,ui in enumerate(tax):  #Find average tension
						taxn[lt[u]] += tax[u]/(1.0*edN)
					tax = taxn
				except ValueError:
					print "Max range: " + str(extMin) +", " + str(extMax)
					for i in _elx:
						sys.stdout.write(str(i)+",")
						print extToTension(i)
				ax, ay = groupAve(tax, rax, 0.1)
				axe.plot(ax,ay,c=c, label=str.format('{0:.4f}', e["Params"][1]))
				axe.plot(rax,tax,'x',c=c)
				
		plt.title("Edge Tension vs Distance from Origin")
		plt.xlabel("Distance (units = 1.0)")
		plt.ylabel(r'Tension $\frac{F}{kT} $')
		box = axe.get_position()
		axe.set_position([box.x0, box.y0, box.width * 0.85, box.height])
		plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		pdf.savefig()
		plt.savefig("Tension.png",dpi=300, format="png")
		plt.close()

		pos = None #initialize layout
		revise = True
		for i,e in enumerate(sbl):
			rax = [i[1] for i in e["EdgeDist"]]
			#tax = [ extToTension(i[1]) for i in e["EdgeLength"]]
			_elx = [i[1] for i in e["SegLen"]]
			tax = []
			for kl in _elx: #convert distance to tension
				newten = extToTension(kl/zeroTensionLength)
				tax.append(newten)
			lt = e["labelTrans"]
			edN = len(e["SegLen"])/e["EdgeCount"]
			taxn = np.zeros(e["EdgeCount"])
			for u,ui in enumerate(tax):  #Find average tension
				taxn[lt[u]] += tax[u]/(1.0*edN)
			tax = taxn

			ebx = [ np.sqrt(i[2]/(i[0]-1.0)/i[0]) for i in e["EdgeLength"]]
			hbx = [ np.sqrt(i[2]/(i[0]-1.0)/i[0]) for i in e["EdgeDist"]]
			plt.errorbar(rax,tax, yerr=ebx, xerr=hbx,fmt='*',ecolor='g')
			pdf.savefig()
			plt.close()

			g, edgecols,edgelist = graphFromTopo(adjmat, hasphantoms, tax, np.amin(tax), np.amax(tax)) #make graph (vertices) plot
			tensionColorMap = []
			for edge in g.edges():  #generate tension color map
				color = []
				for c in edgecols[edge]:
					color.append(c)
				tensionColorMap.append( (int(edge.source()), int(edge.target()), color) )
			print "Dumping Pickle"
			with open("tensioncm_"+str.format('{0:.4f}', e["Params"][1])+".pickle","w") as fh:
				pickle.dump(tensionColorMap, fh)

			while revise:  #loop until you are happy with a layout
				p2 = sfdp_layout(g, cooling_step=0.99)	
				graph_draw(g,p2,vertex_size=5, edge_pen_width=1.5,edge_color=edgecols, output_size=(1000, 1000))
				top = Tkinter.Tk()
				B1 = Tkinter.Button(top, text = "Keep? Close Dialog if not", command = deleteme)
				B1.pack()
				top.mainloop()
				if _KeepLayout:
					pos = p2
					revise = False

			#output tension graph
			if pos == None:
				pos = sfdp_layout(g, cooling_step=0.999)	
			graph_draw(g,pos,vertex_size=0.5, edge_pen_width=3.5,edge_color=edgecols, output_size=(2000, 2000),output="graph_"+str.format('{0:.4f}', e["Params"][1])+".png")
			graph_draw(g,pos,vertex_size=0.5, edge_pen_width=1.5,edge_color=edgecols, output_size=(500, 500),output="graph_"+str.format('{0:.4f}', e["Params"][1])+".pdf")
			#Now add in Logaarithmic color bar
			gimd = matplotlib.image.imread( "graph_"+str.format('{0:.4f}', e["Params"][1])+".png" )
			plt.subplot(1, 1, 1)
			plt.imshow(gimd, norm=genColorMapNorm(np.amin(tax),np.amax(tax)), cmap=colormaps.magma)
			plt.colorbar()
			plt.savefig("graph_"+str.format('{0:.4f}', e["Params"][1])+".png",dpi=300, format="png")
			plt.close()
			#output gkm map
			#kramcolor = colorFromKramer(g,edgelist)
			#graph_draw(g,pos,vertex_size=0.5, edge_pen_width=1.5,edge_color=kramcolor, output_size=(500, 500),output="kramer_"+str.format('{0:.4f}', e["Params"][1])+".pdf")


			"""img = mpimg.imread("mytmp.pdf") #load png
			plt.imshow(img)
			plt.title(r'$\epsilon$ is ' + str.format('{0:.4f}', e["Params"][1]))
			pdf.savefig()
			plt.close()"""

			#str.format('{0:.4f}', e["Params"][1])
			pickleGraph(e["XAxis"],e["Density"], "density_"+str.format('{0:.4f}', e["Params"][1])+".pickle")
			xden1 = e["XAxis"][1:]
			xden2 = e["XAxis"][0:-1]
			xjac = 4.0/3.0*np.pi*(np.power(xden1,3.0) - np.power(xden2,3.0) )
			yden = e["Density"][1:]
			#plt.plot(e["XAxis"],e["Density"])
			plt.plot(xden1, yden/xjac)
			plt.xlabel("Distance from CM (units = 1.0)")
			plt.ylabel("Normalized monomer count")
			pdf.savefig()
			plt.close()
