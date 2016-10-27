import numpy as np
import os
import os.path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f","--func", type=str, help="Specify the functionality of the monomers in the simulation")
parser.add_argument("--surpress", action='store_true', help="Only generate temp_params.h.  Do not call compile")
args = parser.parse_args()

zip_prefix="Dendrimer_5g"
batches = 5#25
batch_lines = 100#0
eps_mult = 1#  The number of duplicate epsilon values.  Useful to split batches into smaller batches

start_eps = 0.0
mid_eps = 0.2  #Seperation between high density region and low density region
end_eps = 2.0

numHD = 20
numLD = 20
dx = (mid_eps - start_eps)/numHD
hd_ax_t = np.linspace(start_eps,mid_eps - dx, numHD - 1)
ld_ax_t = np.linspace(mid_eps,end_eps,numLD + 1)

hd_ax_p = []
ld_ax_p = []

for i in hd_ax_t:
	for b in range(eps_mult):
		hd_ax_p.append(i)

for i in ld_ax_t:
	for b in range(eps_mult):
		ld_ax_p.append(i)

hd_ax = np.array(hd_ax_p)
ld_ax = np.array(ld_ax_p)

eps_ax = np.concatenate( (hd_ax, ld_ax) )  #generate epsilons for runs

seeds = np.empty(0)
while np.size(np.unique(seeds)) != (numLD + numHD)*eps_mult:  #verify all seeds are unique
	seeds = np.array([np.fromfile("/dev/urandom",dtype=np.uint64, count=1)[0] for i in range(eps_mult*(numLD+numHD)) ] )#generate random seeds

print eps_ax
print len(eps_ax)

for i,e in enumerate(eps_ax):
	fp = open("temp_params.h","w")
	fp.write("#define INIT_SEED "+str(seeds[i])+"\n")
	fp.write("#define ZIPNAME \""+zip_prefix+"_eps"+str(e)+"_"+str(seeds[i])+"\"\n" ) 
	fp.write("#define BATCHES "+str(batches)+"\n")
	fp.write("#define LINES "+str(batch_lines)+"\n")
	fp.write("#define MOL_EPSILON "+str(e)+"\n")
	fp.write("#define TOPOLOGY_FILE topo.csv \n")
	if hasattr(args,"func"):
		fp.write("#define FUNCTIONALITY "+str(args.func)+"\n" )
	fp.close()
	if not hasattr(args,"surpress"):
		os.system("./compile.sh")
		os.system("cp build/Point build/Point_"+str(i))




