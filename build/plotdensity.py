import matplotlib.pyplot as plt
import numpy as np

fh = open("density2d_17190222644612017813.csv","rb")

for i in range(40):
	v0 = size = np.fromfile(fh,dtype=np.float64,count=1)
	print (v0,i)
	size = np.fromfile(fh,dtype=np.int32,count=1)
	monomers = np.fromfile(fh, dtype=np.int32,count=1)
	dim = int(np.sqrt(size))
	fullImage = None
	for j in range(monomers):
		image_data = np.fromfile(fh,dtype=np.float64, count=size)
		image = image_data.reshape((dim,dim))
		if fullImage is None :
			fullImage = 0*image
		else:
			fullImage = fullImage + image
		#print j
		#plt.imshow(image)
		#plt.show()
	plt.imshow(fullImage)
	plt.show()
fh.close()