import matplotlib.pyplot as plt
import numpy as np

fh = open("density2d.csv","rb")

for i in range(40):
	print i
	size = np.fromfile(fh,dtype=np.int32,count=1)
	dim = int(np.sqrt(size))
	image_data = np.fromfile(fh,dtype=np.float64, count=size)
	image = image_data.reshape((dim,dim))
	plt.imshow(image)
	plt.show()
fh.close()