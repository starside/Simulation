import matplotlib.pyplot as plt
import numpy as np

fh = open("rg2.csv","rb")
while fh:
	v0 = np.fromfile(fh,dtype=np.float64,count=1)
	size = np.fromfile(fh,dtype=np.int32,count=1)
	hist_data = np.fromfile(fh,dtype=np.int32, count=size)
	label_data = np.fromfile(fh,dtype=np.float64, count=size)
	plt.plot(label_data, hist_data)
	plt.show()
fh.close()