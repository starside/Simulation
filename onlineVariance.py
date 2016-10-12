import numpy as np

N = 5000
M = 5000

def olv(data, nA=0, meanA=0.0, M2A=0.0):
	n = nA
	mean = meanA
	M2 = M2A

	for x in data:
		n += 1
		delta = x - mean
		mean += delta/n
		M2 += delta*(x-mean)

	if n < 2:
		return double('nan')
	else:
		return M2/ (n-1), n, mean, M2

a1 = np.random.normal(size=(N))
a2 = np.random.normal(size=(M))

ac =  np.concatenate((a1,a2) )
mstd1, n1, m1, m21 = olv(a1)
mstd2, n2, m2, m22 = olv(a2)
print olv(ac)

#Number of points, mean, and sum squared x2 
def OLVCombine(n1,m1,m21, n2,m2,m22):
	delta = m2 - m1
	xb = (n1*m1 + n2*m2)/(n1+n2)
	M2X = m21 + m22 + delta*delta*(n1*n2/(n1+n2))
	return M2X/(n1+n2 - 1)