def get_halfrad(lograds, loglums):
	from scipy import interpolate
	import math
	import numpy as np
	print('from halflight_math')
	N=np.ndim(lograds)
	if N==2:
		for n in range(0, len(lograds)):
			logL=loglums[n]
			logr=lograds[n]
	
			maxL=10**np.max(logL)
			halfL=maxL/2
			logL12=np.log10(halfL)
	
			f=interpolate.interp1d(logL,logr, kind='linear', axis=-1)
	
			logr12=f(loghalfL)
	return logr12
	
def my_linregress3(x,y,err):	#This one works!
	import numpy as np
	import math
	import scipy.stats as stats
	import scipy.optimize as opt
	#print('My linear regression!')
	N=len(x)
	Sigma = err**2*np.eye(N)
	Ah = np.vstack([x, np.ones(len(x))]).T
	Ae=np.dot(Ah.T, np.dot(np.linalg.inv(Sigma), Ah))
	be=np.dot(Ah.T, np.dot(np.linalg.inv(Sigma), y))### or y.T
	huh=np.linalg.lstsq(Ae, be)[0]
	
	slope=huh[0]
	inter=huh[1]
	
	S = np.linalg.inv(Ae)
	
	slope_err=np.sqrt(S[0,0])
	
	return slope, inter, slope_err