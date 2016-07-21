from scipy import interpolate
import math
import numpy as np
def get_halfrad(R, L):
	print('from halflight_math')
	#should be in linear
	maxL=np.max(L)
	halfL=maxL/2
			
	f=interpolate.interp1d(L,R, kind='linear', axis=-1)
	
	r12=f(halfL)
	return r12
	
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