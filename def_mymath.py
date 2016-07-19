def halflight(x, y):
	from scipy import interpolate
	import math
	import numpy as np
	f=interpolate.interp1d(y,x, kind='linear', axis=-1)
	half=math.log10(10**np.max(y)/2.0)
	halfrad=f(half)
	return halfrad
	
def my_linregress2(X,Y,err):
	import numpy as np
	import math
	import scipy.stats as stats
	import scipy.optimize as opt
	
	slope, inter, r_value, p_value, std_err = stats.linregress(X,Y)
	
	print('Standard Error in slope before error fit= ', std_err)
	
	p=(slope, inter)
	
	fit = lambda p, x: p[1] + p[0] * x
	errfunc= lambda p,x,y, err,: (y-fit(p, x))/err
	
	out=opt.leastsq(errfunc,p,args=(X,Y,err), full_output=1)
	
	return out[0][0], out[0][1]
	
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