def meanlum2(logdatarr, lograds, Naps, scale=''):
	import numpy as np
	import math
	import matplotlib.pyplot as plt
	import matplotlib.mlab as mlab
	import scipy.stats as stats
	#import as 1D clumps
	print('this is from half-light second')
	print(lograds)
	rads=10**lograds
	print('min= ', np.min(rads), 'max= ', np.max(rads))
	
	if scale=='lindata':
		datarr=10**logdatarr
	else:
		datarr=logdatarr
	
	radmin=1 ###### <======== change this one occasionally
	radmax=80
	bb=np.logspace(math.log10(radmin), math.log10(radmax),num=11, endpoint=True)
	print('bb= ',bb)
	Naps=len(bb)

	hist1, radhists=np.histogram(rads, bins=bb, weights=datarr, density=False)	#hist1=average value in each bin
	hist2, radhists=np.histogram(rads, bins=bb, density=False)	#hist2=number of galaxies in each bin
	#radhists here are in 10^n, not log scale

	means=hist1/hist2
	ind=np.digitize(rads, bins=bb) #says which bin each galaxy is in

	print('hist1', hist1)
	print('hist2', hist2)
	
	
	bin_centers = 2.*(radhists[1:]**3 - radhists[:-1]**3)/(3.*(radhists[1:]**2 - radhists[:-1]**2))
		
	if scale=='lindata':
		logmeans=means
		means=np.log10(np.array(logmeans))
	logbin_centers=np.log10(np.array(bin_centers))
	print('Mean= ', means)
	print('log bincenters= ',logbin_centers)
	return means, logbin_centers, bb
	
def get_errors(logdatarr, lograds, bb, error='', scale=''):
	import numpy as np
	import math
	import matplotlib.pyplot as plt
	import matplotlib.mlab as mlab
	
	if error=='bootstrap_stdv':
		print('This uses the conept of replacement in sampling')
		rads=10**lograds

		Means=[]
		for m in range(500):
			if scale=='lindata':
				datarr=10**logdatarr
				newarr=datarr[np.random.choice(datarr.shape[0], size=len(datarr), replace=True),:]
			else:	
				newarr=logdatarr[np.random.choice(logdatarr.shape[0], size=len(logdatarr), replace=True),:]
		#	print(np.shape(newarr), len(newarr[0]))
			h1, r1=np.histogram(rads, bins=bb, weights=newarr, density=False)
			h2, r2=np.histogram(rads, bins=bb, density=False)
			means=h1/h2
			Means.append(means)
		Means=np.array(Means)
		#print(Means)
		test=''
		if test:
			mymean=np.mean(Means[:,1])
			mystd=np.std(Means[:,1])
			print(mymean, mystd)
			plt.hist(Means[:,1], bins=10, alpha=0.7)
			plt.axvline(mymean, ymin=0, ymax=1, c='r')
			plt.xlabel('Mean Luminosities')
			plt.axvline(mymean+mystd, ymin=0, ymax=1, c='g')
			plt.axvline(mymean-mystd, ymin=0, ymax=1, c='g')
			plt.show()
			
		error=np.std(Means, axis=0)
		if scale=='lindata':
			error=np.std(np.log10(Means), axis=0)
		return  error

	else:
		print('hi')
