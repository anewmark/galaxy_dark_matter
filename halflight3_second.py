def meanlum2(lumden, rads, Naps, grange=[]):
	import numpy as np
	import math
	import matplotlib.pyplot as plt
	import matplotlib.mlab as mlab
	import scipy.stats as stats
	#import as 1D clumps
	print('this is from half-light3')

	print('Regular Rads Check of one row: ', rads[0])
	
	print('min= ', np.min(rads), 'max= ', np.max(rads))
	
	print('regular lumdens= ', lumden[0])
	
	try:
		radmin=grange[0]
	except:
		radmin=1 ###### <======== change this one occasionally
	try:
		radmax=grange[1]
	except:
		radmax=80
	try:
		ns=grange[2]
	except:
		ns=11
		
	bb=np.logspace(math.log10(radmin), math.log10(radmax),num=ns, endpoint=False)
	
	print('bb= ',bb)
	Naps=len(bb)

	hist1, radhists=np.histogram(rads, bins=bb, weights=lumden, density=False)	#hist1=average value in each bin
	hist2, radhists=np.histogram(rads, bins=bb, density=False)	#hist2=number of galaxies in each bin
	#radhists here are in 10^n, not log scale

	means=hist1/hist2

	print('summed lumdens:', hist1)
	print('num in each bin:', hist2)
	print('radhists (hopefully in linear): ', radhists)
	
	#Is this the way to calculate the bin_centers?
	bin_centers = 2.*(radhists[1:]**3 - radhists[:-1]**3)/(3.*(radhists[1:]**2 - radhists[:-1]**2))

	print('Mean= ', means)
	print('bincenters= ',bin_centers)
	
	test=''
		#this will test the stacked profile against the individual galaxies
	if test=='yes':
		print(len(lumden))
		print('log of mean= ', np.log10(means))
		print('log of bincenters= ', np.log10(bin_centers))
		for n in range(len(lumden)):
			plt.plot(np.log10(rads[n]), np.log10(lumden[n]), color='lightgrey', marker='.', zorder=1)
		plt.scatter(np.log10(bin_centers), np.log10(means), zorder=2, color='red', marker='o')
		plt.xlabel('Log Radii')
		plt.ylabel('Log Luminosity/Density')
		plt.title('Mean')
		plt.show()
		
	return means, bin_centers, bb
	
def get_errors(lumden, rads, bb, meanL, error=''):
	import numpy as np
	import math
	import matplotlib.pyplot as plt
	import matplotlib.mlab as mlab
	
	if error=='bootstrap_stdv':
		print('This uses the conept of replacement in sampling')
		print('CHECK: first row of rads should be in linear space= ', rads[0])
		print('CHECK: luminosity density is in linear space= ', lumden[0])

		Means=[]
		for m in range(500):
			#EVERYTHING IS ALREADY IN LINEAR
			#print('lumden.shape[0]= ', lumden.shape[0])
			#print('lumden.shape[1]= ', lumden.shape[1])
			#print('len(lumden)= ', len(lumden))
			
			newarr=lumden[np.random.choice(lumden.shape[0], size=len(lumden), replace=True),:]
			#print('shape of newarr: ', np.shape(newarr), 'length of one row= ', len(newarr[0]))
			h1, r1=np.histogram(rads, bins=bb, weights=newarr, density=False)
			h2, r2=np.histogram(rads, bins=bb, density=False)
			littlemeans=h1/h2
			Means.append(littlemeans)
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
		sigmabs=np.std(Means, axis=0)
		
		logerror=sigmabs/meanL

		return  logerror

	else:
		print('hi')

def medlum2(lumden, rads, Naps, grange=[]):
	import numpy as np
	import math
	import scipy.stats as stats
	import matplotlib.pyplot as plt
	Naps=0.0
	#since in log scale
	#rads in 10logspace
	print('Getting median!')
	print('Regular Rads Check of one row: ', rads[0])
	
	print('min= ', np.min(rads), 'max= ', np.max(rads))
	
	print('regular lumdens= ', lumden[0])
	
	try:
		radmin=grange[0]
	except:
		radmin=1 ###### <======== change this one occasionally
	try:
		radmax=grange[1]
	except:
		radmax=80
	try:
		ns=grange[2]
	except:
		ns=10
	
	bb=np.logspace(math.log10(radmin), math.log10(radmax),num=ns, endpoint=False)
	rads=rads.flatten()
	lumden=lumden.flatten()
	
	print(np.ndim(rads), np.ndim(lumden))
	
	medarr, radhists, ind=stats.binned_statistic(rads, lumden, statistic='median', bins=bb)
		
	bin_centers = 2.*(radhists[1:]**3 - radhists[:-1]**3)/(3.*(radhists[1:]**2 - radhists[:-1]**2))
	
	bintest=''
	if bintest:
		print(len(lumden))
		print('log of median= ', np.log10(medarr))
		print('log of bincenters= ', np.log10(bin_centers))

		for n in range(len(lumden)):
			plt.plot(np.log10(rads[n]), np.log10(lumden[n]), color='lightgrey', marker='.', zorder=1)
		plt.scatter(np.log10(bin_centers), np.log10(medarr), zorder=2, color='red', marker='o')
		plt.xlabel('Log Radii')
		plt.ylabel('Log Luminosity/Density')
		plt.title('Median')
		plt.show()
	
	return medarr, bin_centers, bb