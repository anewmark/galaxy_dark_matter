def TF_meanSB(dataflag, datanot, aper, band, magname, labs,zrange, outdir):
	import matplotlib.pyplot as plt
	import numpy as np
	import math
	import matplotlib.patches as mpatches
	pi=math.pi
	
	if not zrange:
		zrange='None Specified'
	else:
		zrange=str(zrange)
	nflag=len(dataflag)
	numflag=str(nflag) #string version
	nnot=len(datanot)
	numnot=str(nnot) #string version
	clumpflag=[]
	clumpnot=[]
	Naps=len(aper)
	
	if magname is 'mag_cmodel':
		tag='cmodel'
		print(tag)
	elif magname is'mag_aperture0':
		tag='aper'
		print(tag)
	else:
		tag=''
		print('no tag')
	
	#print('mag ',datanot[band+magname+'0'])
	e=math.exp(1)
	errflag=[]
	errnot=[]
	for a in range(0, Naps):
		js=str(a)
		mag=dataflag[band+magname+js]
		sb=mag+2.5*math.log10(4*pi*aper[a]**2)
		Fflag=10**(-0.4*sb) #getting 10^SB
		#	Fflag.append(F)
		Fflagmean=np.mean(Fflag)
		sbavg=-2.5*math.log10(Fflagmean)
		#print(sbavg)
		clumpflag.append(sbavg)
	#*********
		#error
		stdFflag=np.std(Fflag)
		dmflag=abs(-2.5*math.log10(e)*stdFflag/Fflagmean)
		errflag.append(dmflag)
	#*********

		if nnot==0:
			clumpnot=np.zeros(Naps)
		else:
			magn=datanot[band+magname+js]
			sbn=magn+2.5*math.log10(4*pi*aper[a]**2)
			Fnot=10**(-0.4*sbn) #getting 10^SB
			Fnotmean=np.mean(Fnot)
			sbnotavg=-2.5*math.log10(Fnotmean)
			clumpnot.append(sbnotavg)
		#********
			stdFnot=np.std(Fnot)
			#print(stdFnot, Fnotmean)
			dmnot=abs(-2.5*math.log10(e)*stdFnot/Fnotmean)
			errnot.append(dmnot)
		#*********
	print('no flag error ',errflag, 'flagged error ', errnot)

	fig=plt.figure()
	plt.plot(aper, clumpflag, c='r', marker='^')
	plt.errorbar(aper,clumpflag,yerr=errflag, ecolor='m', fmt='none')
	if nnot==0:
		print('No Flagged Galaxies Here!')
	else:
		plt.plot(aper, clumpnot, c='b', marker='*')
		plt.errorbar(aper, clumpnot, yerr=errnot, ecolor='c', fmt='none')
	plt.xlabel('Aperture Radius Log Scale (arcseconds)', fontsize=10)
	plt.xscale('log')
	plt.ylabel('Mean Surface Brightness (mag/arcsec^2)', fontsize=10)
	plt.suptitle("Mean Surface Brightness vs. Aperture Radius in "+band, fontsize=15)
	plt.title("Redshift range "+ zrange, fontsize=12)
	red_patch = mpatches.Patch(color='red', label=labs[1] +'('+numflag+')')
	blue_patch = mpatches.Patch(color='blue', label=labs[2] +'('+numnot+')')
	magenta_patch = mpatches.Patch(color='magenta', label=labs[1] +' error')
	cyan_patch = mpatches.Patch(color='cyan', label=labs[2] +' error')
	plt.legend(handles=[red_patch, magenta_patch, blue_patch, cyan_patch], loc=7,prop={'size':5})
	plt.gca().invert_yaxis()
	#plt.show()
	fig.savefig(outdir+band+'_'+labs[3]+'+err_TFmeanSB.pdf')
	print(outdir+band+'_'+labs[3]+'_TFmeanSB.pdf')
	
	two_plots='hi'
	if two_plots:
		f=plt.figure()
		if not zrange:
				zrange='None Specified'
		else:
			zrange=str(zrange)
		print(labs[3])
		f, (ax0, ax1) = plt.subplots(1,2, sharey=True)
		f.suptitle("Mean Surface Brightness vs. Aperture Radius in "+band)
		ax0=plt.subplot(121)
		ax0.plot(aper, clumpflag, c='r', marker='^')
		ax0.set_xscale('log')
		ax0.set_xlabel('Aperture Radius in Log Scale (arcseconds)', fontsize=6.5)
		ax0.errorbar(aper,clumpflag,yerr=errflag, ecolor='m', fmt='none')
		ax0.set_ylabel('Mean Surface Brightness (mag/arcsec^2)', fontsize=6.5)
		ax0.set_title('Number of Galaxies Not Flagged= '+numflag, fontsize=9)
		plt.ylim(min(clumpflag)-errflag[0],max(clumpflag)+errflag[9])
		plt.tick_params(axis='both', which='major', labelsize=5)
		ax0.invert_yaxis()
		ax0.plot(0,0, label='Redshift is '+zrange, marker='', c='k')
		#ax0.plot(0,0, label=labs[1], marker='', c='k')
		#ax0.plot(0,0,label='Number of Galaxies= '+numflag, marker='', c='k')
		ax0.legend(loc=7,prop={'size':5}, numpoints = 1)
	
		ax1=plt.subplot(122)
		if nnot==0:
			ax1.set_title('Number of Flagged Galaxies= '+numnot, fontsize=9)
			#ax1.plot(0,0,label='No Galaxies were Flagged', marker='', c='k')
			#ax1.legend(loc=9,prop={'size':8}, numpoints = 1)
			#plt.tick_params(axis='both', which='major', labelsize=5)
			ax1.invert_yaxis()
		else:
			ax1.plot(aper, clumpnot, c='b', marker='*')
			ax1.set_xlabel('Aperture Radius in Log Scale (arcseconds)', fontsize=6.5)
			ax1.set_xscale('log')
			ax1.set_ylabel('Mean Surface Brightness (mag/arcsec^2)', fontsize=6.5)
			ax1.set_title('Number of Flagged Galaxies= '+numnot, fontsize=9)
			ax1.errorbar(aper, clumpnot, yerr=errnot, ecolor='c', fmt='none')
			plt.ylim(min(clumpnot)-errnot[0],max(clumpnot)+errnot[9])
			plt.tick_params(axis='both', which='major', labelsize=5)
			ax1.invert_yaxis()
			ax1.plot(0,0, label='Redshift is '+zrange, marker='', c='k')
			#ax1.plot(0,0, label=labs[2], marker='', c='k')
			#ax1.plot(0,0,label='Number of Flagged Galaxies= '+numnot, marker='', c='k')
			ax1.legend(loc=7,prop={'size':5}, numpoints = 1)
	
		#plt.show()
		f.savefig(outdir+band+'_'+labs[3]+'_2TFmeanSB.pdf')
		print(outdir+band+'_'+labs[3]+'_2TFmeanSB.pdf')
	else:
		print('only making combined plots')
	
def meanlum2(datarr, rads, Naps, scale=''):
	import numpy as np
	import math
	import matplotlib.pyplot as plt
	import matplotlib.mlab as mlab
	import scipy.stats as stats
	#import as 1D clumps
	print('Length of galaxy array: ', len(datarr))

	if scale=='linear':
		print('This is a linear scale, forgot to do bin_centers, but didnt do anything for it yet')
	else:
		print('This is a log scale')
		print('min= ', rads.min(), 'max= ', rads.max())

		rads=np.array(rads)

		radmin=1 ###### <======== change this one occasionally
		radmax=80
		#bb=np.logspace(math.log10(radmin), math.log10(radmax),num=11, endpoint=True)
		bb=np.linspace(radmin, radmax, num=11, endpoint=True)
		print('bb= ',bb)
		Naps=len(bb)

		hist1, radhists=np.histogram(rads, bins=bb, weights=datarr, density=False)	#hist1=average value in each bin
		hist2, radhists=np.histogram(rads, bins=bb, density=False)	#hist2=number of galaxies in each bin
		#radhists here are in 10^n, not log scale

		hist=hist1/hist2
		ind=np.digitize(rads, bins=bb) #says which bin each galaxy is in
	
		print('hist1', hist1)
		print('hist2', hist2)
		print('Mean= ', hist)
		
		bin_centers = 2.*(radhists[1:]**3 - radhists[:-1]**3)/(3.*(radhists[1:]**2 - radhists[:-1]**2))
		
		#bin_centers=[math.log10(n) for n in bin_centers]
		
		print('radhists= ',bin_centers)
		return hist, bin_centers, bb
		
def medlum2(datarr, rads):
	import numpy as np
	import math
	import scipy.stats as stats
	#since in log scale
	trads=rads
	rads=[]
	for i in range(len(trads)):
		grad=trads[i]
		rad=[10**n for n in grad]
		rads.append(rad)
	#now have rads in 10^rad scale
	radmin=1
	radmax=80
	bb=np.logspace(math.log10(radmin), math.log10(radmax),num=11, endpoint=True)
	
	
	rads=np.array(rads)
	rads=rads.flatten()
	datarr=datarr.flatten()
	#datarr=np.array(datarr)
	#print(len(rads), len(datarr), len(rads[0,:]), len(datarr[0,:]))
	print(np.ndim(rads), np.ndim(datarr))
	medarr, radhists, ind=stats.binned_statistic(rads, datarr, statistic='median', bins=bb)
		
	bin_centers = 2.*(radhists[1:]**3 - radhists[:-1]**3)/(3.*(radhists[1:]**2 - radhists[:-1]**2))
	bin_centers=[math.log10(n) for n in bin_centers]

	
	
	return medarr, bin_centers, bb
		
def get_error(datarr, rads, bb, error=''):
	import numpy as np
	import math
	import matplotlib.pyplot as plt
	import matplotlib.mlab as mlab
	#rads in in log10, bb is in 10^n
	if error=='stdv':
		bin_stdv=[np.std(datarr[ind==i]) for i in range(1, len(bb))]
		print('std= ', bin_stdv)
		return hist, bin_stdv, bin_centers
	
	elif error=='bootstrap_stdv':
		print('This uses the conept of replacement in sampling')
		trads=rads
		rads=[]

		rads=np.array(rads)
		
		Means=[]
		for m in range(500):
			newarr=datarr[np.random.choice(datarr.shape[0], size=len(datarr), replace=True),:]
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
		return  error

	else:
		return hist, hist2, bin_centers