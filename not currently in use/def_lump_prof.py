#makes single plots
def lump_prof(newdata, bands, aperture, outdir, label, tag):
	import numpy as np
	from def_get_mags import get_zdistmod, get_kcorrect2, aper_and_comov, abs2lum, lumdensity, abs_mag
	from defclump import meanlum
	import matplotlib.pyplot as plt
	
	Naps=len(aperture)
	Ndat=len(newdata)
	redshifts=newdata['Z']
	
	DM= get_zdistmod(newdata, 'Z')

	kcorrect=get_kcorrect2(newdata,'mag_forced_cmodel', '_err', bands, '','hsc_filters.dat',redshifts)

	fig=plt.figure()
	bigLI=[]
	bigrad=[]
	for n in range(0, Ndat):
		#this goes through every galaxy
		LG=[]
		LR=[]
		LI=[]
		LZ=[]
		LY=[]
		string=str(n)
		grays=str(.999-n*0.00015)
		print('gray shade= ', grays)
		radkpc=aper_and_comov(aperture, redshifts[n])
		print('redshifts is ', redshifts[n])
	
		for a in range(0, Naps):	#this goes through every aperture
			ns=str(a)
			print('aperture0',ns)
		
			absg, absr, absi, absz, absy= abs_mag(newdata[n], 'mag_aperture0', kcorrect, DM[n], bands, ns, n) 
	
			Lumg, Lumr, Lumi, Lumz, Lumy=abs2lum(absg, absr, absi, absz, absy)
		
			Lg, Lr, Li, Lz, Ly=lumdensity(Lumg, Lumr, Lumi, Lumz, Lumy, radkpc[a])
		
			LG.append(Lg)
			LR.append(Lr)
			LI.append(Li)
			LZ.append(Lz)
			LY.append(Ly)
		print('LI for ',n,' galaxy is ', LI)
		#for cmodel
		absg, absr, absi, absz, absy= abs_mag(newdata[n], 'mag_forced_cmodel', kcorrect, DM[n], bands, '', n)
		
		bigLI.append(LI)
		bigrad.append(radkpc)
		plt.plot(radkpc, LI, marker='.', color=grays, alpha=0.9, zorder=1)
		#print(np.ndim(radkpc), np.ndim(LI), len(LI), len(radkpc))
	#print('Mins and Maxes: ', min(bigrad), max(bigrad))
	bigLI=np.array(bigLI)
	bigrad=np.array(bigrad)
	bigLI.flatten()
	bigrad.flatten()
	mean, error, radavg=meanlum(bigLI, bigrad,Naps, outdir='',scale='log', error='stdv')

	plt.plot(radavg, mean, marker='.', color='r',zorder=2)
	plt.errorbar(radavg, mean, yerr=error, color='m', fmt='.', zorder=3)
	plt.xlabel('Comoving Distance (kpc)', fontsize=10)
	plt.xscale('log')
	plt.yscale('log')
	plt.ylabel('Luminosity Density (Lsolar/kpc^2)', fontsize=10)
	plt.suptitle('Luminosity Density vs. Comoving Distance in band I', fontsize=15)
	plt.plot(0,0,label='Number of '+label+' = '+str(Ndat), c='k', marker='')
	plt.plot(0,0,label='Standard Deviations: ', c='m')
	for b in range(len(error)):
		bs=str(b)
		errors=round(error[b],4)
		errstr=str(errors)
		plt.plot(0,0,label='Error'+bs+'= '+errstr, c='m', marker='')

	plt.legend(loc=1,prop={'size':5.3})
	
	fig.savefig(outdir+'i'+tag+'_lumdens.pdf')
	print(outdir+'i'+tag+'_lumdens.pdf')
	
#returns values so we can overplot plots
def lump_prof2(newdata, bands, aperture, err=''):
	import numpy as np
	from def_get_mags import get_zdistmod, get_kcorrect2, aper_and_comov, abs2lum, lumdensity, abs_mag
	import math
	from defclump import meanlum2, get_error
	import matplotlib.pyplot as plt
	
	Naps=len(aperture)
	Ndat=len(newdata)
	redshifts=newdata['Z']
	
	DM= get_zdistmod(newdata, 'Z')

	kcorrect=get_kcorrect2(newdata,'mag_forced_cmodel', '_err', bands, '','hsc_filters.dat',redshifts)

	fig=plt.figure()
	bigLI=[]
	bigrad=[]
	for n in range(0, Ndat):
		#this goes through every galaxy
		#LG=[]
		#LR=[]
		LI=[]
		#LZ=[]
		#LY=[]
		string=str(n)
		radkpc=aper_and_comov(aperture, redshifts[n])
		print('redshifts is ', redshifts[n])
	
		for a in range(0, Naps):	#this goes through every aperture
			ns=str(a)
			print('aperture0',ns)
		
			absg, absr, absi, absz, absy= abs_mag(newdata[n], 'mag_aperture0', kcorrect, DM[n], bands, ns, n) 
	
			Lumg, Lumr, Lumi, Lumz, Lumy=abs2lum(absg, absr, absi, absz, absy)
		
			Lg, Lr, Li, Lz, Ly=lumdensity(Lumg, Lumr, Lumi, Lumz, Lumy, radkpc[a])
			
		
			Logli=math.log10(Li)
			#LG.append(Lg)
			#LR.append(Lr)
			LI.append(Logli)
			#LZ.append(Lz)
			#LY.append(Ly)
		print('LI for ',n,' galaxy is ', LI)
		#for cmodel
		absg, absr, absi, absz, absy= abs_mag(newdata[n], 'mag_forced_cmodel', kcorrect, DM[n], bands, '', n)
		
		bigLI.append(LI)
		lograd=[math.log10(radkpc[n]) for n in range(len(radkpc))]
		bigrad.append(lograd)
		#plt.plot(radkpc, LI, marker='.', color=grays, alpha=0.9)
		
	bigLIs=np.array(bigLI)
	bigrads=np.array(bigrad)
	
	bigLIs.flatten()
	bigrads.flatten()
	
	
	mean, radavg, bb=meanlum2(bigLIs, bigrads,Naps,scale='log')
	
	error=get_error(bigLIs, bigrads, bb, error=err)
	
	return bigLIs, bigrad, mean, error, radavg
	
def get_halflight(newdata, bands, aperture, scale=''):
	import numpy as np
	from def_get_mags import get_zdistmod, get_kcorrect2, aper_and_comov, abs2lum, lumdensity, abs_mag
	import math
	from defclump import meanlum2
	from my_def_plots import halflight_plot, scatter_fit
	from scipy import interpolate
	import matplotlib.pyplot as plt
	from def_mymath import halflight
	
	Naps=len(aperture)
	Ndat=len(newdata)
	redshifts=newdata['Z']
	
	DM= get_zdistmod(newdata, 'Z')
	kcorrect=get_kcorrect2(newdata,'mag_forced_cmodel', '_err', bands, '','hsc_filters.dat',redshifts)
	fig=plt.figure()
	bigLI=[]
	bigrad=[]
	bigden=[]
	for n in range(0, Ndat):
		#this goes through every galaxy
		#LG=[]
		#LR=[]
		LI=[]
		LI2=[]
		lumdi=[]
		#LZ=[]
		#LY=[]
		string=str(n)
		radkpc=aper_and_comov(aperture, redshifts[n])
		print('redshifts is ', redshifts[n])
	
		for a in range(0, Naps):	#this goes through every aperture
			ns=str(a)
			print('aperture0',ns)
		
			absg, absr, absi, absz, absy= abs_mag(newdata[n], 'mag_aperture0', kcorrect, DM[n], bands, ns, n) 
	
			Lumg, Lumr, Lumi, Lumz, Lumy=abs2lum(absg, absr, absi, absz, absy)			
			
			Lg, Lr, Li, Lz, Ly=lumdensity(Lumg, Lumr, Lumi, Lumz, Lumy, radkpc[a])
			
			if scale== 'log':
				print('getting logs')
				logLumi=math.log10(Lumi)
				logLi=math.log10(Li)
				LI.append(logLumi)
				lumdi.append(logLi)
			else:
				LI.append(Lumi)
				lumdi.append(Li)
			#LZ.append(Lz)
			#LY.append(Ly)
		print('LI for ',n,' galaxy is ', LI)
		bigLI.append(LI)
		bigden.append(lumdi)
		if scale== 'log':
			lograd=[math.log10(radkpc[n]) for n in range(len(radkpc))]
			bigrad.append(lograd)
		else:
			bigrad.append(radkpc)
		
	
	bigLIs=np.array(bigLI)

	bigrads=np.array(bigrad)
	
	lumdensi=np.array(bigden)
	
	halfrad=[]
	halflum=[]
	for x in range(len(bigrads)):
		#print('Galaxy #', str(x))
		#print('Luminosity= ', bigLIs[x])
		#print('Radii= ', bigrads[x])
		#f=interpolate.interp1d(bigLIs[x],bigrads[x], kind='linear', axis=-1)
			#half=(np.max(bigLIs[x])+np.min(bigLIs[x]))/2.0
		half=math.log10(10**np.max(bigLIs[x])/2.0)
		
		#print('Half Luminosity= ', half)
		#print('Half radius= ',f(half))
		fhalf=halflight(bigrads[x],bigLIs[x])
		
		halflum.append(half)
		halfrad.append(fhalf)
		
	halfrads=np.array(halfrad)
	halflums=np.array(halflum)
	
	print(len(bigrads),len(bigrads[0]))
	
	return halfrads, bigrads, lumdensi
	
def halflight_slopes(hx, x, y, names, plots='', outdir=''):
	print('Calculate the slope of lumdensity(r) outside or r(L1/2)')
	# hx(738,1), x(738,10), y(738,10)
	import numpy as np
	import scipy.stats as stats
	from def_mymath import my_linregress3
	from my_def_plots import scatter_fit, simple_hist
	#create array of all points greater than half value
	
	#print(len(hx), len(x), len(y))
	#print(x[0], hx[0])

	#x=x[:,:-2]	#removes last two columns
	#y=y[:,:-2]
	xcut=[]
	ycut=[]
	
	for i in range(len(x)):
		xrow=x[i]
		yrow=y[i]
		N=len(xrow)
		hhx=hx[i]
		#print(hhx)
		mx=xrow[xrow>=hhx]
		my=yrow[xrow>=hhx]
		#print('mx', mx)
		#print('my', my)
		xcut.append(mx)
		ycut.append(my)
		
	#print(len(xcut), len(xcut[0]))
	slopes=[]
	intercepts=[]
	for n in range(len(xcut)):
		print(len(x[n]), len(xcut[n]))
		#gets slope of outside r(L/2) in Lumdens
		slope, int, r_value, p_value, std_err = stats.linregress(xcut[n],ycut[n])
		#print('m=', slope, 'intercept= ', int)
		print('Error is ', std_err)
		slopes.append(slope)
		intercepts.append(int)
		if plots=='yes':
			#y=slope*(xcut[n]-int)
			yfit=slope*xcut[n]+int
			
			#print(xcut[n],y)
			#print(names[n])
			#scatter_fit(x[n],y[n], xcut[n], yfit, std_err, xscale='', yscale='', title='Half Light Radii v Luminosity Densities', sub=str(names[n]), outdir=outdir)
		elif plots=='no':
			print('No plots made')
	#simple_hist(slopes, xtitle='Slope', title='Frequencies of Slopes')
	return slopes
	
def get_mean_halflight(newdata, bands, aperture, scale=''):
	import numpy as np
	from def_get_mags import get_zdistmod, get_kcorrect2, aper_and_comov, abs2lum, lumdensity, abs_mag
	import math
	from defclump import meanlum2, get_error
	import matplotlib.pyplot as plt
	from def_mymath import halflight
	
	Naps=len(aperture)
	Ndat=len(newdata)
	redshifts=newdata['Z']
	
	DM= get_zdistmod(newdata, 'Z')

	kcorrect=get_kcorrect2(newdata,'mag_forced_cmodel', '_err', bands, '','hsc_filters.dat',redshifts)

	bigLI=[]
	bigrad=[]
	bigden=[]
	for n in range(0, Ndat):
		#this goes through every galaxy

		LI=[]
		LI2=[]
		lumdi=[]

		string=str(n)
		radkpc=aper_and_comov(aperture, redshifts[n])
		#print('redshifts are ', redshifts[n])
	
		for a in range(0, Naps):	#this goes through every aperture
			ns=str(a)
			#print('aperture0',ns)
		
			absg, absr, absi, absz, absy= abs_mag(newdata[n], 'mag_aperture0', kcorrect, DM[n], bands, ns, n) 
	
			Lumg, Lumr, Lumi, Lumz, Lumy=abs2lum(absg, absr, absi, absz, absy)			
			
			Lg, Lr, Li, Lz, Ly=lumdensity(Lumg, Lumr, Lumi, Lumz, Lumy, radkpc[a])
			
			#LG.append(Lg)
			#LR.append(Lr)
			#LI.append(Logli)
			if scale== 'log':
				#print('getting logs')
				logLumi=math.log10(Lumi)
				logLi=math.log10(Li)
				LI.append(logLumi)
				lumdi.append(logLi)
			else:
				LI.append(Lumi)
				lumdi.append(Li)
			#LZ.append(Lz)
			#LY.append(Ly)
		#print('Luminosity(i) for ',n,' galaxy is ', LI)
		bigLI.append(LI)
		bigden.append(lumdi)
		if scale== 'log':
			lograd=[math.log10(radkpc[n]) for n in range(len(radkpc))]
			bigrad.append(lograd)
		else:
			bigrad.append(radkpc)
		
	bigLIs=np.array(bigLI)
	bigrads=np.array(bigrad)
	lumdensi=np.array(bigden)
	bigLIs.flatten()	#luminosity
	bigrads.flatten()	#radii
	lumdensi.flatten()	#luminosity density
	#print('Check: ', bigLIs[0], lumdensi[0])
	
	meanlum, radavg, bb=meanlum2(bigLIs, bigrads,Naps,scale='log')
	meandens, radavg, bb=meanlum2(lumdensi, bigrads,Naps,scale='log')
	
	err='bootstrap_stdv'
	lumdenerr=get_error(lumdensi, bigrads, bb, error=err)
	
	print('Mean Luminosity= ', meanlum)
	print('Mean LumDensity=', meandens)
	print('Mean Radii= ', radavg)
	print('Standard Deviation= ', lumdenerr)
	
	halfrad=halflight(radavg,meanlum)
	
	print('Half Radius is ',halfrad)
	
	return meandens, radavg, lumdenerr, halfrad
	
def halflight_meanslope(hx, x, y, error):
	print('Calculate the slope of lumdensity(r) outside or r(L1/2)')
	# hx(738,1), x(738,10), y(738,10)
	import numpy as np
	import scipy.stats as stats
	from def_mymath import my_linregress3
	x=np.array(x)
	xcut=x[x>=hx]
	ycut=y[x>=hx]
	errcut=error[x>=hx]
	
	print('Radii are= ', x)
	print('R1/2 is= ', hx)
	
	#weights=1/(errcut**2)
	#p, r=np.polynomial.polynomial.polyfit(xcut, ycut, 1, w = weights, full=True)
	#slope, int, r_value, p_value, std_err = stats.linregress(xcut,ycut)
	
	sl3, C3, std_err3=my_linregress3(xcut, ycut, errcut)
	
	#c=p[0]
	#m=p[1]
	#print('Residuals: ', r)
	#print('From polyfit: ', m, c)
	
	#print('From Linregress: ', slope, int, 'Stdv= ', std_err)

	print('From my formula3: ', sl3, C3, 'Stdv= ', std_err3)
	
	
	#y1 = m * xcut + c
	#y2 = slope *xcut +int
	
	y4 = sl3 * xcut + C3
	
	test=''
	if test=='yes':
		import matplotlib.pyplot as plt
		figs=plt.figure()
		plt.scatter(x, y, color='k', marker='o')
		plt.errorbar(x, y, yerr=error, fmt='.',color='k', label=error)
		
		#plt.plot(xcut, y1, color='r', label='Polyfit (weighted)')
		#plt.plot(xcut, y2, color='b',linestyle='--', label='Scipy Linregress (unweighted)')
		plt.plot(xcut, y4, color='g', label='My Weighted Linear Regression')
		
		plt.xlabel('Log Radii (kpc)')
		plt.ylabel('Luminosity Densities (Lsolar/kpc^2)')
		plt.legend(loc=0,prop={'size':7.0})
		plt.show()
	
	return sl3, C3, xcut, ycut, std_err3, errcut