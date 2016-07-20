import numpy as np
from def_get_mags import get_zdistmod, get_kcorrect2, aper_and_comov, abs2lum, lumdensity, abs_mag

import math
from halflight_second import meanlum2, get_errors
from def_halflight_math import get_halfrad


def upper_rad_cut(loglum, lograd, logden, m, proof=False): #this should get rid of galaxies outside 4r1/2
	from def_halflight_math import get_halfrad
	nloglum=[]
	nlograd=[]
	nlogden=[]
	mult=m
	print(len(loglum), len(lograd))
	N=len(lograd)
	for n in range(0,N):
		loglums=loglum[n]
		lograds=lograd[n]
		logdens=logden[n]
		logr12=get_halfrad(lograds,loglums)
		r12=10**logr12
		r412=mult*r12
		logr412=np.log10(r412)
		if proof == True:
			print('logr1/2= ', logr12)
			print('log4r1/2= ', logr412)
			print('The radii are ', lograds)
		print(logr412)
		if np.max(lograds) >= logr412:
			logradi=lograds[(lograds>=logr412)&(lograds<=logr412)]
			print(logradi)
			if len(logradi)>=3:
				nloglum.append(loglums)
				nlograd.append(lograds)
				nlogden.append(logdens)
			else:
				print('not enough data points')
		else:
			print('Upper limit out of range')
	nloglum=np.array(nloglum)
	nlograd=np.array(nlograd)
	nlogden=np.array(nlogden)
	return nloglum, nlograd, nlogden

def get_ind_lums(newdata, bands, aperture, scale=''):
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
	try:
		redshifts=newdata['Z']
		DM= get_zdistmod(newdata, 'Z')
	except:
		redshifts=newdata['Z_2']
		DM= get_zdistmod(newdata, 'Z_2')
	kcorrect=get_kcorrect2(newdata,'mag_forced_cmodel', '_err', bands, '','hsc_filters.dat',redshifts)
	fig=plt.figure()
	bigLI=[]
	bigrad=[]
	bigden=[]
	for n in range(0, Ndat):
		LI=[]
		LI2=[]
		lumdi=[]
		string=str(n)
		radkpc=aper_and_comov(aperture, redshifts[n])
		#print('redshifts is ', redshifts[n])
		for a in range(0, Naps):	#this goes through every aperture
			ns=str(a)
			#print('aperture0',ns)
			absg, absr, absi, absz, absy= abs_mag(newdata[n], 'mag_aperture0', kcorrect, DM[n], bands, ns, n) 
			Lumg, Lumr, Lumi, Lumz, Lumy=abs2lum(absg, absr, absi, absz, absy)			
			Lg, Lr, Li, Lz, Ly=lumdensity(Lumg, Lumr, Lumi, Lumz, Lumy, radkpc[a])
			if scale== 'log':
				#print('getting logs')
				logLumi=math.log10(Lumi)
				logLi=math.log10(Li)
				LI.append(logLumi)
				lumdi.append(logLi)
			else:
				LI.append(Lumi)
				lumdi.append(Li)
		#print('LI for ',n,' galaxy is ', LI)
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
	return bigLIs, bigrads, lumdensi
	
def get_avg_lums(logLs, lograds, logLDs, type='', scale=''):
	print('get_avg_lums is in halflight first')
	sc=scale
	#sc is whether or nto we stack the linear data or log data
	Naps=0.0
	
	if type=='mean':
		meanlum, radavg, bb=meanlum2(logLs, lograds, Naps,scale=sc)
		meandens, radavg, bb=meanlum2(logLDs, lograds,Naps,scale=sc)
	
		err='bootstrap_stdv'
		lumdenerr=get_errors(logLDs, lograds, bb, error=err, scale=sc)
	
		print('Mean Luminosity= ', meanlum)
		print('Mean LumDensity=', meandens)
		print('Binned Radii= ', radavg)
		print('Standard Deviation= ', lumdenerr)
	
		return meanlum, meandens, radavg, lumdenerr
		
	if type== 'median':
		medlum, radavg, bb=medlum2(bigLIs, bigrads)
		medens, radavg, bb=medlum2(lumdensi, bigrads)
		err='bootstrap_stdv'
		lumdenerr=get_error(lumdensi, bigrads, bb, error=err)
	
		print('Median Luminosity= ', medlum)
		print('Median LumDensity=', medens)
		print('Binned Radii= ', radavg)
		print('Standard Deviation= ', lumdenerr)
		
		return medlum, medens, radavg, lumdenerr
		
def get_halflight(logLs, lograds):
	from scipy import interpolate
	import math
	import numpy as np
	print('not from halflight_math')
	N=np.ndim(lograds)
	if N==2:
		logr12=[]
		for n in range(0, len(lograds)):
			logL=logLs[n]
			logr=lograds[n]
			maxL=10**np.max(logL)
			halfL=maxL/2
			logL12=np.log10(halfL)
			f=interpolate.interp1d(logL,logr, kind='linear', axis=-1)
			alogr12=f(logL12)
			logr12.append(alogr12)
		logr12=np.array(logr12)
	else:
		logL=logLs
		logr=lograds
		maxL=10**np.max(logL)
		halfL=maxL/2
		logL12=np.log10(halfL)
		f=interpolate.interp1d(logL,logr, kind='linear', axis=-1)
		logr12=f(logL12)
		
	return logr12
		
def get_slopes(logr12s, lograd, logld, error=None, smax=False):
	import scipy.stats as stats
	from def_halflight_math import my_linregress3
	from my_def_plots import scatter_fit, simple_hist
	print('slopes from halflight_first')
	mult=4
	Ndim=np.ndim(lograd)
	N=len(lograd)
	if error is None:
		print('No error was given')
		error=np.ones((N, len(lograd[0])))

	N=np.ndim(lograd)
	
	logrcut=[]
	logldcut=[]
	errcut=[]
	if N==2:
		for i in range(len(lograd)):
			logrrow=lograd[i]
			logldrow=logld[i]
			errow=error[i]
			logr12=logr12s[i]
			
			if smax== True:
				r12=10**logr12
				r412=mult*r12
				logr412=np.log10(r412)
				
				#print(hhx2, np.max(xrow))
				if np.max(logrrow) >= logr412:
					mlogr=logrrow[(logrrow>=logr12)&(logrrow<=logr412)]
					mlogld=logldrow[(logrrow>=logr12)&(logrrow<=logr412)]
					merr=errow[(logrrow>=logr12)&(logrrow<=logr412)]
					if len(mlogr) >=4:
						logrcut.append(mlogr)
						logldcut.append(mlogld)
						errcut.append(merr)
				else:
					print('Upper Cut is Out of the Radius Range')
			else:
				merr=errow[logrrow>=logr12]
				mlogr=logrrow[logrrow>=logr12]
				mlogld=logldrow[logrrow>=logr12]
				if len(mlogr) >=4:
					logrcut.append(mlogr)
					logldcut.append(mlogld)
					errcut.append(merr)
		slopes=[]
		intercepts=[]
		errs=[]
		for n in range(len(logrcut)):
			slope, int, std_err=my_linregress3(logrcut[n], logldcut[n], errcut[n])
			slopes.append(slope)
			intercepts.append(int)
			errs.append(std_err)
		return slopes, intercepts, errs
	else: #for arrays of 1D *aka* the stacked profile
		lograd=np.array(lograd)
		print('r1/2 limit is ', logr12s)
		print('xrange for stacked is ', lograd)
		if error is None:
			error=np.ones(N)
		if smax== True:
			r12=10**logr12
			r412=mult*r12
			logr412=np.log10(r412)
			print('upper limit is ', logr412)
			if np.max(lograd) <= logr412:
				print('Upper cut is out of the Radius range')
			else:
				logrcut=lograd[(lograd>=logr12)&(lograd<=logr412)]
				logldcut=logld[(lograd>=logr12)&(lograd<=logr412)]
				errcut=error[(lograd>=logr12)&(lograd<=logr412)]
		else:
			logrcut=lograd[lograd>=logr12s]
			logldcut=logld[lograd>=logr12s]
			errcut=error[lograd>=logr12s]
		print('Log Radii are= ', lograd)
		print('LogR1/2 is= ', logr12s)
		
		sl3, C3, std_err3=my_linregress3(logrcut, logldcut, errcut)
		
		return sl3, C3, logrcut, logldcut, std_err3, errcut