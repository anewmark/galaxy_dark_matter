import numpy as np
from def_get_mags import get_zdistmod, get_kcorrect2, aper_and_comov, abs2lum, lumdensity, abs_mag
from scipy import interpolate
import math
import matplotlib.pyplot as plt
from halflight3_second import meanlum2, get_errors, medlum2
from def_halflight_math import get_halfrad

def get_ind_lums(newdata, bands, aperture, my_band='i', scale=''):
	import numpy as np
	from def_get_mags import get_zdistmod, get_kcorrect2, aper_and_comov, abs2lum, lumdensity, abs_mag
	import math
	from defclump import meanlum2
	from my_def_plots import halflight_plot, scatter_fit
	from scipy import interpolate
	import matplotlib.pyplot as plt
	from def_mymath import halflight
	print('getting this in halflight3')
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
	bigL=[]
	bigrad=[]
	bigden=[]
	for n in range(0, Ndat):
		L=[]
		L2=[]
		lumd=[]
		string=str(n)
		radkpc=aper_and_comov(aperture, redshifts[n])
		#print('redshifts is ', redshifts[n])
		for a in range(0, Naps):	#this goes through every aperture
			ns=str(a)
			#print('aperture0',ns)
			absg, absr, absi, absz, absy= abs_mag(newdata[n], 'mag_aperture0', kcorrect, DM[n], bands, ns, n) 
			Lumg, Lumr, Lumi, Lumz, Lumy=abs2lum(absg, absr, absi, absz, absy)			
			Lg, Lr, Li, Lz, Ly=lumdensity(Lumg, Lumr, Lumi, Lumz, Lumy, radkpc[a])
			if my_band=='i':
				if scale== 'log':
					#print('getting logs')
					logLum=math.log10(Lumi)
					logL=math.log10(Li)
					L.append(logLum)
					lumd.append(logL)
				else:
					L.append(Lumi)
					lumd.append(Li)
			elif my_band=='g':
				if scale== 'log':
					#print('getting logs')
					logLum=math.log10(Lumg)
					logL=math.log10(Lg)
					L.append(logLum)
					lumd.append(logL)
				else:
					L.append(Lumg)
					lumd.append(Lg)
		#print('LI for ',n,' galaxy is ', LI)
		bigL.append(L)
		bigden.append(lumd)
		if scale== 'log':
			lograd=[math.log10(radkpc[n]) for n in range(len(radkpc))]
			bigrad.append(lograd)
		else:
			bigrad.append(radkpc)
	bigLs=np.array(bigL)
	bigrads=np.array(bigrad)
	lumdens=np.array(bigden)
	return bigLs, bigrads, lumdens
	
def upper_rad_cut(lum, rad, den, m, proof=False):
	from def_halflight_math import get_halfrad #This is done linearly, so we dont have to do much
	print('this is from halflight3')
	nlum=[]
	nrad=[]
	nden=[]
	mult=m
	#print(len(lum), len(rad)) #make sure they are same length
	N=len(rad)
	for n in range(0,N):
		L=lum[n]
		R=rad[n]
		dens=den[n]
		r12=get_halfrad(R, L) #lower limit
		r412=mult*r12 #upper limit
		if proof == True:
			print('r1/2= ', r12)
			print('4r1/2= ', r412)
			print('The radii are ', R)
		if np.max(R) >= r412:
			rcut=R[(R>=r12)&(R<=r412)]
			if proof == True:
				print('Cut Radius range= ', rcut)
			if len(rcut)>=4:
				nlum.append(L)
				nrad.append(R)
				nden.append(dens)
				print('good')
			else:
				print('not enough data points')
		else:
			print('Upper limit out of range')
		#break
	nlum=np.array(nlum)
	nrad=np.array(nrad)
	nden=np.array(nden)
	return nlum, nrad, nden

def get_avg_lums(Ls, rads, LDs, gr=[], type=''):
	print('get_avg_lums is in halflight3')
	
	Naps=0.0
	
	if type=='mean':
		avglum, radavg, bb=meanlum2(Ls, rads, Naps,grange=gr)
		avgdens, radavg, bb=meanlum2(LDs, rads, Naps,grange=gr)
		
		err='bootstrap_stdv'
		logdenerr=get_errors(LDs, rads, bb, avgdens, error=err)
	
		print('Mean Luminosity= ', avglum)
		print('Mean LumDensity=', avgdens)
		print('Binned Radii= ', radavg)
		print('Standard Deviation (in log)= ', logdenerr)

		
	if type=='med':
		print('getting median')
		avglum, radavg, bb=medlum2(Ls, rads, Naps,grange=gr)
		avgdens, radavg, bb=medlum2(LDs, rads, Naps,grange=gr)
	
		err='bootstrap_stdv'
		logdenerr=get_errors(LDs, rads, bb, avgdens, error=err)
	
		print('Median Luminosity= ', avglum)
		print('Median LumDensity=', avgdens)
		print('Binned Radii= ', radavg)
		print('Standard Deviation (in log)= ', logdenerr)
		#hi=hi
		
	errtest='' #this testplots error bars amongst individual galaxies
	if errtest:
		print(len(LDs))
		print('log of ',type ,' lumdens= ', np.log10(avgdens))
		print('log of radavg= ', np.log10(radavg))
		for n in range(len(LDs)):
			plt.plot(np.log10(rads[n]), np.log10(LDs[n]), color='lightgrey', marker='.', zorder=1)
		plt.scatter(np.log10(radavg), np.log10(avgdens), zorder=2, color='red', marker='o')
		plt.errorbar(np.log10(radavg), np.log10(avgdens), yerr=logdenerr, zorder=4, color='red')
		plt.xlabel('Log Radii')
		plt.ylabel('Log Luminosity Density')
		plt.title(type)
		plt.show()
	return avglum, avgdens, radavg, logdenerr #outputs logmeans and log mean_er
		
def get_halflight2(Ls, rads, mult=4, inn=1, maxout=None, maxin=None):
	import math
	import numpy as np
	print('from halflight3')
	N=np.ndim(rads)
	if N==2:
		r12s=[]
		r412s=[]
		for n in range(0, len(rads)):
			L=Ls[n]
			R=rads[n]
			ar12=get_halfrad(R, L)
			if maxout==None:
				ar412=ar12*mult
				print('This is 4r1/2= ', ar412)
			else:
				ar412=maxout
				print('maximum outer range is ', ar412)
			if maxin==None:
				print('This is r1/2= ', ar12)
				ar12s=inn*ar12
			else:
				ar12s=maxin
			print('This is the inner limit= ', ar12s)
			r12s.append(ar12s)
			r412s.append(ar412)
		r12=np.array(r12s)
		r412=np.array(r412s)
	else:
		r12s=get_halfrad(rads, Ls)
		print('Stacked halflight radius : ', r12s)
		if maxout==None:
			r412=r12s*mult
		else:
			r412=maxout
		if maxin==None:
			r12=inn*r12s
		else:
			r12=maxin
		print('Inner stacked radius: ', r12)
		print('Outer stacked Radius: ', r412)
	return r12, r412

def get_slopes1(r12s, r412s, rad, ld, error=None, scale='', smax=False):
	import scipy.stats as stats
	from def_halflight_math import my_linregress3
	print('slopes from halflight3')
	mult=4
	Ndim=np.ndim(rad)
	P=len(rad)
	print('Shape of radius array= ', np.shape(rad))
	if error is None:
		print('No error was given')
		error=np.ones((P, len(rad[0])))
	rcut=[]
	ldcut=[]
	errcut=[]
	if Ndim==2:
		for i in range(P):
			rrow=rad[i]
			ldrow=ld[i]
			errow=error[i]
			r12=r12s[i]
			r412=r412s[i]
			if smax== True:
				if np.max(rrow) >= r412:
					mr=rrow[(rrow>=r12)&(rrow<=r412)]
					mld=ldrow[(rrow>=r12)&(rrow<=r412)]
					merr=errow[(rrow>=r12)&(rrow<=r412)]
					if len(mr) >=3:
						print('good!')
						rcut.append(mr)
						ldcut.append(mld)
						errcut.append(merr)
					else:
						print('not enough data points')
				else:
					print('Upper Cut is Out of the Radius Range')
			else:
				merr=errow[rrow>=r12]
				mr=rrow[rrow>=r12]
				mld=ldrow[rrow>=r12]
				if len(mr) >=3:
					print('good!')
					rcut.append(mr)
					ldcut.append(mld)
					errcut.append(merr)
		rcut=np.array(rcut)
		ldcut=np.array(ldcut)
		errcut=np.array(errcut)
		slopes=[]
		intercepts=[]
		errs=[]
		if scale=='linear':
			for n in range(len(rcut)):
				slope, int, std_err=my_linregress3(rcut[n], ldcut[n], errcut[n])
				slopes.append(slope)
				intercepts.append(int)
				errs.append(std_err)
		else:
			print('inputting log data for fitting') #we want to fit in log-log
			#error is already in log
			for n in range(len(rcut)):
				logrcut=np.log10(rcut[n])
				logldcut=np.log10(ldcut[n])
				slope, int, std_err=my_linregress3(logrcut, logldcut, errcut[n])
				slopes.append(slope)
				intercepts.append(int)
				errs.append(std_err)
		slopes=np.array(slopes)
		intercepts=np.array(intercepts)
		errs=np.array(errs)
		return slopes, intercepts, errs #these are all in
	else: #for arrays of 1D *aka* the stacked profile
		print('1D array!')
		rad=np.array(rad)
		r12=r12s
		r412=r412s
		if error is None:
			error=np.ones(P)
		if smax== True:
			print('upper limit is ', r412)
			if np.max(rad) < r412:
				print('Upper cut is out of the Radius range')
			else:
				rcut=rad[(rad>=r12)&(rad<=r412)]
				ldcut=ld[(rad>=r12)&(rad<=r412)]
				errcut=error[(rad>=r12)&(rad<=r412)]
		else:
			rcut=rad[rad>=r12]
			ldcut=ld[rad>=r12]
			errcut=error[rad>=r12]
		print('Stacked Radii cut are= ', rcut)
		print('Stacked R1/2 is= ', r12)
		print('What is the error cut? ', errcut)
		try:
			print('Stacked 4R1/2 is= ', r412)
		except:
			pass
		if scale=='linear':
			sl3, C3, std_err3=my_linregress3(rcut, ldcut, errcut)
		else:
			print('inputting log data for fitting')
			rcut=np.log10(rcut)
			ldcut=np.log10(ldcut)
			sl3, C3, std_err3=my_linregress3(rcut, ldcut, errcut)
		return sl3, C3, rcut, ldcut, std_err3, errcut