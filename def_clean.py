import numpy as np
from def_get_mags import get_zdistmod, get_kcorrect2, aper_and_comov, abs2lum, lumdensity, abs_mag
import math
from defclump import meanlum2, medlum2, get_error
import matplotlib.pyplot as plt
from def_mymath import halflight

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
	return bigLIs, bigrads, lumdensi
	
def get_avg_lums(bigLIs, bigrads, lumdensi, type=''):
	bigLIs.flatten()	#luminosity
	bigrads.flatten()	#radii
	lumdensi.flatten()	#luminosity density
	Naps=0.0
	
	if type=='mean':
		meanlum, radavg, bb=meanlum2(bigLIs, bigrads,Naps,scale='log')
		meandens, radavg, bb=meanlum2(lumdensi, bigrads,Naps,scale='log')
	
		err='bootstrap_stdv'
		lumdenerr=get_error(lumdensi, bigrads, bb, error=err)
	
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

def get_halflight(bigLIs, bigrads):
	print('Getting Halfs')
	N=np.ndim(bigrads)
	if N == 2:
		print('Array is 2D')
		halfrad=[]
		for x in range(len(bigrads)):
			#print('Galaxy #', str(x))
			#print('Luminosity= ', bigLIs[x])
			#print('Radii= ', bigrads[x])
	
			half=math.log10(10**np.max(bigLIs[x])/2.0)
			fhalf=halflight(bigrads[x],bigLIs[x])
			halfrad.append(fhalf)
		halfrads=np.array(halfrad)
	else:
		halfrads=halflight(bigrads,bigLIs)
	return halfrads

def get_slopes(lum, hx, x, y, error=None, names=None, smax=False):
	import scipy.stats as stats
	from def_mymath import my_linregress3
	from my_def_plots import scatter_fit, simple_hist
	
	mult=4
	Ndim=np.ndim(x)
	N=len(x)
	
	if error is None:
		print('No error was given')
		error=np.ones((N, len(x[0])))
	
	if names is None:
		print('No names given')
	if Ndim==2:
		xcut=[]
		ycut=[]
		errcut=[]
		for i in range(len(x)):
			xrow=x[i]
			yrow=y[i]
			errow=error[i]
			hhx=hx[i]
			
			#merr=errow[xrow>=hhx]
			#mx=xrow[xrow>=hhx]
			#my=yrow[xrow>=hhx]
			
			if smax is True:
				hhx10=10**hhx
				hhx2s=mult*hhx10
				hhx2=math.log10(hhx2s)
				
				bad=0
				print(hhx2, np.max(xrow))
				if np.max(xrow) >= hhx2:
					mx=xrow[(xrow>=hhx)&(xrow<=hhx2)]
					my=yrow[(xrow>=hhx)&(xrow<=hhx2)]
					merr=errow[(xrow>=hhx)&(xrow<=hhx2)]
					if len(mx) >=4:
						xcut.append(mx)
						ycut.append(my)
						errcut.append(merr)
				else:
					print('Upper Cut is Out of the Radius Range')
					bad += 1
			elif smax is False:
				merr=errow[xrow>=hhx]
				mx=xrow[xrow>=hhx]
				my=yrow[xrow>=hhx]
				if len(mx) >=4:
					xcut.append(mx)
					ycut.append(my)
					errcut.append(merr)
		slopes=[]
		intercepts=[]
		errs=[]
		for n in range(len(xcut)):
			#print(len(x[n]), len(xcut[n]))
			#slope, int, r_value, p_value, std_err = stats.linregress(xcut[n],ycut[n])
			slope, int, std_err=my_linregress3(xcut[n], ycut[n], errcut[n])
			slopes.append(slope)
			intercepts.append(int)
			errs.append(std_err)
		return slopes, intercepts, errs
	else:
		x=np.array(x)
		if error is None:
			error=np.ones(N)
		if smax is True:
			hx10=10**hx
			hx2s=mult*hx10
			hx2=math.log10(hx2s)
			
			if np.max(x) <= hx2:
				print('Upper cut is out of the Radius range')
			else:
				xcut=x[(x>=hx)&(x<=hx2)]
				ycut=y[(x>=hx)&(x<=hx2)]
				errcut=error[(x>=hx)&(x<=hx2)]
		else:
			xcut=x[x>=hx]
			ycut=y[x>=hx]
			errcut=error[x>=hx]
		print('Radii are= ', x)
		print('R1/2 is= ', hx)
		
		sl3, C3, std_err3=my_linregress3(xcut, ycut, errcut)
		#y4 = sl3 * xcut + C3
		return sl3, C3, xcut, ycut, std_err3, errcut
		
def get_slopes1(data, lum, hx, x, y, error=None, names=None, smax=False):
	import scipy.stats as stats
	from def_mymath import my_linregress3
	from my_def_plots import scatter_fit, simple_hist
	get_vir_r= lambda M, ro, delta_c: (M/(4./3.*math.pi*rho*delta_c))**(1./3.)
	mult=5
	Ndim=np.ndim(x)
	N=len(x)
	if error is None:
		print('Is this going through?')
		error=np.ones((N, len(x[0])))
	if names is None:
		print('No names given')
	if Ndim==2:
		xcut=[]
		ycut=[]
		errcut=[]
		bad=0
		good=0
		for i in range(len(x)):
			xrow=x[i]
			yrow=y[i]
			Lum=lum[i]
			errow=error[i]
			hhx=hx[i]
			
			if smax is True:
				from astropy.cosmology import FlatLambdaCDM
				from astropy import units as u
				cosmo = FlatLambdaCDM(H0=100, Om0=0.3)
				redshift=data['Z']
				rhos=cosmo.critical_density(redshift[i]) #in g/cm^3
				solMkpc3= u.solMass / u.kpc**3
				rho=rhos.to(solMkpc3)
				print('Critical Density= ', rho)
				rho=rho.value
				M2L=10
				
				Lmax=np.max(Lum)
				
				Mvir=M2L*10**Lmax
				#Mvir=math.log10(Mvir)
				#Mvir=M2L*Lmax
				print('Virial Mass= ', Mvir)
				
				delta_c=500
				
				r_vir=get_vir_r(Mvir, rho, delta_c)
				r_vir=math.log10(r_vir)
				print('Radii= ', str(xrow))
				print('Virial Radius is ', str(r_vir))
				if np.max(xrow) >= r_vir:
					mx=xrow[(xrow>=hhx)&(xrow<=r_vir)]
					my=yrow[(xrow>=hhx)&(xrow<=r_vir)]
					merr=errow[(xrow>=hhx)&(xrow<=r_vir)]
					if len(mx) >=4:
						xcut.append(mx)
						ycut.append(my)
						errcut.append(merr)
						good +=1
					else:
						print('Not enough data')	
				else:
					print('Upper Cut is Out of the Virial Radius Range')
					#bad=sum(bad+1)
					bad +=1
			elif smax is False:
				merr=errow[xrow>=hhx]
				mx=xrow[xrow>=hhx]
				my=yrow[xrow>=hhx]
				if len(mx) >=4:
					xcut.append(mx)
					ycut.append(my)
					errcut.append(merr)
				else:
					print('Not enough data')
		print('Number of good profiles: ', good)
		print('Number of bad profiles: ', bad)
		print(len(x), len(xcut))
		hi=hi
		slopes=[]
		intercepts=[]
		errs=[]
		for n in range(len(xcut)):
			slope, int, std_err=my_linregress3(xcut[n], ycut[n], errcut[n])
			slopes.append(slope)
			intercepts.append(int)
			errs.append(std_err)
		return slopes, intercepts, errs
	else:
		x=np.array(x)
		if error is None:
			error=np.ones(N)
		if smax is True:
			hx10=10**hx
			hx2s=mult*hx10
			hx2=math.log10(hx2s)
			
			if np.max(x) <= hx2:
				print('Upper cut is out of the Radius range')
			else:
				xcut=x[(x>=hx)&(x<=hx2)]
				ycut=y[(x>=hx)&(x<=hx2)]
				errcut=error[(x>=hx)&(x<=hx2)]
		else:
			xcut=x[x>=hx]
			ycut=y[x>=hx]
			errcut=error[x>=hx]
		
		sl3, C3, std_err3=my_linregress3(xcut, ycut, errcut)
		return sl3, C3, xcut, ycut, std_err3, errcut
	