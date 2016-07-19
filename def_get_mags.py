import numpy as np
import math
from astropy.cosmology import FlatLambdaCDM
def get_zdistmod(data, colname):
	print('We are going to get the distance modulus given a z shift')
	cosmo = FlatLambdaCDM(H0=100, Om0=0.3)
	Z=data[colname]
	DMs= cosmo.distmod(Z)
	DM=DMs.value
	#print('Redshift is ', Z)
	#print(len(DM))
	#print('Quantity of DM is ', DMs)
	return DM

def get_kcorrect(data,mag, magerr, bands, filter, redshift):
	import kcorrect
	
	#magnitudes in different bands
	i=data[bands[2]+mag+aps]
	ierr=data[bands[2]+mag+aps+magerr]
	r=data[bands[1]+mag+aps]
	rerr=data[bands[1]+mag+aps+magerr]
	g=data[bands[0]+mag+aps]
	gerr=data[bands[0]+mag+aps+magerr]
	z=data[bands[3]+mag+aps]
	zerr=data[bands[3]+mag+aps+magerr]
	y=data[bands[4]+mag+aps]
	yerr=data[bands[4]+mag+aps+magerr]
	
	#print(i, r, g, z)
	
	maggiesi=10**(-0.4*i)
	maggiesr=10**(-0.4*r)
	maggiesg=10**(-0.4*g)
	maggiesz=10**(-0.4*z)
	maggiesy=10**(-0.4*y)
	
	magiivar=1./((10**(-0.4*(ierr)))**2)/(0.4*math.log(10.))**2/maggiesi**2
	maggivar=1./((10**(-0.4*(gerr)))**2)/(0.4*math.log(10.))**2/maggiesg**2
	magrivar=1./((10**(-0.4*(rerr)))**2)/(0.4*math.log(10.))**2/maggiesr**2
	magzivar=1./((10**(-0.4*(zerr)))**2)/(0.4*math.log(10.))**2/maggiesz**2
	magyivar=1./((10**(-0.4*(yerr)))**2)/(0.4*math.log(10.))**2/maggiesy**2

	print('is loading the templates taking a while')
	kcorrect.load_templates()
	print('or is loading the filters?')	#<---- it's this
	kcorrect.load_filters(f=filter)
	
	maggies = np.array([maggiesg, maggiesr, maggiesi, maggiesz, maggiesy], dtype='float32')
	maggies_ivar = np.array([maggivar, magrivar, magiivar, magzivar, magyivar], dtype='float32')
	coeffs = kcorrect.fit_nonneg(redshift, maggies, maggies_ivar)
	
	rm = kcorrect.reconstruct_maggies(coeffs)	#reconstructed maggies
	rm0 = kcorrect.reconstruct_maggies(coeffs, redshift=0.)
	kc = -2.5*np.log10(rm[1:]/rm0[1:])
	
	print('kcorrect in g, r, i, z, y= ',kc)
	#in order g, r, i, z, y
	return kc
	
def abs_mag(data, mag, kcorrect, DM, bands, aps, n):
#	print('Gets absolute magnitude of a single galaxy')
	#bands=['kg', 'kr', 'ki','kz', 'ky']
	if aps:
		#get apparent magnitude
		im=data[bands[2]+mag+aps]
		rm=data[bands[1]+mag+aps]
		gm=data[bands[0]+mag+aps]
		zm=data[bands[3]+mag+aps]
		ym=data[bands[4]+mag+aps]
	else:
		im=data[bands[2]+mag]
		rm=data[bands[1]+mag]
		gm=data[bands[0]+mag]
		zm=data[bands[3]+mag]
		ym=data[bands[4]+mag]
	#print('apparent mags in g, r, i, z, y= ', gm, rm, im, zm, ym)
	#M=m-DM-K
	#print('Using K corrects: ', kcorrect[n])
	absi=im-DM-kcorrect[n][2]
	absr=rm-DM-kcorrect[n][1]
	absg=gm-DM-kcorrect[n][0]
	absz=zm-DM-kcorrect[n][3]
	absy=ym-DM-kcorrect[n][4]
	
	#print('abs mags in g, r, i, z, y= ', absg, absr, absi, absz, absy)
	#print('Absolute mag i:', absi,'= ', i,'- ',DM, '- ', kcorrect[2])
	return absg, absr, absi, absz, absy

def abs2lum(absg, absr, absi, absz, absy):
	print('Getting Luminosities in Solar Luminosities')
	sung = 5.07487
	sunr = 4.64305
	suni = 4.52725
	sunz = 4.50954
	suny = 4.50291
	
	lumg=10**(-0.4*(absg-sung))
	lumr=10**(-0.4*(absr-sunr))
	lumi=10**(-0.4*(absi-suni))
	lumz=10**(-0.4*(absz-sunz))
	lumy=10**(-0.4*(absy-suny))
	
	#print('Luminosity in i: ', lumi)
		
	return lumg, lumr, lumi, lumz, lumy

def aper_and_comov(aper, redshift):
	pi=math.pi
	cosmo = FlatLambdaCDM(H0=100, Om0=0.3)
	
	Dcs=cosmo.comoving_distance(redshift)	#get distance in Mpc
	
	Dc=Dcs.value
	Dckpc=Dc*1000
	
	aperkpc=[]
	
	N=len(aper)
	
	for n in range(0, N):
		new=aper[n]/3600.0*pi/180.0*Dckpc
		aperkpc.append(new)
	#	print('New aperture length is ', new, 'kpc')
	
	return aperkpc

def lumdensity(lg, lr, li, lz, ly, radap):	
	pi=math.pi
	
	ldg=lg/(4*pi*radap**2)
	ldr=lr/(4*pi*radap**2)
	ldi=li/(4*pi*radap**2)
	ldz=lz/(4*pi*radap**2)
	ldy=ly/(4*pi*radap**2)
	
	return ldg, ldr, ldi, ldz, ldy
	
def get_kcorrect2(data,mag, magerr, bands, aps, filter, redshift):
	import kcorrect
	
	N=len(redshift)
	
	zshift=np.array(redshift,dtype='float32')
	
	kcorrect.load_templates()
	print('It takes a while to load the filters...')	#<---- it's this
	kcorrect.load_filters(f=filter)
	
	#magnitudes in different bands
	kc=[]
	for p in range(0,N):
	
		i=data[bands[2]+mag+aps][p]
		ierr=data[bands[2]+mag+aps+magerr][p]
		r=data[bands[1]+mag+aps][p]
		rerr=data[bands[1]+mag+aps+magerr][p]
		g=data[bands[0]+mag+aps][p]
		gerr=data[bands[0]+mag+aps+magerr][p]
		z=data[bands[3]+mag+aps][p]
		zerr=data[bands[3]+mag+aps+magerr][p]
		y=data[bands[4]+mag+aps][p]
		yerr=data[bands[4]+mag+aps+magerr][p]
	
		#print(i, r, g, z)
	
		maggiesi=10**(-0.4*i)
		maggiesr=10**(-0.4*r)
		maggiesg=10**(-0.4*g)
		maggiesz=10**(-0.4*z)
		maggiesy=10**(-0.4*y)
	
		magiivar=1./((10**(-0.4*(ierr)))**2)/(0.4*math.log(10.))**2/maggiesi**2
		maggivar=1./((10**(-0.4*(gerr)))**2)/(0.4*math.log(10.))**2/maggiesg**2
		magrivar=1./((10**(-0.4*(rerr)))**2)/(0.4*math.log(10.))**2/maggiesr**2
		magzivar=1./((10**(-0.4*(zerr)))**2)/(0.4*math.log(10.))**2/maggiesz**2
		magyivar=1./((10**(-0.4*(yerr)))**2)/(0.4*math.log(10.))**2/maggiesy**2
	
		maggies = np.array([maggiesg, maggiesr, maggiesi, maggiesz, maggiesy], dtype='float32')
		maggies_ivar = np.array([maggivar, magrivar, magiivar, magzivar, magyivar], dtype='float32')
		coeffs = kcorrect.fit_nonneg(redshift[p], maggies, maggies_ivar)
	
		rm = kcorrect.reconstruct_maggies(coeffs)	#reconstructed maggies
		rm0 = kcorrect.reconstruct_maggies(coeffs, redshift=0.)
		kcs = -2.5*np.log10(rm[1:]/rm0[1:])
	
		#print('kcorrect in g, r, i, z, y= ',kcs)
		kc.append(kcs)
	#print(kc)
	return kc