print('Will get absolute magnitudes from this')

import astropy.table as table 
import numpy as np
from defcuts import *
from def_get_mags import *
from my_def_plots import *
from defflags import many_flags
from defclump import meanlum
import matplotlib.pyplot as plt

indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'

outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'

datatab = table.Table.read(indir+ 'LOWZ_HSCGAMA15_apmgs.fits')

bands=['g', 'r', 'i','z', 'y']

parm=['flags_pixel_saturated_center','flags_pixel_edge','flags_pixel_interpolated_center','flags_pixel_cr_center','flags_pixel_suspect_center', 'flags_pixel_clipped_any','flags_pixel_bad']

daperture=[1.01,1.51,2.02,3.02,4.03,5.71,8.40,11.8,16.8,23.5]

aperture=[x*0.5 for x in daperture]

#get rid of cuts
#mag_cmodel
ne=[99.99, 199.99, 0.0]

mincut=0.1
maxcut=''
cutdata=not_cut(datatab, bands, 'mag_aperture00', ne)



#get rid of flagged galaxies
for b in range(0, len(bands)-1):
	newdata=many_flags(cutdata, parm, bands[b])	#I think flags are only in band i
	cutdata=newdata

newdata=cutdata

Naps=len(aperture)

objID=newdata['object_id']
redshift=newdata['Z']
Ndat=len(redshift)
DM= get_zdistmod(newdata, 'Z')

kcorrect=get_kcorrect2(newdata,'mag_aperture0', '_err', bands, '0', 'hsc_filters.dat', redshift)

#newdata=nan_cut(newdata, kcorrect)
#print(kcorrect[1])		
#print(kcorrect[1][0])	first number is row, second is element in row
#print(kcorrect[2])
#print(kcorrect[2][2])
#for n in range(0, Ndat,10):

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
	radkpc=aper_and_comov(aperture, redshift[n])
	print('Redshift is ', redshift[n])
	for a in range(0, Naps):
		#this goes through every aperture
		ns=str(a)
		print('aperture0',ns)
	
		#get magnitude
		absg, absr, absi, absz, absy= abs_mag(newdata[n], 'mag_aperture0', kcorrect, DM[n], bands, ns, n) 
	
		Lumg, Lumr, Lumi, Lumz, Lumy=abs2lum(absg, absr, absi, absz, absy)
		
		Lg, Lr, Li, Lz, Ly=lumdensity(Lumg, Lumr, Lumi, Lumz, Lumy, radkpc[a])
		
		LG.append(Lg)
		LR.append(Lr)
		LI.append(Li)
		LZ.append(Lz)
		LY.append(Ly)
	print('LI for ',n,' galaxy is ', LI)
	bigLI.append(LI)
	bigrad.append(radkpc)
	plt.plot(radkpc, LI, marker='.', color=grays, alpha=0.9)
#print('Mins and Maxes: ', min(bigrad), max(bigrad))
mean, error, radavg=meanlum(bigLI, bigrad, Naps, scale='log')
plt.plot(radavg, mean, marker='.', color='r')
plt.errorbar(radavg, mean, yerr=error, color='m', fmt='.')
plt.xlabel('Comoving Distance (kpc)', fontsize=10)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Luminosity Density (Lsolar/kpc^2)', fontsize=10)
plt.suptitle('Luminosity Density vs. Comoving Distance in band I', fontsize=15)
plt.plot(0,0,label='Number of galaxies='+str(Ndat), c='k', marker='')
plt.plot(0,0,label='Standard Deviations: ', c='m')
for b in range(0, Naps):
	bs=str(b)
	errors=round(error[b],4)
	errstr=str(errors)
	plt.plot(0,0,label='Aperture'+bs+'= '+errstr, c='m', marker='')

plt.legend(loc=1,prop={'size':5.5})
plt.savefig(outdir+'istack_lumdens_prof.pdf')


#will eventually need comoving