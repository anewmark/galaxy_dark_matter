print('Will plot single galaxy luminosity density profiles')

import astropy.table as table 
from defcuts import *
from def_get_mags import *
from my_def_plots import *
from defflags import many_flags
import matplotlib.pyplot as plt


indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'

outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/single_plot/'

datatab = table.Table.read(indir+ 'LOWZ_HSCGAMA15_apmgs.fits')

bands=['g', 'r', 'i','z', 'y']

parm=['flags_pixel_saturated_center','flags_pixel_edge','flags_pixel_interpolated_center','flags_pixel_cr_center','flags_pixel_suspect_center', 'flags_pixel_clipped_any','flags_pixel_bad']

daperture=[1.01,1.51,2.02,3.02,4.03,5.71,8.40,11.8,16.8,23.5]

aperture=[x*0.5 for x in daperture]

#get rid of cuts
mincut=0.1
maxcut=''
cutdatag, crange=out_cut(datatab, bands[0], 'mag_aperture00',mincut, maxcut)
cutdatai, crange=out_cut(cutdatag, bands[2], 'mag_aperture00',mincut, maxcut)
cutdatar, crange=out_cut(cutdatai, bands[1], 'mag_aperture00',mincut, maxcut)
cutdatay, crange=out_cut(cutdatar, bands[4], 'mag_aperture00', mincut, maxcut)
cutdataz, crange=out_cut(cutdatay, bands[3], 'mag_aperture00',mincut, maxcut)

ne=[199.99, 99.99]
cutdata1=not_cut(cutdataz, bands[1], 'mag_aperture00', ne)
cutdata=not_cut(cutdata1, bands[0], 'mag_aperture00', ne)

#get rid of flagged galaxies
newdata=many_flags(cutdata, parm, 'i')	#I think flags are only in band i

Naps=len(aperture)

objID=newdata['object_id']
redshift=newdata['Z']
Ndat=len(redshift)
DM= get_zdistmod(newdata, 'Z')

kcorrect=get_kcorrect2(newdata,'mag_aperture0', '_err', bands, '0', 'hsc_filters.dat', redshift)

#for n in range(0, Ndat,10):

for n in range(0,Ndat):
	#this goes through every galaxy
	name=objID[n]
	name=str(name)
	LG=[]
	LR=[]
	LI=[]
	LZ=[]
	LY=[]
	
	radkpc=aper_and_comov(aperture, redshift[n])
	for a in range(0, Naps):
		#this goes through every aperture
		ns=str(a)
		print(ns)
	
		#get magnitude
		absg, absr, absi, absz, absy= abs_mag(newdata[n], 'mag_aperture0', kcorrect, DM[n], bands, ns, n) 
	
		Lumg, Lumr, Lumi, Lumz, Lumy=abs2lum(absg, absr, absi, absz, absy)
		
		Lg, Lr, Li, Lz, Ly=lumdensity(Lumg, Lumr, Lumi, Lumz, Lumy, radkpc[a])
		
		LG.append(Lg)
		LR.append(Lr)
		LI.append(Li)
		LZ.append(Lz)
		LY.append(Ly)
		
	#	break
	#creating luminosity densities for the apertures at each band
	print('Galaxy # ', n)
	lum_comov_plot(LG, LR, LI, LZ, LY, radkpc, name, outdir)

#will eventually need comoving