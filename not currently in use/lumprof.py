print('This program should make plots of aperture size v aperture maginosity')

indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.text as txt
from astropy.io import fits
import astropy.table as table 
import matplotlib.patches as mpatches
from sympy import *
import os, sys

datatab = table.Table.read(indir+ 'LOWZ_HSCGAMA15_apmgs.fits')

#datatab=fits.open(indir+ 'LOWZ_HSCGAMA15_apmgs.fits') <-- probably better to use table
#hdulist=fits.open(indir+ 'LOWZ_HSCGAMA15_apmgs.fits')
#hdulist.info()
#print(datatab['Z'])
#datatab['imag_aperture00'] gives first aperture maginosity in band i
#
aperture=[3.0*0.168,4.5*0.168,6.0*0.168,9.0*0.168,12.0*0.168,17.0*0.168,25.0*0.168,35.0*0.168,50.0*0.168,70.0*0.168]
pi=math.pi
Naps=len(aperture)
N=len(datatab)

band= ['r','i','g', 'z', 'y']
band=band[0]

no_flag=0
no_sat=0
no_edge=1

if no_flag:
	flag=('No Flag', 'no_flag')
	num=len(datatab)
	Nnew=num
	num=str(num)
if no_sat:
	datatab_flag = datatab[(datatab[band+'flags_pixel_saturated_center']==False) & (datatab[band+'mag_aperture00']<50)& (datatab[band+'mag_aperture01']<50)& (datatab[band+'mag_aperture02']<50) & (datatab[band+'mag_aperture03']<50) & (datatab[band+'mag_aperture04']<50) & (datatab[band+'mag_aperture05']<50) & (datatab[band+'mag_aperture06']<50) & (datatab[band+'mag_aperture07']<50) & (datatab[band+'mag_aperture08']<50) & (datatab[band+'mag_aperture09']<50)]
	datatab=datatab_flag
	num=len(datatab)
	Nnew=num
	num=str(num)
	flag=('No Saturated Centers', 'no_satcen')
if no_edge:
	datatab_flag = datatab[(datatab[band+'flags_pixel_edge']==False) & (datatab[band+'mag_aperture00']<50)& (datatab[band+'mag_aperture01']<50)& (datatab[band+'mag_aperture02']<50) & (datatab[band+'mag_aperture03']<50) & (datatab[band+'mag_aperture04']<50) & (datatab[band+'mag_aperture05']<50) & (datatab[band+'mag_aperture06']<50) & (datatab[band+'mag_aperture07']<50) & (datatab[band+'mag_aperture08']<50) & (datatab[band+'mag_aperture09']<50)]
		#this should hopefully get rid of random bad pixels
	datatab=datatab_flag
	num=len(datatab)
	Nnew=num
	num=str(num)
	flag=('No Edge Galaxies', 'no_edge')



outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/single_plot/'
if not os.path.exists(outdir):
	os.mkdir(outdir)

	
for i in range(0,Nnew):	#this goes through every galaxy
	objID=datatab['object_id']
	name=objID[i]
	name=str(name)
	print(Nnew, N)
	print(name)
	magap=[]
	SB=[]
	magerr=[]
	SBerr=[]
	for j in range(0,Naps):	#this goes through every aperture for maginosity
		j=str(j)
		mag=datatab[band+'mag_aperture0'+j][i]
		merr=datatab[band+'mag_aperture0'+j+'_err'][i]
		magap.append(mag)
		magerr.append(merr)
	for j in range(0,Naps):	#this goes through every aperture for Sb
		js=str(j)
		mag=datatab[band+'mag_aperture0'+js][i]
		merr=datatab[band+'mag_aperture0'+js+'_err'][i]
		sb=mag+2.5*math.log10(4*pi*aperture[j]**2)
		sberr=diff(sb)
		#F=math.pow(10, -mag/2.5)
		#sb=-2.5*math.log10(F/(4*pi*aperture[j]**2))
		#Ferr=math.pow(10, -merr/2.5)
		#sberr=-2.5*math.log10(Ferr/(4*pi*aperture[j]**2))
		SB.append(sb)
		SBerr.append(sberr)
	fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)
	ax0.scatter(aperture, magap, c='k', marker='^')
	ax0.errorbar(aperture, magap, yerr=magerr, marker='', mfc='red', mec='green', ms=5, mew=1)
	ax0.invert_yaxis()
	ax0.set_xlabel('Aperture Radius (arcseconds)', fontsize=9)
	ax0.set_ylabel('Aperture Magnitude', fontsize=9)
	ax0.set_title("Luminosity Profiles vs. Aperture Radius in "+band+" ("+name+')', fontsize=11)
	ax0.set_xlim(xmin=0, xmax=max(aperture))
	plt.plot(0,0,label=flag[0]+'('+num+')', marker='', c='k')
	plt.legend(loc=4,prop={'size':4})
	ax1.scatter(aperture, SB, c='k', marker='^')
	ax1.errorbar(aperture, SB, yerr=SBerr, marker='', mfc='red', mec='green', ms=5, mew=1)
	ax1.invert_yaxis()
	ax1.set_xlabel('Aperture Radius (arcseconds)', fontsize=9)
	ax1.set_ylabel('Surface Brightness(mag/arcsec^2)', fontsize=9)
	ax1.set_title("Surface Brightness vs. Aperture Radius in "+band+" ("+name+')', fontsize=11)

	#plt.show()
	fig.savefig(outdir+flag[1]+'_'+band+name+'_lumprof.pdf')
	#break