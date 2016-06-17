print('This program should plot average flux vs. aperture size for different bands')

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

aperture=[3.0*0.168,4.5*0.168,6.0*0.168,9.0*0.168,12.0*0.168,17.0*0.168,25.0*0.168,35.0*0.168,50.0*0.168,70.0*0.168]
pi=math.pi
objID=datatab['object_id']
Naps=len(aperture)
N=len(datatab)


bands= ['r','g', 'i', 'y', 'z']

for band in bands:
	#getting the no flagged plots in z range[0.15,0.25]
	Z=datatab['Z']
	zmin=0.15
	zmax=0.25

	dataclump=datatab[(Z>=zmin) & (Z<zmean)]
	range=[zmin, zmean]
	range=str(range)

	nclump=len(dataclump1)
	nclump2=len(dataclump2)
	nclump3=len(dataclump3)
	nclump4=len(dataclump4)
	SBclump1=[]
	SBclump2=[]
	SBclump3=[]
	SBclump4=[]
	for j in range(0,Naps):
		js=str(j)
		magclump1=dataclump1[band+'mag_aperture0'+js]
		sbclump1=magclump1+2.5*math.log10(4*pi*aperture[j]**2)
		sbsum1=np.sum(sbclump1)
		sbavg1=sbsum1/nclump1
		SBclump1.append(sbavg1)
		magclump2=dataclump2[band+'mag_aperture0'+js]
		sbclump2=magclump2+2.5*math.log10(4*pi*aperture[j]**2)
		sbsum2=np.sum(sbclump2)
		sbavg2=sbsum2/nclump2
		SBclump2.append(sbavg2)
		magclump3=dataclump3[band+'mag_aperture0'+js]
		sbclump3=magclump3+2.5*math.log10(4*pi*aperture[j]**2)
		sbsum3=np.sum(sbclump3)
		sbavg3=sbsum1/nclump3
		SBclump3.append(sbavg3)
		magclump4=dataclump4[band+'mag_aperture0'+js]
		sbclump4=magclump4+2.5*math.log10(4*pi*aperture[j]**2)
		sbsum4=np.sum(sbclump4)
		sbavg4=sbsum4/nclump4
		SBclump4.append(sbavg4)
	f=plt.figure()
	plt.subplots_adjust(bottom=0.075, top=.97, right=.9)
	ax0=plt.subplot(221)
	plt.scatter(aperture, SBclump1, c='k', marker='^')
	plt.xlabel('Aperture Radius (arcseconds)', fontsize=5)
	plt.ylabel('Averaged Surface Brightness(mag/arcsec^2)', fontsize=5)
	plt.title("Averaged Surface Brightness vs. Aperture Radius in "+band, fontsize=7)
	plt.ylim(min(SBclump1)-1,max(SBclump1)+1)
	ax0.invert_yaxis()
	plt.plot(0,0, label='Zrange is '+range1)
	plt.legend(loc=1,prop={'size':6})
	ax1=plt.subplot(222)
	plt.scatter(aperture, SBclump2, c='k', marker='^')
	plt.xlabel('Aperture Radius (arcseconds)', fontsize=5)
	plt.ylabel('Averaged Surface Brightness(mag/arcsec^2)', fontsize=5)
	plt.title("Averaged Surface Brightness vs. Aperture Radius in "+band, fontsize=7)
	plt.ylim(min(SBclump2)-1,max(SBclump2)+1)
	ax1.invert_yaxis()
	plt.plot(0,0, label='Zrange is '+range2)
	plt.legend(loc=1,prop={'size':6})
	ax2=plt.subplot(223)
	plt.scatter(aperture, SBclump3, c='k', marker='^')
	plt.xlabel('Aperture Radius (arcseconds)', fontsize=5)
	plt.ylabel('Averaged Surface Brightness(mag/arcsec^2)', fontsize=5)
	plt.title("Averaged Surface Brightness vs. Aperture Radius in "+band, fontsize=7)
	plt.ylim(min(SBclump3)-1,max(SBclump3)+1)
	ax2.invert_yaxis()
	plt.plot(0,0, label='Zrange is '+range3)
	plt.legend(loc=1,prop={'size':6})
	ax3=plt.subplot(224)
	plt.scatter(aperture, SBclump4, c='k', marker='^')
	plt.xlabel('Aperture Radius (arcseconds)', fontsize=5)
	plt.ylabel('Averaged Surface Brightness(mag/arcsec^2)', fontsize=5)
	plt.title("Averaged Surface Brightness vs. Aperture Radius in "+band, fontsize=7)
	plt.ylim(min(SBclump4)-1,max(SBclump4)+1)
	plt.plot(0,0, label='Zrange is '+range4)
	plt.legend(loc=1,prop={'size':6})
	ax3.invert_yaxis()
	plt.show()
	f.savefig(outdir+band+'_zclump.pdf')
