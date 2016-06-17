# this will plot the magnitudes for saturated and unsaturated galaxies in i

indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'

outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/single_plot/'

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.text as txt
from astropy.io import fits
import astropy.table as table 
import matplotlib.patches as mpatches
from sympy import *
import os, sys

data = table.Table.read(indir+ 'LOWZ_HSCGAMA15_apmgs.fits')


aperture=[3.0*0.168,4.5*0.168,6.0*0.168,9.0*0.168,12.0*0.168,17.0*0.168,25.0*0.168,35.0*0.168,50.0*0.168,70.0*0.168]
pi=math.pi
Naps=len(aperture)

print(len(data))
band='i'

dataflags = data[(data[band+'flags_pixel_saturated_center']==False)&(data[band+'flags_pixel_edge']==False)&(data[band+'flags_pixel_cr_center']==False)&(data[band+'flags_pixel_bad']==False)&(data[band+'flags_pixel_suspect_center']==False)&(data[band+'flags_pixel_clipped_any']==False)] 
dataflag=dataflags[dataflags['imag_aperture00']>1.0]
#datanot= data[(data[band+'flags_pixel_saturated_center']==True)&(data[band+'flags_pixel_edge']==True)&(data[band+'flags_pixel_cr_center']==True)&(data[band+'flags_pixel_bad']==True)&(data[band+'flags_pixel_suspect_center']==True)&(data[band+'flags_pixel_clipped_any']==True)]

datanot=data[dataflag==False]

numflag=len(dataflag)
strnumflag=str(numflag)
numnot=len(datanot)
strnumnot=str(numnot)
print(numflag, numnot)

#for aperture00


fig=plt.figure()
#for n in range(0, numflag):	#will print all of the galaxies in this range
#    flagmag=dataflag['imag_aperture00']
    #print('flagmag is ', flagmag)
#    plt.plot(aperture[0], flagmag, c='r', marker='^')
#for n in range(0, numnot):
#    flagnot=datanot['imag_aperture00'][n]
    #print('flagnot is ',flagnot)
#    plt.plot(aperture[0], flagnot, c='b', marker='*')

nbins=8
flagmag=dataflag['imag_aperture00']
flagnot=datanot['imag_aperture00']
plt.hist(flagmag,color='red',alpha=0.8)
plt.hist(flagnot,color='blue',alpha=0.8)
plt.xlabel('Aperture Magnitude', fontsize=10)
plt.ylabel('Frequency', fontsize=10)
plt.suptitle("Frequency of Different Aperture Magnitudes", fontsize=15)
#red_patch = mpatches.Patch(color='red', label='UnSaturated Galaxies' +'('+strnumflag+')')
plt.plot(0,0,label='No Saturated, No Edge, No Cr Cen, No Suspect, No Bad Pixels, No Clipped' +'('+strnumflag+')', c='r', marker='^')
#plt.plot(0,0,label='Unsaturated Range('+minf+' to '+maxf+')', c='r')
#blue_patch = mpatches.Patch(color='blue', label='Saturated Galaxies' +'('+strnumnot+')')
plt.plot(0,0,label='Saturated, Edge, Cr Cen, Suspect, Bad Pixels, Clipped' +'('+strnumnot+')', c='b', marker='*')
#plt.plot(0,0,label='Saturated Range('+minn+' to '+maxn+')', c='b')
#plt.legend(handles=[red_patch, blue_patch], loc=7,prop={'size':5})
plt.legend(loc=2,prop={'size':5})
#plt.gca().invert_yaxis()
plt.show()

