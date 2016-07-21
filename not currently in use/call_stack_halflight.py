print('Finding Half Light Radii and Corresponding Slopes of Stacked Galaxies')

import astropy.table as table 
import numpy as np
from defcuts import *
from def_get_mags import *
from def_lump_prof import *
from my_def_plots import *
from defflags import *
from defclump import * 

indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'

datatab = table.Table.read(indir+ 'LOWZ_HSCGAMA15_apmgs+cmodmag.fits')

bands=['g', 'r', 'i','z', 'y']

parm=['flags_pixel_saturated_center','flags_pixel_edge','flags_pixel_interpolated_center','flags_pixel_cr_center','flags_pixel_suspect_center', 'flags_pixel_clipped_any','flags_pixel_bad']

Flags=['flags_pixel_bright_object_center', 'brobj_cen_flag-', 'No Bright Ojbect Centers', 'Only Bright Object Centers', 'brobj_cen_flag']

daperture=[1.01,1.51,2.02,3.02,4.03,5.71,8.40,11.8,16.8,23.5]

aperture=[x*0.5 for x in daperture]

ne=[99.99, 199.99, 0.0]

mincut=0.1
maxcut=''
cutdata=not_cut(datatab, bands, 'mag_forced_cmodel', ne)

#get rid of flagged galaxies
for b in range(0, len(bands)-1):
	newdata=many_flags(cutdata, parm, bands[b])	#flags not in y?
	cutdata=newdata

bandi=['i']

outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'

Flagdat, Notdat,lab= TFflag(bandi,Flags, newdata)

objid=Flagdat['object_id_1']

def my_halflight(Flagdat, Notdat):
	meandensF, radavgF, lumdenerrF, halfradF = get_mean_halflight(Flagdat, bands, aperture, scale='log')

	mF, cF, xcutF, ycutF, st_errF, eF = halflight_meanslope(halfradF, radavgF, meandensF, lumdenerrF)


	meandensN, radavgN, lumdenerrN, halfradN = get_mean_halflight(Notdat, bands, aperture, scale='log')

	mN, cN, xcutN, ycutN, st_errN, eN = halflight_meanslope(halfradN, radavgN, meandensN, lumdenerrN)

	ynewF = mF * xcutF + cF

	ynewN = mN * xcutN + cN
	
	def lum_mult_fit(x1, x2, y1, y2, xcut1, xcut2, yfit1, yfit2, error1, error2, outdir=''):
		print('Make Scatter Plots')
		import matplotlib.pyplot as plt
		import numpy as np
		import math
		f=plt.figure()
		plt.scatter(x1, y1, color='r', marker='o',label='Not Flagged Galaxies')
		plt.plot(xcut1, yfit1, color='m', label='Fitted Not Flagged Galaxies: slope= '+str(mF)+' +- '+str(st_errF))
		plt.errorbar(x1, y1, yerr=error1, fmt='.',color='r')	
	
		plt.scatter(x2, y2, color='b', marker='o',label='Flagged Galaxies')
		plt.plot(xcut2, yfit2, color='c', label='Fitted Flagged Galaxies: slope= '+str(mN)+' +- '+str(st_errN))
		plt.errorbar(x2, y2, yerr=error2, fmt='.',color='b')
	
		plt.xlabel('Log Radii (kpc)')
		plt.ylabel('Luminosity Densities (Lsolar/kpc^2)')
		plt.title('Average Luminosity Densities v Radii')
		plt.xlim(math.log10(1), math.log10(80))
		plt.ylim(6,8.6)
		#plt.ylim(np.min(x1)-2*error1, np.max(x1)+2*error1)
		plt.legend(loc=0,prop={'size':6.0})

		#a = plt.axes([.18, .18, .3, .3], axisbg='y')
		#plt.scatter(xcut1, ycutF, color='r', marker='.')
		#plt.plot(xcut1, yfit1, color='m')
		#plt.errorbar(xcut1, ycutF, yerr=eF, fmt='.',color='r')
		#plt.scatter(xcut2, ycutN, color='b', marker='.')
		#plt.plot(xcut2, yfit2, color='c')
		#plt.errorbar(xcut2, ycutN, yerr=eN, fmt='.',color='b')
		#plt.plot(0,0, color='m',label='Standard Error of Not Flagged Slope: '+str(st_errF))
		#plt.plot(0,0, color='c',label='Standard Error of Flagged Slope: '+str(st_errN))
		#plt.xlim(0.7, 1.9)
		#plt.ylim(6.1,7.8)
		#plt.title('Radii > R1/2', fontsize=6.0)
		#plt.tick_params(axis='both', which='major', labelsize=5.0)
		#plt.legend(loc=0,prop={'size':3.0})
		
		outdirs=outdir+'TF.pdf'
		f.savefig(outdirs)
		print(outdirs)
		
	def dist_mean(Flagdat, Notdat):
		import matplotlib.pyplot as plt
		import numpy as np
		import math
		figs=plt.figure()
		halfrad, rad, lumdens=get_halflight(Flagdat, bands, aperture, scale='log')
		halfradN, radN, lumdensN=get_halflight(Notdat, bands, aperture, scale='log')
		m1=halflight_slopes(halfrad, rad, lumdens, objid, plots='no', outdir=outdir)
		m2=halflight_slopes(halfradN, radN, lumdensN, objid, plots='no', outdir=outdir)
		bs=np.linspace(-2.0,-1.4,num=15, endpoint=False)
		plt.hist(m1, bs, color='red', label='Not Flagged Galaxies', alpha=0.8)
		plt.hist(m2,bs, color='blue', label='Flagged Galaxies', alpha=0.8)
		plt.axvline(x=mF, color='magenta', label='Not Flagged Galaxies: slope= '+str(mF)+'  +- ' +str(st_errF))
		plt.axvline(x=mN, color='cyan', label='Flagged Galaxies: slope= '+str(mN)+' +- '+str(st_errN))
		plt.xlabel('Slopes', fontsize=10)
		plt.legend(loc=0,prop={'size':7.0})
		plt.ylabel('Frequency', fontsize=10)
		plt.show()
	
	dist_mean(Flagdat, Notdat)
	
	lum_mult_fit(radavgF, radavgN, meandensF, meandensN, xcutF, xcutN, ynewF, ynewN, lumdenerrF, lumdenerrN, outdir=outdir)
	
	
my_halflight(Flagdat, Notdat)
