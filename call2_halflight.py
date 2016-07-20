print('Finding Half Light Radii and Corresponding Slopes of Stacked Galaxies')

import astropy.table as table 
import numpy as np
from defcuts import *
from defflags import *
from halflight_first import *
from def_get_mags import *
from def_halflight_math import *

bands=['g', 'r', 'i','z', 'y']
daperture=[1.01,1.51,2.02,3.02,4.03,5.71,8.40,11.8,16.8,23.5]
aperture=[x*0.5 for x in daperture]

ty='mean'

stax=True
if stax==False:
	tag=''
else:
	tag='uplim'
txtdist= 'Figure2'
txtslope='Figure1'

outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/2'+ty+tag
doutdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/distribution/2'+ty+tag
Flags=['flags_pixel_bright_object_center', 'brobj_cen_flag-', 'No Bright Ojbect Centers', 'Only Bright Object Centers', 'brobj_cen_flag']

indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'
bigdata = table.Table.read(indir+ 'LOWZ_HSCGAMA15_apmgs+cmodmag.fits')
def do_cuts(datatab):
	parm=['flags_pixel_saturated_center','flags_pixel_edge','flags_pixel_interpolated_center','flags_pixel_cr_center','flags_pixel_suspect_center', 'flags_pixel_clipped_any','flags_pixel_bad']
	ne=[99.99, 199.99, 0.0]
	mincut=0.1
	maxcut=''
	cutdata=not_cut(datatab, bands, 'mag_forced_cmodel', ne)
	for b in range(0, len(bands)-1):
		newdata=many_flags(cutdata, parm, bands[b])	#flags not in y?
		cutdata=newdata
	
	return newdata
def get_TF(data):
	bandi=['i']
	Flag, Not,lab= TFflag(bandi,Flags, data)
	return Flag, Not	
	
newdata=do_cuts(bigdata)

Flagdat, Notdat=get_TF(newdata)

def my_halflight2(dat1):
	loglum, lograd, loglumd= get_ind_lums(dat1, bands, aperture, scale='log')
	
	if stax==True:
		loglum, lograd, loglumd= upper_rad_cut(loglum, lograd, loglumd, 4, proof=False)
	#print('length of radius array is ', len(lograd))
	
	mloglum,  mlogdens, mlograd, mlogerr= get_avg_lums(loglum, lograd, loglumd, type=ty, scale='lindata')
	
	logr12s= get_halflight(loglum, lograd)
	
	logr12= get_halflight(mloglum, mlograd)
	
	Ms, cs, errs= get_slopes(logr12s, lograd, loglumd, error=None, smax=stax)
	M, c, logrcut, logldcut, sterr, errcut =get_slopes(logr12, mlograd, mlogdens, error=mlogerr, smax=stax)
	print(sterr)
	
	cutmlogld = M * logrcut + c
	
	ind=[loglum, loglumd, lograd, logr12s]
	means=[mloglum,mlogdens,mlograd,logr12, mlogerr]
	ind_slope=[Ms, cs, errs]
	mean_slopes=[M, c, logrcut, logldcut, cutmlogld, sterr, errcut]
	#logrcut and logldcut are for lines of best fit
	
	return ind, means, ind_slope, mean_slopes
	
inds1, means1, ind_slope1, mean_slopes1=my_halflight2(Flagdat)
inds2, means2, ind_slope2, mean_slopes2=my_halflight2(Notdat)

print('mean radii= ', means1[2], means2[2])
def my_graphs(inds1, means1, ind_slope1, mean_slopes1, inds2, means2, ind_slope2, mean_slopes2):
	import matplotlib.pyplot as plt
	import numpy as np
	import math
	#ind=[loglum, loglumd, lograd, logr12s]
	#means=[mloglum,mlogdens,lograd,logr12, mlogerr]
	#ind_slope=[Ms, cs, errs]
	#mean_slopes=[M, c, logrcut, logldcut, cutmlogld, sterr, errcut]
	
	def lum_mult_fit(x1, x2, y1, y2, xcut1, xcut2, yfit1, yfit2, sterr1, sterr2 , m1, m2, error1, error2, outdir=''):
		print('Make Scatter Plots')
		f=plt.figure()
		plt.scatter(x1, y1, color='r', marker='o',label='Not Flagged Galaxies')
		plt.plot(xcut1, yfit1, color='m', label='Fitted Not Flagged Galaxies: slope= '+str(np.round(m1,2))+' +- '+str(np.round(sterr1,2)))
		plt.errorbar(x1, y1, yerr=error1, fmt='.',color='r')	

		plt.scatter(x2, y2, color='b', marker='o',label='Flagged Galaxies')
		plt.plot(xcut2, yfit2, color='c', label='Fitted Flagged Galaxies: slope= '+str(np.round(m2,2))+' +- '+str(np.round(sterr2,2)))
		plt.errorbar(x2, y2, yerr=error2, fmt='.',color='b')

		plt.xlabel('Log Radii (kpc)')
		plt.ylabel('Luminosity Densities (Lsolar/kpc^2)')
		plt.title('Average Luminosity Densities v Radii')
		#plt.xlim(math.log10(1), math.log10(80))
		#plt.ylim(6,8.6)
		plt.legend(loc=0,prop={'size':6.0})
		f.text(0.05, 0.05, txtslope, color='red', weight='bold')
		outdirs=outdir+tag+'TF.pdf'
		#plt.show()
		f.savefig(outdirs)
		print(outdirs)

	def dist_mean(m1s, m2s, m1, m2, sterr1, sterr2, KS=False):

		figs=plt.figure()
		bs=np.linspace(-2.0,-1.4,num=15, endpoint=False)
		n1, b1, p1= plt.hist(m1s, bs, color='red', label='Not Flagged Galaxies ('+str(len(m1s))+')', alpha=0.8)
		n2, b2, p2= plt.hist(m2s,bs, color='blue', label='Flagged Galaxies ('+str(len(m2s))+')', alpha=0.8)
		
		ts=''
		if KS==True:
			M=m1s+m2s
			import scipy
			D, p=scipy.stats.ks_2samp(m1s,m2s)
			plt.plot(0,0, c='green', marker='*', label='K-S test is '+str(D))
			plt.xlim(np.min(M),-1.4)
			ts='KS'
		
		print('Standard Deviation (Not Flagged): ', str(np.std(m1s)))
		print('Standard Deviation (Flagged): ', str(np.std(m2s)))
		
		plt.axvline(x=m1, color='magenta', label='Not Flagged Galaxies: slope= '+str(np.round(m1,2))+'  +- ' +str(np.round(sterr1,2)))
		plt.axvline(x=m2, color='cyan', label='Flagged Galaxies: slope= '+str(np.round(m2,2))+' +- '+str(np.round(sterr2,2)))
		plt.xlabel('Slopes', fontsize=10)
		plt.legend(loc=0,prop={'size':6.5})
		plt.ylabel('Frequency', fontsize=10)
		plt.title('With '+ty+' Slopes')
	
		outdirs=doutdir+'slopedist.pdf'
		#figs.text(0.03, 0.03, txtdist, color='red', weight='bold')
		#plt.show()
		figs.savefig(outdirs)
		print(outdirs)
		
	def all_lumprof(lum1s, lum2s, rad1s, rad2s, mrad1, mrad2, mden1, mden2, error1, error2):
		f=plt.figure()
		#print(len(mrad1)) #these are the mean radii
		#print(len(mrad2))
		#print(len(mden1))
		#print(len(mden2))
		for n in range(len(lum1s)):
			plt.plot(rad1s[n], lum1s[n],color='lightgrey', marker='.')
		for n in range(len(lum2s)):
			plt.plot(rad2s[n], lum2s[n],color='lightgrey', marker='.')
		plt.scatter(mrad1, mden1, color='red', marker='o',label='Not Flagged Galaxies ('+str(len(lum1s))+')', zorder=3)
		plt.scatter(mrad2,mden2,color='blue', marker='o',label='Flagged Galaxies ('+str(len(lum1s))+')', zorder=3)
		plt.xlabel('Log Radii (kpc)')
		plt.ylabel('Luminosity Densities (Lsolar/kpc^2)')
		plt.title('Average Luminosity Densities v Radii')
		plt.legend(loc=0,prop={'size':6.0})
		outdirs=outdir+tag+'all_lumprof.pdf'
		#plt.show()
		f.savefig(outdirs)
		print(outdirs)
		
	dist_mean(ind_slope1[0],ind_slope2[0],mean_slopes1[0],mean_slopes2[0],mean_slopes1[5], mean_slopes2[5], KS=False)
	
	all_lumprof(inds1[1], inds2[1], inds1[2], inds2[2], means1[2], means2[2], means1[1], means2[1],means1[4], means2[4])
	
	lum_mult_fit(means1[2], means2[2], means1[1], means2[1], mean_slopes1[2], mean_slopes2[2], mean_slopes1[4], mean_slopes2[4], mean_slopes1[5], mean_slopes2[5], mean_slopes1[0], mean_slopes2[0],means1[4], means2[4], outdir=outdir)
		
my_graphs(inds1, means1, ind_slope1, mean_slopes1, inds2, means2, ind_slope2, mean_slopes2)