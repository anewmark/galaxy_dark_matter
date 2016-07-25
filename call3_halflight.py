print('doing this linearly this time!')
import astropy.table as table 
import numpy as np
from defcuts import *
from defflags import *
import matplotlib.pyplot as plt
from halflight3_first import *
from def_get_mags import *
from def_halflight_math import *

bands=['g', 'r', 'i','z', 'y']
daperture=[1.01,1.51,2.02,3.02,4.03,5.71,8.40,11.8,16.8,23.5]
aperture=[x*0.5 for x in daperture]

ty='mean'
#ty='med'

stax=True
if stax==False:
	tag=''
else:
	tag='uplim'
txtdist= 'Figure2'
txtslope='Figure1'

outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/3'+ty+tag
doutdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/distribution/3'+ty+tag
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
Notdat=Notdat[Notdat['Z']>0.2]   ## WANT X ZCUT FOR NOTDAT ONLY

def my_halflights(dat1):
	lum, rad, ld= get_ind_lums(dat1, bands, aperture, scale='linear')
	#print(lum[0], rad[0], lumd[0]) #this confirms it is all linear
	if stax==True:
		print('stax is true and doing upper cut now')
		lum, rad, ld= upper_rad_cut(lum, rad, ld, 4, proof=False)
	
	#print('min radius in lin= ', np.min(rad), 'max radius in lin= ', np.max(rad))
	binrange=[2.2,65,10]
	#Flagdat #points in each bin: [592 733 655 717 732 640 667 681 697 610]
	#Notdat #points in each bin:[177 163 187 183 169 165 170 169 185  94]
	mlum,  mdens, mrad, mlogerr= get_avg_lums(lum, rad,ld, gr=binrange, type=ty)
	r12s, r412s= get_halflight2(lum, rad, mult=4)
	r12, r412= get_halflight2(mlum, mrad, mult=4)
	
	halftest=''
	if halftest:
		import matplotlib.pyplot as plt
		plt.hist(r12s, 10,color='green', alpha=.8, label= 'r1/2')
		plt.hist(r412s,10, color='blue', alpha=.8, label= '4r1/2')
		plt.axvline(r12, color='magenta', label='stacked r1/2')
		plt.axvline(r412, color='cyan', label='stacked 4r1/2')
		plt.xlabel('Radii (kpc)', fontsize=10)
		plt.legend(loc=0,prop={'size':6.5})
		plt.show()
		
	#everything at this point comes out as  logs	
	Ms, cs, errs= get_slopes1(r12s, r412s,rad, ld, error=None, scale='log', smax=stax)
	M, c, logrcut, logldcut, sterr, errcut =get_slopes1(r12, r412, mrad, mdens, error=mlogerr,  scale='log',smax=stax)
	
	logldfit = M * logrcut + c #in log
	
	print('the lum dens cut: ', logldcut)
	print('the lum dens fit: ', logldfit)
	
	slopetest='' #distribution of slopes
	if slopetest:
		import matplotlib.pyplot as plt
		plt.hist(Ms, 10,color='green', alpha=.8, label= 'individual slopes')
		plt.axvline(M, color='blue', label='stacked slope')
		plt.xlabel('Slope Distribution', fontsize=10)
		plt.legend(loc=0,prop={'size':6.5})
		plt.show()
		
	goodgod='' #how slope fits to stacked image
	if goodgod:
		import matplotlib.pyplot as plt
		plt.scatter(np.log10(mrad), np.log10(mdens), color='red', marker='o')
		plt.xlabel('Log Radii (kpc)', fontsize=10)
		plt.ylabel('Log Luminosity Density')
		plt.plot(logrcut, logldfit, color='magenta')
		plt.show()
		
	#error is already in log
	print('min: ', np.min(Ms), 'max: ', np.max(Ms), 'mean: ', np.mean(Ms))
	print('stacked m: ', M)
	
	ind=[np.log10(lum), np.log10(ld), np.log10(rad), np.log10(r12s)]
	means=[np.log10(mlum),np.log10(mdens),np.log10(mrad),np.log10(r12), mlogerr]
	ind_slope=[Ms, cs, errs]
	mean_slopes=[M, c, logrcut, logldcut, logldfit, sterr, errcut]
	return ind, means, ind_slope, mean_slopes
	
inds1, means1, ind_slope1, mean_slopes1=my_halflights(Flagdat)
inds2, means2, ind_slope2, mean_slopes2=my_halflights(Notdat)

def z_distr(data1, data2):
	plt.hist(data1, color='red',alpha=.8, label='Not Flagged')
	plt.hist(data2, color='blue', alpha=0.8, label='Flagged')
	plt.xlabel('Redshift Distribution')
	plt.legend()
	plt.show()

#z_distr(Flagdat['Z'], Notdat['Z'])

def my_graphs(inds1, means1, ind_slope1, mean_slopes1, inds2, means2, ind_slope2, mean_slopes2):
	import matplotlib.pyplot as plt
	import numpy as np
	import math
	#ind=[loglum, loglumd, lograd, logr12s]
	#means=[mloglum,mlogdens,lograd,logr12, mlogerr]
	#ind_slope=[Ms, cs, errs]
	#mean_slopes=[M, c, logrcut, logldcut, cutmlogld, sterr, errcut]
	
	def lum_mult_fit(lum1s, lum2s, rad1s, rad2s,x1, x2, y1, y2, xcut1, xcut2, yfit1, yfit2, sterr1, sterr2 , m1, m2, error1, error2, outdir=''):
		print('Make Scatter Plots')
		import matplotlib.pyplot as plt
		import numpy as np
		import math
		f=plt.figure()
		for n in range(len(lum1s)):
			plt.plot(rad1s[n], lum1s[n],color='lightgrey', marker='.')
		for n in range(len(lum2s)):
			plt.plot(rad2s[n], lum2s[n],color='lightgrey', marker='.')
		plt.scatter(x1, y1, color='r', marker='o',label='Not Flagged Galaxies ('+str(len(rad1s))+')')
		plt.plot(xcut1, yfit1, color='m', label='Not Flagged mean slope= '+str(round(m1,3))+' +- '+str(round(sterr1,3)))
		plt.errorbar(x1, y1, yerr=error1, fmt='.',color='r', zorder=4)	

		plt.scatter(x2, y2, color='b', marker='o',label='Flagged Galaxies ('+str(len(rad2s))+')')
		plt.plot(xcut2, yfit2, color='c', label='Flagged mean slope= ' +str(round(m2,3))+' +- '+str(round(sterr2,3)))
		plt.errorbar(x2, y2, yerr=error2, fmt='.',color='b', zorder=4)

		plt.xlabel('Log Radii (kpc)')
		plt.ylabel('Luminosity Densities (Lsolar/kpc^2)')
		plt.title('Average Luminosity Densities v Radii')
		#plt.xlim(math.log10(1), math.log10(80))
		#plt.ylim(6,8.6)
		plt.legend(loc=0,prop={'size':6.0})
		#f.text(0.05, 0.05, txtslope, color='red', weight='bold')
		outdirs=outdir+'TF.pdf'
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
		plt.scatter(mrad2,mden2,color='blue', marker='o',label='Flagged Galaxies ('+str(len(lum2s))+')', zorder=3)
		plt.xlabel('Log Radii (kpc)')
		plt.ylabel('Luminosity Densities (Lsolar/kpc^2)')
		plt.title('Average Luminosity Densities v Radii')
		plt.legend(loc=0,prop={'size':6.0})
		outdirs=outdir+'all_lumprof.pdf'
		#plt.show()
		f.savefig(outdirs)
		print(outdirs)
		
	def just_one(lums,  rads, mlum,  mrad):
		f=plt.figure()
		for n in range(len(lums)):
			plt.plot(rads[n], lums[n],color='lightgrey', marker='.')
		plt.scatter(mrad,mlum,color='blue', marker='o',label='Flagged Galaxies ('+str(len(lums))+')', zorder=3)
		plt.xlabel('Log Radii (kpc)')
		plt.ylabel('Log Luminosities (Lsolar)')
		plt.title('Average Luminosity Densities v Radii')
		plt.legend(loc=0,prop={'size':6.0})
		plt.show()
	
	dist_mean(ind_slope1[0],ind_slope2[0],mean_slopes1[0],mean_slopes2[0],mean_slopes1[5], mean_slopes2[5], KS=False)
	
	#all_lumprof(inds1[1], inds2[1], inds1[2], inds2[2], means1[2], means2[2], means1[1], means2[1],means1[4], means2[4])
	
	
	lum_mult_fit(inds1[1], inds2[1], inds1[2], inds2[2],means1[2], means2[2], means1[1], means2[1], mean_slopes1[2], mean_slopes2[2], mean_slopes1[4], mean_slopes2[4], mean_slopes1[5], mean_slopes2[5], mean_slopes1[0], mean_slopes2[0],means1[4], means2[4], outdir=outdir)
	#just_one(inds2[0],inds2[2], means2[0], means2[2])
		
my_graphs(inds1, means1, ind_slope1, mean_slopes1, inds2, means2, ind_slope2, mean_slopes2)