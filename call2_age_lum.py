import astropy.table as table 
import numpy as np
from defcuts import *
from defflags import *
from halflight_first import *
from def_get_mags import *
from def_halflight_math import *

from def_ages import *

ty='mean'

stax=True
if stax==False:
	tag=''
else:
	tag='uplim'

txtdist= ''
txtslope=''

outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/2'+ty+tag
doutdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/distribution/2'+ty+tag

indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'
	
DATA=table.Table.read(indir+'med_vespa_LOWZ.fits')


bands=['g', 'r', 'i','z', 'y']
daperture=[1.01,1.51,2.02,3.02,4.03,5.71,8.40,11.8,16.8,23.5]
aperture=[x*0.5 for x in daperture]

def do_cuts(datatab):
	parm=['flags_pixel_saturated_center','flags_pixel_edge','flags_pixel_interpolated_center','flags_pixel_cr_center','flags_pixel_suspect_center', 'flags_pixel_clipped_any','flags_pixel_bad']
	ne=[99.99, 199.99, 0.0]
	mincut=0.1
	maxcut=''
	cutdata=not_cut(datatab, bands, 'mag_forced_cmodel', ne)

	for b in range(0, len(bands)-1):
		newdata=many_flags(cutdata, parm, bands[b])	#flags not in y?
		cutdata=newdata
	bandi=['i']
	
	return newdata
	
DATA=do_cuts(DATA)

def get_agebin_dat(Data, hm):
	print('Running Age Bin')
	runid=Data['RUNID']
	runIDs, count=np.unique(runid, return_counts=True)
	ndata=Data[(runid==5) | (runid==1)]
	runtest=''
	if runtest=='yes':
		test1=Data[runid==1]
		test5=Data[runid==5]
		gal1, count1=np.unique(test1['SPECOBJID'],return_counts=True)
		gal5, count5=np.unique(test5['SPECOBJID'],return_counts=True)
		print('Number of galaxies in run1= ', len(gal1))
		print('Howmany times they repeat: ', count1)
		print('Number of galaxies in run5= ', len(gal5))
		print('Howmany times they repeat: ', count5)
		gal, count=np.unique(ndata['SPECOBJID'],return_counts=True)
		print('Number of galaxies in run==1,5: ', len(gal))
		print('How many times they repeat: ', count)
	newdata, notdat=mass_frac_cut1(ndata, hm, get_opp=True)

	return newdata, notdat
	
hh=0.7
newdata, datanot= get_agebin_dat(DATA, hh)

starts1=newdata['AGESTART']
starts2=datanot['AGESTART']
data1=newdata[starts1==9.04]
data2=datanot[starts2==9.04]

def my_halflight2(dat1):
	loglum, lograd, loglumd= get_ind_lums(dat1, bands, aperture, scale='log')
	
	if stax==True:
		print('hi')
		loglum, lograd, loglumd= upper_rad_cut(loglum, lograd, loglumd, 4, proof=False)
	#print('length of radius array is ', len(lograd))
	
	logr12s, logr412s= get_halflight2(loglum, lograd, 4)
	print('min radius in lin= ', np.min(10**lograd), 'max radius in lin= ', np.max(10**lograd))
	print('min r1/2 is ', np.min(10**logr12s),'max 4r1/2 is ', np.max(10**logr412s))
	
	#mloglum,  mlogdens, mlograd, mlogerr= get_avg_lums(loglum, lograd, loglumd, gr=[0.7,71,11], type=ty, scale='lindata')
	mloglum,  mlogdens, mlograd, mlogerr= get_avg_lums(loglum, lograd, loglumd, gr=[0.7,71,11], type=ty, scale='lindata')
	
	logr12s, logr412s= get_halflight2(loglum, lograd, 4)
	
	logr12, logr412= get_halflight2(mloglum, mlograd, 4)
	
	#for n in range(len(logr12s)):
	#	print(np.round(logr12s[n],3), np.round(logr412s[n],3), 'lograd= ', np.round(lograd[n],3))
	
	Ms, cs, errs= get_slopes1(logr12s, logr412s,lograd, loglumd, error=None, smax=stax)
	M, c, logrcut, logldcut, sterr, errcut =get_slopes1(logr12, logr412, mlograd, mlogdens, error=mlogerr, smax=stax)
	
	cutmlogld = M * logrcut + c
	
	ind=[loglum, loglumd, lograd, logr12s]
	means=[mloglum,mlogdens,mlograd,logr12, mlogerr]
	ind_slope=[Ms, cs, errs]
	mean_slopes=[M, c, logrcut, logldcut, cutmlogld, sterr, errcut]
	#logrcut and logldcut are for lines of best fit
	
	return ind, means, ind_slope, mean_slopes
	
	
def my_graphs(inds1, means1, ind_slope1, mean_slopes1, inds2, means2, ind_slope2, mean_slopes2):
	
	per=[str(hh*100), '%']
	per=''.join(per)
	tag1=['Number of Galaxies= '+str(len(inds1[0])), 'Galaxies w/ Mass Fractions >'+per,'Mass Fractions > '+per]
	tag2=['Number of Galaxies= '+str(len(inds2[0])), 'Galaxies w/ Mass Fractions <'+per,'Mass Fractions < '+per]
	#inds=[lum1, lumd1, rad1, hrad1]
	#means=[mlum1,mdens1,mrad1,mhrad1, merr1]
	#ind_slope=[m1s, c1s, err1s]
	#mean_slopes=[m1, c1, radcut1, dencut1, ynew1,sterr1, errcut1]
	
	def lum_mult_fit(x1, x2, y1, y2, xcut1, xcut2, yfit1, yfit2, sterr1, sterr2 , m1, m2, error1, error2, outdir=''):
		print('Make Scatter Plots')
		import matplotlib.pyplot as plt
		import numpy as np
		import math
		f=plt.figure()
		plt.scatter(x1, y1, color='r', marker='o',label=tag1[1]+' ('+str(len(inds1[0]))+')')
		plt.plot(xcut1, yfit1, color='m', label='(>'+str(per)+') mean slope= '+str(round(m1,2))+' +- '+str(round(sterr1,2)))
		plt.errorbar(x1, y1, yerr=error1, fmt='.',color='r')	

		plt.scatter(x2, y2, color='b', marker='o',label=tag2[1]+' ('+str(len(inds2[0]))+')')
		plt.plot(xcut2, yfit2, color='c', label='(<'+str(per)+') mean slope= ' +str(round(m2,2))+' +- '+str(round(sterr2,2)))
		plt.errorbar(x2, y2, yerr=error2, fmt='.',color='b')

		plt.xlabel('Log Radii (kpc)')
		plt.ylabel('Luminosity Densities (Lsolar/kpc^2)')
		plt.title('Average Luminosity Densities v Radii')
		plt.legend(loc=0,prop={'size':6.0})
		#f.text(0.05, 0.05, txtslope, color='red', weight='bold')
		outdirs=outdir+'lumage.pdf'
		#plt.show()
		f.savefig(outdirs)
		print(outdirs)

	def dist_mean(m1s, m2s, m1, m2, sterr1, sterr2, KS=False):
		import matplotlib.pyplot as plt
		import numpy as np
		import math
		figs=plt.figure()
		bs=np.linspace(-2.0,-1.4,num=15, endpoint=False)
		n1, b1, p1= plt.hist(m1s, bs, color='red', label=tag1[1]+ ' ('+str(len(m1s))+')', alpha=0.65, zorder=2, normed=1)
		n2, b2, p2= plt.hist(m2s,bs, color='blue', label=tag2[1]+ ' ('+str(len(m2s))+')', alpha=0.65,zorder=2,normed=1)
		ts=''
		if KS==True:
			M=m1s+m2s
			import scipy
			D, p=scipy.stats.ks_2samp(m1s,m2s)
			plt.plot(0,0, c='green', marker='*', label='K-S test is '+str(np.round(D,3)))
			plt.xlim(np.min(M),-1.3)
			ts='KS'
	
		#print('Standard Deviation ('+tag1[2]+'): ', str(round(np.std(m1s),2)))
		#print('Standard Deviation ('+tag2[2]+'): ', str(round(np.std(m2s),2)))
		
		plt.axvline(x=m1, color='magenta',label='(>'+str(per)+') mean slope= '+str(round(m1,2))+' +- '+str(round(sterr1,2)), zorder=3)
		plt.axvline(x=m2, color='cyan', label='(<'+str(per)+') mean slope= '+str(round(m2,2))+' +- '+str(round(sterr2,2)), zorder=3)
		plt.xlabel('Slopes', fontsize=10)
		
		plt.legend(loc=0,prop={'size':6.5})
		plt.ylabel('Frequency', fontsize=10)
		plt.title('With '+ty+' Slopes')
		
		
		outdirs=doutdir+ts+'slope_agedist.pdf'
		#figs.text(0.03, 0.03, txtdist, color='red', weight='bold')
		#plt.show()
		figs.savefig(outdirs)
		print(outdirs)
		
	def slopevLmax(m1, m2, L1, L2):
		import matplotlib.pyplot as plt
		N1=len(m1)
		N2=len(m2)
		Lmax1=[np.max(L1[n]) for n in range(N1)]
		Lmax2=[np.max(L2[n]) for n in range(N2)]
		#gives us Lmax
	
		fs=plt.figure()
		plt.scatter(m1, Lmax1, color='red', label='Not Flagged Galaxies')
		plt.scatter(m2, Lmax2, color='blue', label='Flagged Galaxies')
		plt.xlabel('Slopes')
		plt.ylabel('Max Luminosities (Lsolar)')
		plt.title('Max Luminosities v Slopes')
		plt.legend(loc=0,prop={'size':7.0})
		plt.show()
		
	def all_lumprof(lum1s, lum2s, rad1s, rad2s, x1, x2, y1, y2, error1, error2):
		f=plt.figure()
		print(x1)
		print(x2)
		print(y1)
		print(y2)
		for n in range(len(lum1s)):
			plt.plot(rad1s[n], lum1s[n],color='lightgrey', marker='.')
		for n in range(len(lum2s)):
			plt.plot(rad2s[n], lum2s[n],color='lightgrey', marker='.')
		plt.scatter(x1, y1, color='red', marker='o',label='#'+tag1[1]+': '+ str(len(inds1[0])), zorder=3)
		plt.scatter(x2,y2,color='blue', marker='o',label='#'+tag2[1]+': '+str(len(inds2[0])), zorder=3)
		plt.xlabel('Log Radii (kpc)')
		plt.ylabel('Luminosity Densities (Lsolar/kpc^2)')
		plt.title('Average Luminosity Densities v Radii')
		plt.legend(loc=0,prop={'size':6.0})
		#plt.show()
		outdirs=outdir+'allage_lumprof.pdf'
		#plt.show()
		f.savefig(outdirs)
		print(outdirs)
			
	all_lumprof(inds1[1], inds2[1], inds1[2], inds2[2], means1[2], means2[2], means1[1], means2[1],means1[4], means2[4])
	#slopevLmax(ind_slope1[0],ind_slope2[0], inds1[1], inds2[1])
	dist_mean(ind_slope1[0],ind_slope2[0],mean_slopes1[0],mean_slopes2[0],mean_slopes1[5], mean_slopes2[5], KS=False)
	
	lum_mult_fit(means1[2], means2[2], means1[1], means2[1], mean_slopes1[2], mean_slopes2[2], mean_slopes1[4], mean_slopes2[4], mean_slopes1[5], mean_slopes2[5], mean_slopes1[0], mean_slopes2[0],means1[4], means2[4], outdir=outdir)
	
inds1, means1, ind_slope1, mean_slopes1=my_halflight2(data1)
inds2, means2, ind_slope2, mean_slopes2=my_halflight2(data2)		
		
my_graphs(inds1, means1, ind_slope1, mean_slopes1, inds2, means2, ind_slope2, mean_slopes2)

flagtest=''
if flagtest':
	Flag1=['flags_pixel_bright_object_center', 'brobj_cen_flag-', 'No Bright Ojbect Centers', 'Only Bright Object Centers', 'brobj_cen_flag']
	
	Flag2=['flags_pixel_bright_object_any', 'brobj_all_flag-', 'No Bright Ojbects', 'Only Bright Objects', 'brobj_all_flag']
	bandi='i'
	_, flag1,lab= TFflag(bandi,Flag1, data1)
	_,flag2, lab= TFflag(bandi,Flag1, data2)
	_, flag3,lab= TFflag(bandi,Flag2, data1)
	_,flag4, lab= TFflag(bandi,Flag2, data2)
	print('Total Objects in older= ', len(data1))
	print('Bright Object Centers in older= ', len(flag1))
	print('Bright Objects in older= ', len(flag3))
	
	print('Total Objects in younger= ', len(data2))
	print('Bright Object Centers in younger= ', len(flag2))
	print('Bright Objects in younger= ', len(flag4))
	
	
	
	
	
	
	
	
	
	
	
	
	
	
#not in use currently	
def my_halflight(dat1):
	loglum, lograd, loglumd= get_ind_lums(dat1, bands, aperture, scale='log')
	if stax==True:
		loglum, lograd, loglumd= upper_rad_cut(loglum, lograd, loglumd, 4, proof=False)
	#print('length of radius array is ', len(lograd))
	
	mloglum,  mlogdens, mlograd, mlogerr= get_avg_lums(loglum, lograd, loglumd, gr=[1,80,11],type=ty, scale='lindata')
	
	logr12s= get_halflight(loglum, lograd)
	
	logr12= get_halflight(mloglum, mlograd)
	
	Ms, cs, errs= get_slopes(logr12s, lograd, loglumd, error=None, smax=stax)
	M, c, logrcut, logldcut, sterr, errcut =get_slopes(logr12, mlograd, mlogdens, error=mlogerr, smax=stax)
	
	cutmlogld = M * logrcut + c
	
	ind=[loglum, loglumd, lograd, logr12s]
	means=[mloglum,mlogdens,mlograd,logr12, mlogerr]
	ind_slope=[Ms, cs, errs]
	mean_slopes=[M, c, logrcut, logldcut, cutmlogld, sterr, errcut]
	#logrcut and logldcut are for lines of best fit
	
	return ind, means, ind_slope, mean_slopes