print('doing this linearly this time (ages)!')
import astropy.table as table 
import numpy as np
from defcuts import *
from defflags import *
import matplotlib.pyplot as plt
from halflight3_first import *
from def_get_mags import *
from def_halflight_math import *
from def_ages import mass_frac_cut1
from my_style import get_presentation
get_presentation()
bands=['g', 'r', 'i','z', 'y']
daperture=[1.01,1.51,2.02,3.02,4.03,5.71,8.40,11.8,16.8,23.5]
aperture=[x*0.5 for x in daperture]

types=['mean', 'med']
ty=types[0]

comps=['oy', 'oF', 'yF']
Toy=comps[0]

stax=True
if stax==False:
	tag=''
else:
	tag='uplim'
txtdist= 'Figure2'
txtslope='Figure1'

my_band='g'
maxin=None
if my_band=='i':
	inner=2.9 #do 3 for cut
	if inner==1:
		bnd=''
	else:
		bnd='+'
if my_band=='g':
	inner=2.3	 #do 3 for cut
	if inner==1:
		bnd='g'
	else:
		bnd='g+'

Flags=['flags_pixel_bright_object_center', 'brobj_cen_flag-', 'No Bright Ojbect Centers', 'Only Bright Object Centers', 'brobj_cen_flag']

indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'
bigdata = table.Table.read(indir+ 'med_vespa_LOWZ.fits')
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
	
DATA=do_cuts(bigdata)
	
print(np.min(DATA['Z_2']))
DATA=DATA[DATA['Z_2']>0.2]
def get_TF(data):
	bandi=my_band
	Flag, Not,lab= TFflag(bandi,Flags, data)
	return Flag, Not
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
	
hh=0.60
per=[str(np.round(hh*100,2)), '$\%$']
per=''.join(per)
odata, ydata= get_agebin_dat(DATA, hh)

starts1=odata['AGESTART']
starts2=ydata['AGESTART']

dataold=odata[starts1==9.04]
datayoung=ydata[starts2==9.04]

Nold, Fold=get_TF(dataold)
Nyoung, Fyoung=get_TF(datayoung)

def my_halflights(dat1,binrange):
	lum, rad, ld= get_ind_lums(dat1, bands, aperture, my_band=my_band, scale='linear')
	#print(lum[0], rad[0], lumd[0]) #this confirms it is all linear
	x=6.0
	if stax==True:
		print('stax is true and doing upper cut now')
		#lum, rad, ld= upper_rad_cut(lum, rad, ld, x, proof=False)
		
	mlum,  mdens, mrad, mlogerr= get_avg_lums(lum, rad,ld, gr=binrange, type=ty)
	r12s, r412s= get_halflight2(lum, rad, mult=x, inn=inner, maxout=np.max(mrad), maxin=maxin)
	r12, r412= get_halflight2(mlum, mrad, mult=x, inn=inner, maxout=np.max(mrad), maxin=maxin)
	
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

def my_graphs(inds1, means1, ind_slope1, mean_slopes1, inds2, means2, ind_slope2, mean_slopes2, subtitle, T, tag1, tag2):
	outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/3'+bnd+T+ty+tag
	doutdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/distribution/3'+bnd+T+ty+tag

	#inds=[lum1, lumd1, rad1, hrad1]
	#means=[mlum1,mdens1,mrad1,mhrad1, merr1]
	#ind_slope=[m1s, c1s, err1s]
	#mean_slopes=[m1, c1, radcut1, dencut1, ynew1,sterr1, errcut1]
	
	def lum_mult_fit(lum1s, lum2s, rad1s, rad2s,x1, x2, y1, y2, xcut1, xcut2, yfit1, yfit2, sterr1, sterr2 , m1, m2, error1, error2, outdir=''):
		print('Make Scatter Plots')
		import matplotlib.pyplot as plt
		plt.rc('text', usetex=True)
		import numpy as np
		import math

		fig=plt.figure()
		ax1, ax2 = plt.subplots(sharex=True)
		plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0) 
		left1, bottom1, width1, height1 = [0.15, 0.3, 0.8, 0.6]
		left2, bottom2, width2, height2 = [0.15, 0.1, 0.8, 0.15]
		ax1 = fig.add_axes([left1, bottom1, width1, height1])
		ax2 = fig.add_axes([left2, bottom2, width2, height2])
		for n in range(len(lum1s)):
			ax1.plot(rad1s[n], lum1s[n],color='lightgrey')
		for n in range(len(lum2s)):
			ax1.plot(rad2s[n], lum2s[n],color='lightgrey')
		ax1.scatter(x1, y1, color='r', marker='o',label=tag1[1]+' ('+str(len(inds1[0]))+')')
		ax1.plot(xcut1, yfit1, color='m', label=tag1[2]+' Stacked Slope')
		ax1.errorbar(x1, y1, yerr=error1, fmt='.',color='r', zorder=4)	

		ax1.scatter(x2, y2, color='b', marker='o',label=tag2[1]+' ('+str(len(inds2[0]))+')')
		ax1.plot(xcut2, yfit2, color='c', label=tag2[2]+' Stacked Slope')
		ax1.errorbar(x2, y2, yerr=error2, fmt='.',color='b', zorder=4)

		#ax1.set_xlabel('Log Radii (kpc)')
		ax1.set_ylabel(r'Log $\rho_L$ (L$_\odot$/kpc$^2$)')
		ax1.set_xlim(0,2.5)
		ax1.legend(loc=0, prop={'size':14.0})
		
		errlum=np.sqrt(error1**2+error2**2)
		difflum=y2-y1
		ax2.scatter(x2, difflum, c='k', marker='o', zorder=2)
		ax2.errorbar(x2, difflum, yerr=errlum, fmt='None')
		ax2.axhline(y=0.0,  color = 'g', linestyle='--')
		ax2.set_ylabel(r'$\Delta\rho_L$')
		#lumrange=[np.round(np.min(difflum),3), 0.0, np.round(np.max(difflum),3)] 
		#ax2.set_yticklabels(lumrange)
		ax2.set_xlabel('Log Radii (kpc)')
		ax2.set_xlim(0,2.5)
		ax2.set_ylim(-0.1, 0.1)
		
		ax1.set_title(subtitle, fontsize=18)
		
		plt.setp(ax1.get_xticklabels(), visible=False)
		plt.setp(ax2.get_yticklabels([-1]), visible=True, fontsize=17)
		plt.subplots_adjust(hspace=.0)
		outdirs=outdir+'lumage.pdf'
		fig.savefig(outdirs)
		print('Slope range: ', 10**np.min(xcut2), 10**np.max(xcut2))
		print(outdirs)	

	def dist_mean(m1s, m2s, m1, m2, sterr1, sterr2, KS=False):
		import matplotlib.pyplot as plt
		import numpy as np
		import math
		figs=plt.figure() 
		bs=np.linspace(-1.9,-1.4,num=20, endpoint=False)
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
		M=np.concatenate((m1s, m2s))
		plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0) 
		plt.axvline(x=m1, color='magenta',label=tag1[2]+' Stacked Slope', zorder=3)
		#plt.plot(0,0, color='magenta', label=tag1[2]+'Median Slope = '+str(np.round(np.median(m1s),3)))
		plt.axvline(x=m2, color='cyan', label=tag2[2]+' Stacked Slope', zorder=3)
		#plt.plot(0,0, color='cyan', label=tag2[2]+'Median Slope = '+str(np.round(np.median(m2s),3)))
		plt.xlabel('Slopes')
		#plt.xlim(-1.9, -1.4)
		plt.legend(loc=2, prop={'size':14.0})
		plt.ylabel('Frequency')
		plt.xlim(np.min(M), np.max(M))
		#plt.suptitle('With '+ty+' Slopes', fontsize=16)
		plt.suptitle(subtitle, fontsize=18)
		
		outdirs=doutdir+ts+'slope_agedist.pdf'
		
		print('NF: ', np.round(m1,3),' $\pm$ ',np.round(sterr1,3))
		print('F: ', np.round(m2,3),' $\pm$ ',np.round(sterr2,3))
		print('NF median: ', np.median(m1s), 'F median: ', np.median(m2s))
		figs.tight_layout()
		figs.savefig(outdirs)
		print(outdirs)
	
	dist_mean(ind_slope1[0],ind_slope2[0],mean_slopes1[0],mean_slopes2[0],mean_slopes1[5], mean_slopes2[5], KS=False)
	
	lum_mult_fit(inds1[1], inds2[1], inds1[2], inds2[2],means1[2], means2[2], means1[1], means2[1], mean_slopes1[2], mean_slopes2[2], mean_slopes1[4], mean_slopes2[4], mean_slopes1[5], mean_slopes2[5], mean_slopes1[0], mean_slopes2[0],means1[4], means2[4], outdir=outdir)

	
if Toy=='oy':
	if ty=='med':
		binrange=[2,80,20] #for med
	elif ty=='mean':
		binrange=[2,80,20] #for mean
	inds1, means1, ind_slope1, mean_slopes1=my_halflights(dataold, binrange)
	inds2, means2, ind_slope2, mean_slopes2=my_halflights(datayoung, binrange)
	sub=['Band ',my_band,': Old ($>$60$\%$) vs. Young ($<$60$\%$)'] 
	sub=''.join(sub)	
	t1=['# Older LRGs= '+str(len(inds1[0])), 'Older','Older']
	t2=['# Younger LRGs= '+str(len(inds2[0])), 'Younger','Younger']
	
	
if Toy=='oF':
	if ty=='med':
		binrange=[2,80,20] #for med
	elif ty=='mean':
		binrange=[2,80,20] #for mean
	inds1, means1, ind_slope1, mean_slopes1=my_halflights(Nold, binrange)
	inds2, means2, ind_slope2, mean_slopes2=my_halflights(Fold, binrange)
	sub='Older: Bright Center Objects Flag'			
	t1=['# of LRGs= '+str(len(inds1[0])), 'Not Bright','Not Bright']
	t2=['# of LRGs= '+str(len(inds2[0])), 'Bright','Bright']
	
	
if Toy=='yF':
	if ty=='med':
		binrange=[2,80,20] #for med
	elif ty=='mean':
		binrange=[2,80,20] #for mean
	inds1, means1, ind_slope1, mean_slopes1=my_halflights(Nyoung, binrange)
	inds2, means2, ind_slope2, mean_slopes2=my_halflights(Fyoung, binrange)
	sub='Younger: Bright Center Objects Flag'		
	t1=['# of LRGs= '+str(len(inds1[0])), 'Not Bright','Not Bright']
	t2=['# of LRGs= '+str(len(inds2[0])), 'Bright','Bright']
	
	
my_graphs(inds1, means1, ind_slope1, mean_slopes1, inds2, means2, ind_slope2, mean_slopes2, sub, Toy, t1, t2)

flagtest=''
if flagtest:
	Flag1=['flags_pixel_bright_object_center', 'brobj_cen_flag-', 'No Bright Ojbect Centers', 'Only Bright Object Centers', 'brobj_cen_flag']
	
	Flag2=['flags_pixel_bright_object_any', 'brobj_all_flag-', 'No Bright Ojbects', 'Only Bright Objects', 'brobj_all_flag']
	bandi='i'
	_, flag1,lab= TFflag(bandi,Flag1, dataold)
	_,flag2, lab= TFflag(bandi,Flag1, datayoung)
	_, flag3,lab= TFflag(bandi,Flag2, dataold)
	_,flag4, lab= TFflag(bandi,Flag2, datayoung)
	print('Total Objects in older= ', len(dataold))
	print('Bright Object Centers in older= ', len(flag1))
	print('Bright Objects in older= ', len(flag3))
	
	print('Total Objects in younger= ', len(datayoung))
	print('Bright Object Centers in younger= ', len(flag2))
	print('Bright Objects in younger= ', len(flag4))