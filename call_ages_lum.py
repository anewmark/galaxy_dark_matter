import astropy.table as table 
#from defcuts import *
import math
import numpy as np
from def_ages import *

from defcuts import *
from def_get_mags import *
from def_clean import *
from my_def_plots import *
from defflags import *
from defclump import * 

ty='mean'

tag=''

txtdist= 'Figure2'
txtslope='Figure1'

outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'+ty+tag
doutdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/distribution/'+ty+tag

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

print('After all of these cuts, there are ', len(data1), 'older galaxies left')
print('After all of these cuts, there are ', len(data2), 'younger galaxies left')

def my_halflight2(dat1):
	lum1, rad1, lumd1= get_ind_lums(dat1, bands, aperture, scale='log')
	
	def upper_rad_cut(lum, rad, den): #this should get rid of galaxies outside 4r1/2
		from def_mymath import halflight
		nlum=[]
		nrad=[]
		nden=[]
		mult=4
		for x in range(len(rad)):
			lums=lum[x]
			rads=rad[x]
			dens=den[x]
			half=math.log10(10**np.max(lums)/2.0)
			hhx=halflight(rads,lums)
			
			hhx10=10**hhx
			hhx2s=mult*hhx10
			hhx2=math.log10(hhx2s)
			if np.max(rads) >= hhx2:
				mx=rads[(rads>=hhx)&(rads<=hhx2)]
				if len(mx)>=4:
					nlum.append(lums)
					nrad.append(rads)
					nden.append(dens)
				else:
					print('not enough data points')
			else:
				print('Upper limit out of range')
		nlum=np.array(nlum)
		nrad=np.array(nrad)
		nden=np.array(nden)
		return nlum, nrad, nden
	
		#print(len(lum1))	
	lum1, rad1, lumd1=upper_rad_cut(lum1, rad1, lumd1)
		#print(len(lum1))
		#print('Min rad= ', 10**np.min(rad1), 'max rad= ', 10**np.max(rad1))
	
	mlum1, mdens1, mrad1, merr1= get_avg_lums(lum1, rad1, lumd1, type=ty)
	
	hrad1= get_halflight(lum1, rad1)
	mhrad1= get_halflight(mlum1, mrad1)
	
	test_half=''
	if test_half=='go':
		bs=np.linspace(0,1.2,num=10, endpoint=False)
		plt.hist(hrad1,bins=bs, label='Total # Galaxies: '+str(len(hrad1)), color='red', alpha=0.9, zorder=2)
		plt.axvline(mhrad1, color='blue',label='Stacked half-light radius', zorder=3)
		plt.xlabel('Log10 Half-Light Radii')
		plt.legend(loc=0,prop={'size':6.5})
		plt.show()
	
	m1s, c1s, err1s= get_slopes(lum1, hrad1, rad1, lumd1, error=None, names=None, smax=True)
		
	m1, c1, radcut1, dencut1, sterr1, errcut1 =get_slopes(mlum1, mhrad1, mrad1, mdens1, error=merr1, names=None, smax=False)
	
	ynew1 = m1 * radcut1 + c1
	
	inds=[lum1, lumd1, rad1, hrad1]
	means=[mlum1,mdens1,mrad1,mhrad1, merr1]
	ind_slope=[m1s, c1s, err1s]
	mean_slopes=[m1, c1, radcut1, dencut1, ynew1,sterr1, errcut1]
	
	return inds, means, ind_slope, mean_slopes

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
		f.text(0.05, 0.05, txtslope, color='red', weight='bold')
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
			plt.plot(0,0, c='green', marker='*', label='K-S test is '+str(D))
			plt.xlim(np.min(M),-1.4)
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
		
		
			
	#slopevLmax(ind_slope1[0],ind_slope2[0], inds1[1], inds2[1])
	dist_mean(ind_slope1[0],ind_slope2[0],mean_slopes1[0],mean_slopes2[0],mean_slopes1[5], mean_slopes2[5], KS=False)
	
	lum_mult_fit(means1[2], means2[2], means1[1], means2[1], mean_slopes1[2], mean_slopes2[2], mean_slopes1[4], mean_slopes2[4], mean_slopes1[5], mean_slopes2[5], mean_slopes1[0], mean_slopes2[0],means1[4], means2[4], outdir=outdir)
	
inds1, means1, ind_slope1, mean_slopes1=my_halflight2(data1)
inds2, means2, ind_slope2, mean_slopes2=my_halflight2(data2)		
		
my_graphs(inds1, means1, ind_slope1, mean_slopes1, inds2, means2, ind_slope2, mean_slopes2)

test='flagged'
if test=='flagged':
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
	