print('Finding Half Light Radii and Corresponding Slopes of Stacked Galaxies')

import astropy.table as table 
import numpy as np
from defcuts import *
from def_get_mags import *
from def_clean import *
from my_def_plots import *
from defflags import *
from defclump import * 

#Variables that can be changed

ty='mean'

tag='outcut'

txtdist= 'Figure2'
txtslope='Figure1'

outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'+ty+tag
doutdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/distribution/'+ty+tag
Flags=['flags_pixel_bright_object_center', 'brobj_cen_flag-', 'No Bright Ojbect Centers', 'Only Bright Object Centers', 'brobj_cen_flag']

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
	
	mlum1, mdens1, mrad1, merr1= get_avg_lums(lum1, rad1, lumd1, type=ty)
	
	hrad1= get_halflight(lum1, rad1)
	
	mhrad1= get_halflight(mlum1, mrad1)
	
	m1s, c1s, err1s= get_slopes(lum1, hrad1, rad1, lumd1, error=None, names=None, smax=True)
		
	m1, c1, radcut1, dencut1, sterr1, errcut1 =get_slopes(mlum1, mhrad1, mrad1, mdens1, error=merr1, names=None, smax=True)
	
	ynew1 = m1 * radcut1 + c1
	
	inds=[lum1, lumd1, rad1, hrad1]
	means=[mlum1,mdens1,mrad1,mhrad1, merr1]
	ind_slope=[m1s, c1s, err1s]
	mean_slopes=[m1, c1, radcut1, dencut1, ynew1,sterr1, errcut1]
	
	return inds, means, ind_slope, mean_slopes
	
def my_graphs(inds1, means1, ind_slope1, mean_slopes1, inds2, means2, ind_slope2, mean_slopes2):

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
		plt.scatter(x1, y1, color='r', marker='o',label='Not Flagged Galaxies')
		plt.plot(xcut1, yfit1, color='m', label='Fitted Not Flagged Galaxies: slope= '+str(m1)+' +- '+str(sterr1))
		plt.errorbar(x1, y1, yerr=error1, fmt='.',color='r')	

		plt.scatter(x2, y2, color='b', marker='o',label='Flagged Galaxies')
		plt.plot(xcut2, yfit2, color='c', label='Fitted Flagged Galaxies: slope= '+str(m2)+' +- '+str(sterr2))
		plt.errorbar(x2, y2, yerr=error2, fmt='.',color='b')

		plt.xlabel('Log Radii (kpc)')
		plt.ylabel('Luminosity Densities (Lsolar/kpc^2)')
		plt.title('Average Luminosity Densities v Radii')
		#plt.xlim(math.log10(1), math.log10(80))
		#plt.ylim(6,8.6)
		plt.legend(loc=0,prop={'size':6.0})
		f.text(0.05, 0.05, txtslope, color='red', weight='bold')
		outdirs=outdir+'TF.pdf'
		#plt.show()
		f.savefig(outdirs)
		print(outdirs)

	def dist_mean(m1s, m2s, m1, m2, sterr1, sterr2):
		import matplotlib.pyplot as plt
		import numpy as np
		import math
		figs=plt.figure()
		bs=np.linspace(-2.0,-1.4,num=15, endpoint=False)
		n1, b1, p1= plt.hist(m1s, bs, color='red', label='Not Flagged Galaxies ('+str(len(m1s))+')', alpha=0.8)
		n2, b2, p2= plt.hist(m2s,bs, color='blue', label='Flagged Galaxies ('+str(len(m2s))+')', alpha=0.8)
		
		print('Standard Deviation (Not Flagged): ', str(np.std(m1s)))
		print('Standard Deviation (Flagged): ', str(np.std(m2s)))
		
		plt.axvline(x=m1, color='magenta', label='Not Flagged Galaxies: slope= '+str(m1)+'  +- ' +str(sterr1))
		plt.axvline(x=m2, color='cyan', label='Flagged Galaxies: slope= '+str(m2)+' +- '+str(sterr2))
		plt.xlabel('Slopes', fontsize=10)
		plt.legend(loc=0,prop={'size':6.5})
		plt.ylabel('Frequency', fontsize=10)
		plt.title('With '+ty+' Slopes')
	
		outdirs=doutdir+'slopedist.pdf'
		figs.text(0.03, 0.03, txtdist, color='red', weight='bold')
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
	dist_mean(ind_slope1[0],ind_slope2[0],mean_slopes1[0],mean_slopes2[0],mean_slopes1[5], mean_slopes2[5])
	
	lum_mult_fit(means1[2], means2[2], means1[1], means2[1], mean_slopes1[2], mean_slopes2[2], mean_slopes1[4], mean_slopes2[4], mean_slopes1[5], mean_slopes2[5], mean_slopes1[0], mean_slopes2[0],means1[4], means2[4], outdir=outdir)	


indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'
datatab = table.Table.read(indir+ 'LOWZ_HSCGAMA15_apmgs+cmodmag.fits')

print(np.shape(datatab))
bands=['g', 'r', 'i','z', 'y']
parm=['flags_pixel_saturated_center','flags_pixel_edge','flags_pixel_interpolated_center','flags_pixel_cr_center','flags_pixel_suspect_center', 'flags_pixel_clipped_any','flags_pixel_bad']

daperture=[1.01,1.51,2.02,3.02,4.03,5.71,8.40,11.8,16.8,23.5]
aperture=[x*0.5 for x in daperture]
ne=[99.99, 199.99, 0.0]
mincut=0.1
maxcut=''
cutdata=not_cut(datatab, bands, 'mag_forced_cmodel', ne)
for b in range(0, len(bands)-1):
	newdata=many_flags(cutdata, parm, bands[b])	#flags not in y?
	cutdata=newdata
bandi=['i']

Flagdat, Notdat,lab= TFflag(bandi,Flags, newdata)
objid=Flagdat['object_id_1']

inds1, means1, ind_slope1, mean_slopes1=my_halflight2(Flagdat)
inds2, means2, ind_slope2, mean_slopes2=my_halflight2(Notdat)

my_graphs(inds1, means1, ind_slope1, mean_slopes1, inds2, means2, ind_slope2, mean_slopes2)