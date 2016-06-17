print('This program should plot average flux vs. aperture size for different bands')

indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'
outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'
if not os.path.exists(outdir):
	os.mkdir(outdir)
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
print(Naps)

bands= ['r', 'i', 'g','z', 'y']
band=bands[0]
##THE DIFFERENT FLAGS
no_sat=1	#
no_edge=0	#
no_interp=0	#
no_crcen=0	#
no_badpix=0	#
no_sus=0	#
no_clip=0	#
no_all=0	#
###THE DIFFERENT FLAGS

#getting the no flagged plots in z range[0.15,0.25]
Z=datatab['Z']
zmin=0.15
zmax=0.25

#dataclump=datatab
datatab=datatab[(Z>=zmin) & (Z<=zmax)]
nclump=len(datatab)
ndat=str(nclump)
rangez=[zmin, zmax]
rangez=str(rangez)
clump=[]
clumpflag=[]
clumpnot=[]

if no_sat:
	datatab_flag= datatab[(datatab[band+'flags_pixel_saturated_center']==False) & (datatab[band+'mag_aperture00']<50)& (datatab[band+'mag_aperture01']<50)& (datatab[band+'mag_aperture02']<50) & (datatab[band+'mag_aperture03']<50) & (datatab[band+'mag_aperture04']<50) & (datatab[band+'mag_aperture05']<50) & (datatab[band+'mag_aperture06']<50) & (datatab[band+'mag_aperture07']<50) & (datatab[band+'mag_aperture08']<50) & (datatab[band+'mag_aperture09']<50)]
	#print(datatab_flag.colnames)
	Z=datatab_flag['Z']
	clump_flag=datatab_flag[(Z>=zmin) & (Z<=zmax)]
	nclumpflag=len(clump_flag)
	nflag=str(nclumpflag)
	datatab_not = datatab[(datatab[band+'flags_pixel_saturated_center']==True) & (datatab[band+'mag_aperture00']<50)& (datatab[band+'mag_aperture01']<50)& (datatab[band+'mag_aperture02']<50) & (datatab[band+'mag_aperture03']<50) & (datatab[band+'mag_aperture04']<50) & (datatab[band+'mag_aperture05']<50) & (datatab[band+'mag_aperture06']<50) & (datatab[band+'mag_aperture07']<50) & (datatab[band+'mag_aperture08']<50) & (datatab[band+'mag_aperture09']<50)]
	Z=datatab_not['Z']
	clump_not=datatab_not[(Z>=zmin) & (Z<=zmax)]
	nclumpnot=len(clump_not)
	nnot=str(nclumpnot)
	labs=['No Flags', 'No Saturated Galaxies', 'Only Saturated Galaxies', 'sat_flag']
	print('numbers are', nclump, nclumpflag, nclumpnot)

if no_edge:
	outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	datatab_flag = datatab[datatab[band+'flags_pixel_edge']==False ]
	Z=datatab_flag['Z']
	clump_flag=datatab_flag[(Z>=zmin) & (Z<=zmax)]
	nclumpflag=len(clump_flag)
	nflag=str(nclumpflag)
	datatab_not = datatab[datatab[band+'flags_pixel_edge']==True]
	Z=datatab_not['Z']
	clump_not=datatab_not[(Z>=zmin) & (Z<=zmax)]
	nclumpnot=len(clump_not)
	nnot=str(nclumpnot)
	labs=['No Flags', 'No Edge Galaxies', 'Only Edge Galaxies', 'edge_flag']
	print('numbers are', nclump, nclumpflag, nclumpnot)
	
if no_interp:
	outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	datatab_flag = datatab[datatab[band+'flags_pixel_interpolated_center']==False ]
	Z=datatab_flag['Z']
	clump_flag=datatab_flag[(Z>=zmin) & (Z<=zmax)]
	nclumpflag=len(clump_flag)
	nflag=str(nclumpflag)
	datatab_not = datatab[datatab[band+'flags_pixel_interpolated_center']==True]
	Z=datatab_not['Z']
	clump_not=datatab_not[(Z>=zmin) & (Z<=zmax)]
	nclumpnot=len(clump_not)
	nnot=str(nclumpnot)
	labs=['No Flags', 'No Interpolated Centers', 'Only Interpolated Centers', 'interp_flag']
	print('numbers are', nclump, nclumpflag, nclumpnot)

if no_crcen:
	outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	datatab_flag = datatab[datatab[band+'flags_pixel_cr_center']==False ]
	Z=datatab_flag['Z']
	clump_flag=datatab_flag[(Z>=zmin) & (Z<=zmax)]
	nclumpflag=len(clump_flag)
	nflag=str(nclumpflag)
	datatab_not = datatab[datatab[band+'flags_pixel_cr_center']==True]
	Z=datatab_not['Z']
	clump_not=datatab_not[(Z>=zmin) & (Z<=zmax)]
	nclumpnot=len(clump_not)
	nnot=str(nclumpnot)
	labs=['No Flags', 'No Interpolated Centers', 'Only Interpolated Centers', 'crcen_flag']
	print('numbers are', nclump, nclumpflag, nclumpnot)
	
if no_badpix:
	outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	datatab_flag = datatab[datatab[band+'flags_pixel_bad']==False ]
	Z=datatab_flag['Z']
	clump_flag=datatab_flag[(Z>=zmin) & (Z<=zmax)]
	nclumpflag=len(clump_flag)
	nflag=str(nclumpflag)
	datatab_not = datatab[datatab[band+'flags_pixel_bad']==True]
	Z=datatab_not['Z']
	clump_not=datatab_not[(Z>=zmin) & (Z<=zmax)]
	nclumpnot=len(clump_not)
	nnot=str(nclumpnot)
	labs=['No Flags', 'No Bad Pixels', 'Only Bad Pixels', 'badpix_flag']
	print('numbers are', nclump, nclumpflag, nclumpnot)

if no_sus:
	outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	datatab_flag = datatab[datatab[band+'flags_pixel_suspect_center']==False ]
	Z=datatab_flag['Z']
	clump_flag=datatab_flag[(Z>=zmin) & (Z<=zmax)]
	nclumpflag=len(clump_flag)
	nflag=str(nclumpflag)
	datatab_not = datatab[datatab[band+'flags_pixel_suspect_center']==True]
	Z=datatab_not['Z']
	clump_not=datatab_not[(Z>=zmin) & (Z<=zmax)]
	nclumpnot=len(clump_not)
	nnot=str(nclumpnot)
	labs=['No Flags', 'No Suspect Center', 'Only Suspect Center', 'suspect_flag']
	print('numbers are', nclump, nclumpflag, nclumpnot)

if no_clip:
	outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	datatab_flag = datatab[datatab[band+'flags_pixel_clipped_any']==False ]
	Z=datatab_flag['Z']
	clump_flag=datatab_flag[(Z>=zmin) & (Z<=zmax)]
	nclumpflag=len(clump_flag)
	nflag=str(nclumpflag)
	datatab_not = datatab[datatab[band+'flags_pixel_clipped_any']==True]
	Z=datatab_not['Z']
	clump_not=datatab_not[(Z>=zmin) & (Z<=zmax)]
	nclumpnot=len(clump_not)
	nnot=str(nclumpnot)
	labs=['No Flags', 'No Clipped (Any)', 'Only Clipped (Any)', 'clip_flag']
	print('numbers are', nclump, nclumpflag, nclumpnot)
	
if no_all:
	outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	datatab_flag = datatab[(datatab[band+'flags_pixel_clipped_any']==False)&(datatab[band+'flags_pixel_suspect_center']==False)& (datatab[band+'flags_pixel_bad']==False)& (datatab[band+'flags_pixel_interpolated_center']==False)&(datatab[band+'flags_pixel_interpolated_center']==False)&(datatab[band+'flags_pixel_edge']==False)&(datatab[band+'flags_pixel_saturated_center']==False)]
	Z=datatab_flag['Z']
	clump_flag=datatab_flag[(Z>=zmin) & (Z<=zmax)]
	nclumpflag=len(clump_flag)
	nflag=str(nclumpflag)
	datatab_not = datatab[datatab[band+'flags_pixel_clipped_any']==True]
	Z=datatab_not['Z']
	clump_not=datatab_not[(Z>=zmin) & (Z<=zmax)]
	nclumpnot=len(clump_not)
	nnot=str(nclumpnot)
	labs=['No Flags', 'No Galaxies Flagged', 'Only Galaxies Flagged)', 'all_flag']
	print('numbers are', nclump, nclumpflag, nclumpnot)

for a in range(0,Naps):
	js=str(a)
	SB10tot=[]
	for n in range(0,nclump):
		magclump=datatab[band+'mag_aperture0'+js]
		#print(magclump[n])
		sbclump=magclump+2.5*math.log10(4*pi*aperture[a]**2)
		SB10=math.pow(10,sbclump[n]) #getting 10^SB
		SB10tot.append(SB10)
	SB10sum=np.sum(SB10tot)
	#print('SB10sum is ', SB10sum)
	SB10avg=SB10sum/nclump
	#print('SB10avg is ', SB10avg)
	sbavg=math.log10(SB10avg)
	#print('sbavg is ', sbavg)
	#sbavg=sbsum/nclump
	clump.append(sbavg)
	#break
	Flxflag=[]
	for n in range(0,nclumpflag):
		magclump=datatab[band+'mag_aperture0'+js]
		sbclump=magclump+2.5*math.log10(4*pi*aperture[a]**2)
		F=math.pow(10,sbclump[n])
		Flxflag.append(F)
	Flxflagsum=np.sum(Flxflag)
	sbflagavg=math.log10(Flxflagsum/nclumpflag)
	#sbflagavg=sbflagsum/nclumpflag
	clumpflag.append(sbflagavg)
	
	if nclumpnot==0:
		clumpnot=np.zeros(10)
	else:
		Flxnot=[]
		for n in range(0,nclumpnot):
			magclump=datatab[band+'mag_aperture0'+js]
			sbclump=magclump+2.5*math.log10(4*pi*aperture[a]**2)
			F=math.pow(10,sbclump[n])
			Flxnot.append(F)
		Flxnotsum=np.sum(Flxnot)
		sbnotsum=math.log10(Flxnotsum)
		sbnotavg=sbnotsum/nclumpnot
		clumpnot.append(sbnotavg)
#print(clump, clumpflag, clumpnot)
	


f=plt.figure()
plt.subplots_adjust(bottom=0.1, top=.9, right=.93, left=.07)

ax0=plt.subplot(131)
plt.scatter(aperture, clump, c='r', marker='^')
plt.xlabel('Aperture Radius (arcseconds)', fontsize=4)
plt.ylabel('Averaged Surface Brightness (mag/arcsec^2)', fontsize=4)
plt.title("Averaged Surface Brightness vs. Aperture Radius in "+band, fontsize=7)
plt.ylim(min(clump),max(clump))
plt.tick_params(axis='both', which='major', labelsize=4)
ax0.invert_yaxis()
plt.plot(0,0, label='Redshift is '+rangez, marker='', c='k')
plt.plot(0,0, label=labs[0], marker='', c='k')
plt.plot(0,0,label='Number of Galaxies= '+ndat, marker='', c='k')
plt.legend(loc=7,prop={'size':4}, numpoints = 1)

ax1=plt.subplot(132)
plt.scatter(aperture,clumpflag, c='b', marker='+')
plt.xlabel('Aperture Radius (arcseconds)', fontsize=4)
plt.ylabel('Averaged Surface Brightness (mag/arcsec^2)', fontsize=4)
plt.title("Averaged Surface Brightness vs. Aperture Radius in "+band, fontsize=7)
plt.ylim(min(clumpflag),max(clumpflag))
plt.tick_params(axis='both', which='major', labelsize=4)
ax1.invert_yaxis()
plt.plot(0,0, label='Redshift is '+rangez, marker='', c='k')
plt.plot(0,0, label=labs[1], marker='', c='k')
plt.plot(0,0,label='Number of Galaxies= ' +nflag, marker='', c='k')
plt.legend(loc=7,prop={'size':4})

ax2=plt.subplot(133)
plt.scatter(aperture, clumpnot, c='g', marker='*')
plt.xlabel('Aperture Radius (arcseconds)', fontsize=4)
plt.ylabel('Averaged Surface Brightness (mag/arcsec^2)', fontsize=4)
plt.title("Averaged Surface Brightness vs. Aperture Radius in "+band, fontsize=7)
plt.ylim(min(clumpnot),max(clumpnot))
plt.tick_params(axis='both', which='major', labelsize=4)
ax2.invert_yaxis()
plt.plot(0,0, label='Redshift is '+rangez, marker='', c='k')
plt.plot(0,0, label=labs[2], marker='', c='k')
plt.plot(0,0,label='Number of Galaxies= ' +nnot, marker='', c='k')
plt.legend(loc=7,prop={'size':4})
plt.show()
f.savefig(outdir+band+'_'+labs[3]+'_zclump.pdf')

#for same graph
fig=plt.figure()
plt.plot(aperture, clump, c='r', marker='^')
plt.plot(aperture,clumpflag, c='b',  marker='+')
plt.plot(aperture, clumpnot, c='g',marker='*')
plt.xlabel('Aperture Radius (arcseconds)', fontsize=10)
plt.ylabel('Averaged Surface Brightness (mag/arcsec^2)', fontsize=10)
plt.suptitle("Averaged Surface Brightness vs. Aperture Radius in "+band, fontsize=15)
plt.title("Redshift range "+ rangez, fontsize=12)
red_patch = mpatches.Patch(color='red', label=labs[0] +'('+ndat+')')
blue_patch = mpatches.Patch(color='blue', label=labs[1] +'('+nflag+')')
green_patch = mpatches.Patch(color='green', label=labs[2] +'('+nnot+')')
#plt.plot(0,0, label='Zrange is '+rangez, marker='*')
plt.legend(handles=[red_patch, blue_patch, green_patch], loc=7,prop={'size':5})
#plt.legend(loc=7,prop={'size':5})
plt.gca().invert_yaxis()
plt.show()
fig.savefig(outdir+band+'_'+labs[3]+'_zcomboclump.pdf')