print('Plot mass fraction as a function of Age')
import astropy.table as table 
#from defcuts import *
import math
import numpy as np
from def_get_mags import get_zdistmod
from def_ages import *
indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'

#if need to rewrite deg_vespa
def my_rad2deg(): 
	datatab = table.Table.read(indir+'vespa_HR_GAMA15.fits')
	def rad2deg(data, ra_col):
		import astropy.units as u
		rarad=datatab[ra_col]
		rarad.unit='radian'
		rarad=rarad.quantity
		col=rarad.to('deg')
		return col
	RA=rad2deg(datatab,'RA')
	DEC=rad2deg(datatab,'DEC')
	datatab['RA']=RA
	datatab['DEC']=DEC
	datatab.write('deg_vespa_HR_GAMA15.fits')
	
data=table.Table.read(indir+'nvespa_LOWZ_GAMA15-1.fits')

def mass_frac_cut(datas):
	ageend=datas['AGEEND']
	print('I THINK I DID IT WOOHOO')
	data1=datas[(ageend==np.max(ageend))&(datas['MASS']/datas['M_STELLAR']>=0.8)]

	umass,ind, inv,count=np.unique(data1['SPECOBJID','M_STELLAR'], return_index=True, return_counts=True,return_inverse=True)
	#a=np.in1d(datas['SPECOBJID','MASS','M_STELLAR'], data1['SPECOBJID','MASS','M_STELLAR'], assume_unique=False)
	a=np.in1d(datas['SPECOBJID','M_STELLAR'], data1['SPECOBJID','M_STELLAR'], assume_unique=False)
	newdat=datas[a==True]
	#print(newdat['SPECOBJID','AGESTART','MASS', 'M_STELLAR'])
	return newdat

def age_bin(datas, tag=[]):
	print(datas['SPECOBJID','AGESTART','MASS', 'M_STELLAR'])
	mass=datas['MASS']
	Tmass=datas['M_STELLAR']
	agestart=datas['AGESTART']
	ageend=datas['AGEEND']
	metal=datas['Z_1']
	
	mass_fraction=mass/Tmass
	print('max= ', np.max(mass_fraction))
	agebins=(ageend+agestart)/2
	ageranges=(ageend+agestart)
	agebin=no_repeats(agebins)
	agerange=no_repeats(ageranges)
	
	start=no_repeats(agestart)
	end=no_repeats(ageend)
	
	
	#stack_mf, abin=stack_mass(agebins, mass_fraction, agebin)
	
	def stack_mass2(x,y,bb):
		#bb is agestarts
		mass_sum=[]
		for n in bb:
			stmass=mass[x==n]
			stm=Tmass[x==n]
			sums=np.mean(stmass/stm)
			mass_sum.append(sums)
			print('mean=', sums)
		return mass_sum	
	
	stack_mf=stack_mass2(agestart, mass_fraction, start)
	
	#print(len(abin), len(agebin))
	
	
	#age_plot(abin, stack_mf, start, end, tag=tag)
	age_plot1(agebin, stack_mf, start, end, tag=tag)
	
	#xa=10**-4 #according to figure 16 in VESPA paper for red galaxies

#newdata=mass_frac_cut(data)

try:
	uname,ind, inv,count=np.unique(newdata['SPECOBJID'], return_index=True, return_counts=True,return_inverse=True)

	tagn=['+mf_cut80_', 'Number of Galaxies= '+str(len(uname)), 'Only Galaxies with Mass Fractions > 80%']
	#age_bin(newdata, tag=tagn)
	data=newdata
except:
	uname,ind, inv,count=np.unique(data['SPECOBJID'], return_index=True, return_counts=True,return_inverse=True)
	tagn=['+all','Number of Galaxies= '+str(len(uname)),'All Galaxies']

age_bin(data, tag=tagn)

