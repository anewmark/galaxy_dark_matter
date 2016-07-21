print('Manipulating tables')
import astropy.table as table 
#from defcuts import *
import math
import numpy as np

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
	
data=table.Table.read(indir+'deg_vespa_HR_GAMA15.fits')
def rename_tab(data):
	names=data['SPECOBJID']
	starts=data['AGESTART']
	ustart=np.unique(starts)
	#print(un)
	for i in ustart:	#goes through each age start
		name=names[starts==i]
		uname, ind, inv,count=np.unique(name, return_index=True, return_counts=True,return_inverse=True )
		print(ind)
		for j in uname:	#goes through each individual name in each age start
			print(j)
			break
		break

rename_tab(data)