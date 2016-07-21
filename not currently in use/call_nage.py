print('Plot mass fraction as a function of Age')
import astropy.table as table 
#from defcuts import *
import math
import numpy as np
from def_get_mags import get_zdistmod
from def_ages import *
from def_age_plots import *
indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'
	
Data=table.Table.read(indir+'small_vespa_LOWZ.fits')

runid=Data['RUNID']
runIDs, count=np.unique(runid, return_counts=True)

run=1
data=Data[runid==run] #only looking at first runID

	#hidden
	#test_sum(data, 'MASS', 'M_STELLAR', 'SPECOBJID')
	#ndata=data_cut(data, 'MASS', 'M_STELLAR')

def age_bin(datas, tag=[]):
	#print(datas['SPECOBJID','AGESTART','MASS', 'M_STELLAR'])
	mass=datas['MASS']
	Tmass=datas['M_STELLAR']
	agestart=datas['AGESTART']
	ageend=datas['AGEEND']
		
	TM=get_tot_mass1(datas, 'MASS','SPECOBJID')

	mass_fraction=mass/Tmass
	
	print('MF max= ', np.max(mass_fraction), 'MF min= ', np.min(mass_fraction))
	#mass_fraction=mass/Tmass
	
	agebins=(ageend+agestart)/2
	ageranges=(ageend+agestart)
	agebin=no_repeats(agebins)
	agerange=no_repeats(ageranges)
	
	start=no_repeats(agestart)
	end=no_repeats(ageend)
	
	def stack_mass2(ages,mfs,bb):
		#bb, x is agestarts
		mass_sum=[]
		std_err=[]
		for bin in bb:
			print('Age Start= ', bin)
			#stmass=mass[x==bin]
			#stm=Tmass[x==bin]
			#print('Number of galaxies in this bin= ', len(stmass))
			mf=mfs[ages==bin]
			#print(len(mf))
			#mf=stmass/stm
			means=np.mean(mf)
			#mean1=np.sum(mf)/len(stm) #used to check to make sure means=mean1
			err=np.std(mf)/len(mf)
			std_err.append(err) #errors on the mean in each individual bin
			mass_sum.append(means)
			print('mean=', means)
			#break
		return mass_sum, std_err
	
	stack_mf, errors=stack_mass2(agestart, mass_fraction, start)
	
		#errs=divide_error(mass, Tmass, datas['MASS_ERROR'], datas['M_STELLAR_ERROR'])
	
		#age_plot(abin, stack_mf, start, end, tag=tag) <-- DO NOT USE
	age_plot1(agebin, stack_mf, start, end, errors, tag=tag)

ndata=data
	
#newdata=mass_frac_cut(ndata, 0.8)
try:
	uname,ind, inv,count=np.unique(newdata['SPECOBJID'], return_index=True, return_counts=True,return_inverse=True)

	tagn=['+mf_cut80_', 'Number of Galaxies= '+str(len(uname)), 'Only Galaxies with Mass Fractions > 80%','RunID= '+str(run)]
	#age_bin(newdata, tag=tagn)
	ndata=newdata
except:
	uname,ind, inv,count=np.unique(ndata['SPECOBJID'], return_index=True, return_counts=True,return_inverse=True)
	tagn=['+all','Number of Galaxies= '+str(len(uname)),'All Galaxies','RunID= '+str(run)]

age_bin(ndata, tag=tagn)



