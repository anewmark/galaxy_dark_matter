print('Plot mass fraction as a function of Age')
import astropy.table as table 
#from defcuts import *
import math
import numpy as np
from def_get_mags import get_zdistmod
from def_ages import *
from def_age_plots import *
indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'
	
DATA=table.Table.read(indir+'small_vespa_LOWZ.fits')

print(np.ndim(DATA))
def get_agebin(Data,hm, plots=False):
	
	runid=Data['RUNID']
	runIDs, count=np.unique(runid, return_counts=True)

	run=5
	ndata=Data[runid==run] #only looking at first runID
	def age_bin(datas, tag=[]):
		#print(datas['SPECOBJID','AGESTART','MASS', 'M_STELLAR'])
		mass=datas['MASS']
		Tmass=datas['M_STELLAR']
		agestart=datas['AGESTART']
		ageend=datas['AGEEND']

		mass_fraction=mass/Tmass
	
		print('MF max= ', np.max(mass_fraction), 'MF min= ', np.min(mass_fraction))
	
	
		agebins=(ageend+agestart)/2
		ageranges=(ageend+agestart)
		agebin=no_repeats(agebins)
		agerange=no_repeats(ageranges)
	
		start=no_repeats(agestart)
		end=no_repeats(ageend)

		stack_mf, errors=stack_mass3(agestart, mass, Tmass, start)
	
			#errs=divide_error(mass, Tmass, datas['MASS_ERROR'], datas['M_STELLAR_ERROR'])
	
			#age_plot(abin, stack_mf, start, end, tag=tag) <-- DO NOT USE
		if plots==True:
			age_plot1(agebin, stack_mf, start, end, errors, tag=tag)
	
	newdata=mass_frac_cut1(ndata, hm, get_opp=False)

	per=[str(hm*100), '%']
	per=''.join(per)

	try:
		uname,ind, inv,count=np.unique(newdata['SPECOBJID'], return_index=True, return_counts=True,return_inverse=True)

		tagn=['+mf_cut80_', 'Number of Galaxies= '+str(len(uname)), 'Only Galaxies with Mass Fractions > '+per,'RunID= '+str(run)]
		#age_bin(newdata, tag=tagn)
		ndata=newdata
	except:
		uname,ind, inv,count=np.unique(ndata['SPECOBJID'], return_index=True, return_counts=True,return_inverse=True)
		tagn=['+all','Number of Galaxies= '+str(len(uname)),'All Galaxies','RunID= '+str(run)]

	age_bin(ndata, tag=tagn)
	
	print(tagn[1])

get_agebin(DATA, 0.585)