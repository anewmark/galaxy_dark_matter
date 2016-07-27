import astropy.table as table 
#from defcuts import *
import math
import numpy as np
import matplotlib.pyplot as plt
from def_ages import mass_frac_cut1, no_repeats, stack_mass3
from my_style import get_presentation
indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'

get_presentation()
DATA=table.Table.read(indir+'small_vespa_LOWZ.fits')

def get_agebin_datold(Data, hm):
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
	#gives data where it masses formed in diff bins
	return newdata, notdat

hm=0.6

per=[str(hm*100), '%']
per=''.join(per)

odata, ydata= get_agebin_datold(DATA, hm)

def age_bin(datas):
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

		return stack_mf, errors, agebin, start, end
		
mf_old, err_old, agebin, start, end=age_bin(odata)
mf_young, err_young,_,_,_=age_bin(ydata)
num1=np.unique(odata['SPECOBJID'])
num2=np.unique(ydata['SPECOBJID'])

def plot_stack_ages(agebin, mf1, mf2, err1, err2,start, end, num1, num2):
	doutdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/distribution/'
	start=np.array(start)
	end=np.array(end)
	arange=end-start
	arange=arange

	bwidth=arange
	print('bar width=', bwidth)
	print(len(bwidth))
	f=plt.figure()
	#plt.style.use('presentation')
	labels=np.append(start, end[len(end)-1])
	#print('labels', labels)
	#xspace=np.logspace(np.log10(np.min(labels)), np.log10(np.max(labels)), num=len(labels))
	label=[str(n) for n in labels]
	xsp=[math.log10(n) for n in labels]
	#print(xsp)
	
	#print(label)
	#plt.bar(x, y, width=bwidth, align='center', color='None')
	plt.bar(agebin, mf1, width=bwidth, align='center', color='red', alpha=0.85, label='Older Galaxies (Mass Fractions >'+per+')')
	plt.bar(agebin, mf2, width=bwidth, align='center', color='blue', alpha=0.7,label='Younger Galaxies (Mass Fractions <'+per+')', zorder=2)
	plt.errorbar(agebin, mf1, yerr=err1,label='Standard Error on the Mean (Older)', fmt='.', color='m')
	plt.errorbar(agebin, mf2, yerr=err2,label='Standard Error on the Mean (Younger)', fmt='.', color='c')
	plt.xlabel('Lookback Time (Gyr)')
	plt.xscale('log')
	plt.yscale('log')
	plt.xticks(xsp, label=label, rotation='vertical')

	plt.ylabel('Stacked Mass Fractions')
	plt.title('Age vs. Mass Fractions')
	#plt.xlim(np.min(x)-bwidth[0], np.max(x)+bwidth[len(bwidth)-1]/2.0)
	plt.xlim(np.min(start), np.max(end))


	plt.plot(0,0, label='# Older Galaxies= '+str(len(num1)), c='r')
	plt.plot(0,0, label='# Younger Galaxies=' +str(len(num2)),c='b')
	plt.legend(loc=2,prop={'size':8.0})
	
	plt.show()
	outdirs=doutdir+'oy_agebin.pdf'
	f.savefig(outdirs)
	print(outdirs)

	
plot_stack_ages(agebin, mf_old, mf_young, err_old, err_young,start, end, num1, num2)
