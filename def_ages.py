import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import math

def mass_frac_cut(datas, mn):
	ageend=datas['AGEEND']
	print('I THINK I DID IT WOOHOO')
	
	TM=get_tot_mass(datas, 'MASS','SPECOBJID')
	#data1=datas[(ageend==np.max(ageend))&(datas['MASS']/datas['M_STELLAR']>=mn)]
	
	data1=datas[(ageend==np.max(ageend))&(datas['MASS']/TM>=mn)]
	
	a=np.in1d(datas['SPECOBJID','M_STELLAR'], data1['SPECOBJID','M_STELLAR'], assume_unique=False)
	newdat=datas[a==True]
	#print(newdat['SPECOBJID','AGESTART','MASS', 'M_STELLAR'])
	print('Length of data after mf mincut: ', len(newdat))
	return newdat
	
def mass_frac_cut1(datas, mn, get_opp=False):  #original way
	ageend=datas['AGEEND']
	print('I THINK I DID IT WOOHOO')
	#print(datas['SPECOBJID','AGESTART','MASS', 'M_STELLAR']) #check to make sure it has more than data1
	data1=datas[(ageend==np.max(ageend))&(datas['MASS']/datas['M_STELLAR']>=mn)]

	a=np.in1d(datas['SPECOBJID','M_STELLAR'], data1['SPECOBJID','M_STELLAR'], assume_unique=False)
	newdat=datas[a==True]
	
	#databad=newdat[newdat['MASS']>newdat['M_STELLAR']]
	if get_opp==True:
		notdat=datas[a==False]
		print('Length of data after mf mincut: ', len(newdat))
		print('Length of other data: ', len(notdat))
		return newdat, notdat
	else:
		print('Length of data after mf mincut: ', len(newdat))
		return newdat
	
def divide_error(dat1, dat2, er1, er2):
	#want to find error of m+-dm/Tm+-dTm
	div_err= lambda m, TM, merr, TMerr: m/TM*math.sqrt((merr/m)**2+(TMerr/TM)**2)
	
	errs=[]
	for i in range(len(dat1)):
		err1=div_err(dat1[i], dat2[i], er1[i], er2[i])
		errs.append(err1)
	#print('errors= ', errs)
	where_nans= np.isnan(errs)
	print('Find nans', where_nans)
	errs[where_nans]=0
	print('errors= ', errs)
	return(errs)
def no_repeats(x):
	print('Length with repeats: ', len(x))
	xs=[]
	for i in x:
		if i not in xs:
			xs.append(i)
	print('Length without repeats: ', len(xs))
	return xs
def common_element(list1, list2):
	result=[]
	for el in list1:
		if el in list2:
			result.append(el)
	result=np.array(result)
	return result
def test_sum(data, col1, col2, namecol):
	name=data[namecol]
	onename=no_repeats(name)
	badgal=[]
	for n in onename:
		newd=data[name==n]
		#print(len(newd))		#why is this 64 and not 16--> 4 trials
		nmass=newd[col1]
		nTmass=newd[col2]
		sing,ind, inv, counts=np.unique(nTmass, return_index=True, return_counts=True,return_inverse=True)
		if len(sing)==2:
			totnums=sing*2
		elif len(sing)==1:
			totnums=sing*4
		elif len(sing)==3:
			totnums=sing
		elif len(sing)==4:
			totnums=sing
		else:
			print('Number of different total masses is actually: ', len(sing))
		sumtotnums=np.sum(totnums)
		sum1=np.sum(nmass)
		print('Galaxy is ', n, 'with mass= ', sum1, ' of total mass= ', sumtotnums)
		if sum1>sumtotnums:
			badgal.append(n)
		#break

	print('These are the bad galaxies: ', len(badgal), len(onename))
def data_cut(data, col1, col2):
	m=data[col1]
	Tm=data[col2]
	databad=data[m>Tm]
	datagood=data[m<=Tm]
	#datagood=datas[mass<=Tmass]
	print('Number of bad lines: ', len(databad))
	print('Number of good lines: ', len(datagood))
	print('Number of total lines: ', len(data))
	
	#print(databad['SPECOBJID','AGESTART', 'AGEEND', col1, col2])
	return datagood
		
def get_tot_mass(data, col1, namecol):
	#gives the sum of the masses in each age bin
	Totmass=[]
	names=data[namecol]
	name=np.unique(names)
	for n in name:
		#newd=data[names==names[n]]
		newd=data[names==n]
		nmass=newd[col1]
		numdat=len(newd)

		div=numdat/16 #total number divided by 16 bins
		#print('divide= ', div)
		M=np.sum(nmass)
		MM=M/div
		#Totmass.append(MM)
		for n in range(numdat):
			Totmass.append(MM)
		#break
	TM=np.array(Totmass)
	#print(TM)
	return TM
		
def get_tot_mass1(data, col1, namecol):
	#gives the sum of the masses in each age bin
	Totmass=[]
	names=data[namecol]
	name, c=np.unique(names,return_counts=True)
	#print(c)
	#print(len(name))
	way=1
	if way==1:
		for n in range(len(names)):
			newd=data[names==names[n]]
			#print(names[n])
			nmass=newd[col1]
			numdat=len(newd)
			#print(newd)
			#print(numdat)
			div=numdat/16 #total number divided by 16 bins
			#print('divide= ', div)
			#print('masses= ', nmass)
			M=np.sum(nmass)
			MM=M/div
			Totmass.append(MM)
			#break
		TM=np.array(Totmass)
	else:
		for n in name:
			newd=data[names==n]
			nmass=newd[col1]
			numdat=len(newd)
			div=numdat/16 
			M=np.sum(nmass)
			MM=M/div
			for i in range(numdat):
				Totmass.append(MM)
	TM=np.array(Totmass)
	#print(TM)
	return TM
	
def stack_mass3(ages,mass, Tmass,bb):
	#bb, x is agestarts
	mass_sum=[]
	std_err=[]
	for bin in bb:
		print('Age Start= ', bin)
		stmass=mass[ages==bin]
		stm=Tmass[ages==bin]
		mf=stmass/stm
		#print(len(mf))
		#print(mf)
		means=np.mean(mf)
		#mean1=np.sum(mf)/len(stm) #used to check to make sure means=mean1
		err=np.std(mf)/len(mf)
		std_err.append(err) #errors on the mean in each individual bin
		mass_sum.append(means)
		print('mean=', means)
		#break
	return mass_sum, std_err
