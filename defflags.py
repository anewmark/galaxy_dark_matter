print('Flag functions are here')

def many_flags(data, params, band):	#input more than one flag
    print('Initial length is ',len(data))
    nparam=len(params)
    for n in range(0, nparam):
        datas=data[data[band+params[n]]==False]
        print()
        data=datas
        print('Length after', params[n], 'is ', len(data))
    new_data=data
    return new_data
		
		
def TFflag(band, flgparams, data):  #add datazcut here
	dataname= flgparams[0]	#column name
	lab1=flgparams[1]	#No Flags
	lab2=flgparams[2]	#no flagged galaxies
	lab3=flgparams[3]	#Only Flagged Galaxies
	lab4=flgparams[4]	#outdir flag name
	#get rid of outlier data

	datatab_flag= data[data[band+dataname]==False] #gets rid of flagged

	nflag=len(datatab_flag)
	print('Number of OK galaxies= ', nflag)
	datatab_not=data[data[band+dataname]==True] #only flagged

	nnot=len(datatab_not)
	print('Number of flagged galaxies= ', nnot)
	labels=[lab1, lab2, lab3, lab4]
	return datatab_flag, datatab_not, labels