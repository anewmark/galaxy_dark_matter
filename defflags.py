print('Flag functions are here')

def many_flags(data, params, band):	#input more than one flag
    print('Initial length is ',len(data))
    nparam=len(params)
    for n in range(0, nparam):
        datas=data[data[band+params[n]]==False]
        print()
        data=datas
        print('Length after', band, params[n], 'is ', len(data))
    new_data=data
    return new_data
		
		
def TFflag(band, flgparams, data):  
	dataname= flgparams[0]	#column name
	lab1=flgparams[1]	#No Flags
	lab2=flgparams[2]	#no flagged galaxies
	lab3=flgparams[3]	#Only Flagged Galaxies
	lab4=flgparams[4]	#outdir flag name
	#get rid of outlier data
	dataflag=data
	datanot=data
	for n in range(0, len(band)):
		datatab_flag= data[dataflag[band[n]+dataname]==False] #gets rid of flagged
		datatab_not=data[datanot[band[n]+dataname]==True] #only flagged
		
		print('Number of OK galaxies in band ',band[n],'= ', len(datatab_flag))
		print('Number of flagged galaxies in band ' ,band[n],'= ', len(datatab_not))
		dataflag=datatab_flag
		datanot=datatab_not
	
	datatab_flag=dataflag
	datatab_not=datanot
	
	labels=[lab1, lab2, lab3, lab4]
	
	return datatab_flag, datatab_not, labels