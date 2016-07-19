print('Cut Functions are here')

def z_cut(data, zmin, zmax):
	#Z=data['Z']
	print('Length before z cut is ', len(data['Z']))
	if zmin:
		datamin=data[data['Z']>=zmin]
		print('length after zmin cut', len(datamin))
	else:
		datamin=data
		print('No ZMin Given')
		zmin='-Infinity'
	if zmax:
		datamax=datamin[datamin['Z']<=zmax]
		print('length after zmax cut= ', len(datamax))
	else:
		datamax=datamin
		print('No ZMax Given')
		zmax='Infinity'
	Zcut=datamax
	#Zcut=datatable[(Z>=zmin)&(Z<=zmax)]
	zlen=len(Zcut)
	print('Length after zcut is ', zlen)
	rangez=[zmin, zmax]
	return Zcut, rangez #this returns new datatable within Z range
    
def out_cut(data, band, colname, mincut, maxcut):	#cuts the outliers
	print('Length before cut is ', len(data))
	
	if mincut:
		datamin=data[data[band+colname]>mincut]
		print('length after min cut', len(datamin))
	else:
		datamin=data
		mincut='-Infinity'
		print('No Min Given')
	if maxcut:
		datamax=datamin[datamin[band+colname]>maxcut]
		print('length after max cut= ', len(datamax))
	else:
		datamax=datamin
		maxcut='Infinity'
		print('No Max Given')
	datacut=datamax
	print('Length after cut is ', len(datacut))
	rangec=str([mincut, maxcut])
	return datacut, rangec

#if different bands
def not_cut(data, band, colname, ne):
	N=len(ne)
	nbands=len(band)
	print('Length before notequal cuts= ', len(data))
	for i in range(0,nbands):
		for n in range(0,N):
			datas=data[data[band[i]+colname]!=ne[n]]
			data=datas
	new_data=data
	print('Length after notequal cuts= ', len(new_data))
	return new_data
	
#if no bands
def not_cut1(data, colname, ne):
	N=len(ne)
	print('Length before notequal cuts= ', len(data))
	
	if N>=2:
		for n in range(0,N):
			datas=data[data[colname]!=ne[n]]
			data=datas
	else:
		print('array of 1')
		data=data[data[colname]!=ne]
	new_data=data
	print('Length after notequal cuts= ', len(new_data))
	return new_data