print('Cut Functions are here')

def z_cut(datatable, zmin, zmax):
    Z=datatable['Z']
    print('Length before z cut is ', len(Z))
    Zcut=datatable[(Z>=zmin)&(Z<=zmax)]
    zlen=len(Zcut)
    print('Length after zcut is ', zlen)
    rangez=[zmin, zmax]
    return Zcut, rangez #this returns new datatable within Z range
    
def out_cut(data, band, colname, mincut, maxcut):	#cuts the outliers
	print('Length before cut is ', len(data))
	if mincut:
		datamin=data[data[band+colname]>mincut]
		print(len(datamin))
	else:
		datamin=data
		print('No Min Given')
	if maxcut:
		datamax=datamin[datamin[band+colname]>maxcut]
		print(len(datamax))
	else:
		datamax=datamin
		print('No Max Given')
	datacut=datamax
	print('Length after cut is ', len(datacut))
	return datacut