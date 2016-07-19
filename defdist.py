print('Distribution plotting goes here')

def hist_mag(data1, data2, flags, band, col_name, outdir, crange):	#single aperture
	import matplotlib.pyplot as plt
	numflag=len(data1)
	strnumflag=str(numflag)
	numnot=len(data2)
	strnumnot=str(numnot)
	print(numflag, numnot)
	
	if not crange:
		crange='None'
	else:
		crange=str(crange)
	#tag=''
	if col_name is 'mag_cmodel':
		tag='cmodel'
		print(tag)
		mag1=data1[band+col_name]
		mag2=data2[band+col_name]
		xtitle='Aperture Magnitude'
		mytitle="Frequency of Different Aperture Magnitudes"
	elif col_name is'mag_aperture00':
		tag='ap00'
		print(tag)
		mag1=data1[band+col_name]
		mag2=data2[band+col_name]
		xtitle='Aperture Magnitude'
		mytitle="Frequency of Different Aperture Magnitudes"
	elif col_name is 'Z':
		tag='z'
		print('This is a histogram of redshifts')
		mag1=data1[col_name]
		mag2=data2[col_name]
		xtitle='Redshifts'
		mytitle="Frequency of Different Redshifts"
		#NOTE, mag1/2 represents Z shift, not magnitudes
	else:
		tag=''
		print('no tag')
	fig=plt.figure()

	nbins=8
	
	#print(min(mag1), max(mag1))
	if len(mag2)>0:
		if min(mag1)>min(mag2):
			mingr=min(mag2)
		else:
			mingr=min(mag1)
		if max(mag1)>max(mag2):
			maxgr=max(mag1)
		else:
			maxgr=max(mag2)
	else:
		mingr=min(mag1)
		maxgr=max(mag1)
		
	
	plt.hist(mag1,20,color='red', alpha=0.8)
	plt.hist(mag2,20,color='blue',alpha=0.7)
	plt.xlabel(xtitle, fontsize=10)
	plt.ylabel('Frequency', fontsize=10)
	plt.suptitle(mytitle, fontsize=15)
	plt.title('Range is '+crange, fontsize=10)
	plt.xlim(mingr, maxgr)
	plt.plot(0,0,label=flags[1] +'('+strnumflag+')', c='r', marker='^')
	plt.plot(0,0,label=flags[2] +'('+strnumnot+')', c='b', marker='*')
	plt.plot(0,0,label=col_name, c='k')
	plt.legend(loc=2,prop={'size':5})
	#plt.show()
	fig.savefig(outdir+band+tag+'_'+flags[3]+'_magdist.pdf')
	print(outdir+band+tag+'_'+flags[3]+'_magdist.pdf')
