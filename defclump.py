print("Clump plotting is here")


def clump_plot(aper,dataflag, datanot, band, magname, labs, outdir,zrange):
	import matplotlib.pyplot as plt
	import numpy as np
	import math
	import matplotlib.patches as mpatches
	pi=math.pi
	Naps=len(aper)
	nflag=len(dataflag)
	numflag=str(nflag) #string version
	nnot=len(datanot)
	numnot=str(nnot) #string version
	clumpflag=[]
	clumpnot=[]
	for a in range(0, Naps):
		js=str(a)
		SB10flag=[]
		for n in range(0,nflag):
			magclump=dataflag[band+magname+js]
			sbclump=magclump+2.5*math.log10(4*pi*aper[a]**2)
			SB10=math.pow(10,sbclump[n]) #getting 10^SB
			SB10flag.append(SB10)
		SB10sum=np.sum(SB10flag)
		#print('SB10sum is ', SB10sum)
		SB10avg=SB10sum/nflag
		#print('SB10avg is ', SB10avg)
		sbavg=math.log10(SB10avg)
		clumpflag.append(sbavg)
		#print('sbclumpflag is ',sbclump)
		if nnot==0:
			clumpnot=np.zeros(Naps)
		else:
			SB10not=[]
			for n in range(0,nnot):
				magclump=datanot[band+magname+js]
				sbclump=magclump+2.5*math.log10(4*pi*aper[a]**2)
				SB10=math.pow(10,sbclump[n]) #getting 10^SB
				SB10not.append(SB10)
			SB10notsum=np.sum(SB10not)
			#print('SB10sum is ', SB10sum)
			SB10notavg=SB10notsum/nnot
			#print('SB10avg is ', SB10avg)
			sbnotavg=math.log10(SB10notavg)
			clumpnot.append(sbnotavg)
		#print('sbclumpnot is ',sbclump)
	#print(nflag, 'clumpflag is ', clumpflag)
	#print(nnot, 'clumpnot is ', clumpnot)
		
	#print(clumpflag, clumpnot)
	if not zrange:
		zrange='None Specified'
	else:
		zrange=str(zrange)
	f, (ax0, ax1) = plt.subplots(1,2, sharey=True)
	f.suptitle("Mean Surface Brightness vs. Aperture Radius in "+band)
	ax0=plt.subplot(121)
	ax0.plot(aper, clumpflag, c='r', marker='^')
	ax0.set_xlabel('Aperture Radius (arcseconds)', fontsize=6.5)
	ax0.set_ylabel('Mean Surface Brightness (mag/arcsec^2)', fontsize=6.5)
	ax0.set_title('Number of Galaxies Not Flagged= '+numflag, fontsize=9)
	plt.ylim(min(clumpflag),max(clumpflag))
	plt.tick_params(axis='both', which='major', labelsize=5)
	ax0.invert_yaxis()
	ax0.plot(0,0, label='Redshift is '+zrange, marker='', c='k')
	#ax0.plot(0,0, label=labs[1], marker='', c='k')
	#ax0.plot(0,0,label='Number of Galaxies= '+numflag, marker='', c='k')
	ax0.legend(loc=7,prop={'size':5}, numpoints = 1)
		
	ax1=plt.subplot(122)
	if nnot==0:
		ax1.set_title('Number of Flagged Galaxies= '+numnot, fontsize=9)
		#ax1.plot(0,0,label='No Galaxies were Flagged', marker='', c='k')
		#ax1.legend(loc=9,prop={'size':8}, numpoints = 1)
		#plt.tick_params(axis='both', which='major', labelsize=5)
		ax1.invert_yaxis()
	else:
		ax1.plot(aper, clumpnot, c='b', marker='*')
		ax1.set_xlabel('Aperture Radius (arcseconds)', fontsize=6.5)
		ax1.set_ylabel('Mean Surface Brightness (mag/arcsec^2)', fontsize=6.5)
		ax1.set_title('Number of Flagged Galaxies= '+numnot, fontsize=9)
		plt.ylim(min(clumpnot),max(clumpnot))
		plt.tick_params(axis='both', which='major', labelsize=5)
		ax1.invert_yaxis()
		ax1.plot(0,0, label='Redshift is '+zrange, marker='', c='k')
		#ax1.plot(0,0, label=labs[2], marker='', c='k')
		#ax1.plot(0,0,label='Number of Flagged Galaxies= '+numnot, marker='', c='k')
		ax1.legend(loc=7,prop={'size':5}, numpoints = 1)

	#plt.show()
	f.savefig(outdir+band+'_'+labs[3]+'_zclump.pdf')

	#single plot
	fig=plt.figure()
	plt.plot(aper, clumpflag, c='r', marker='^')
	if nnot==0:
		print('No Flagged Galaxies Here!')
	else:
		plt.plot(aper, clumpnot, c='b', marker='*')
	plt.xlabel('Aperture Radius (arcseconds)', fontsize=10)
	plt.ylabel('Mean Surface Brightness (mag/arcsec^2)', fontsize=10)
	plt.suptitle("Mean Surface Brightness vs. Aperture Radius in "+band, fontsize=15)
	plt.title("Redshift range "+ zrange, fontsize=12)
	red_patch = mpatches.Patch(color='red', label=labs[1] +'('+numflag+')')
	blue_patch = mpatches.Patch(color='blue', label=labs[2] +'('+numnot+')')
	plt.legend(handles=[red_patch, blue_patch], loc=7,prop={'size':5})
	plt.gca().invert_yaxis()
	#plt.show()
	fig.savefig(outdir+band+'combo_'+labs[3]+'_combozclump.pdf')


def TF_meanSB(dataflag, datanot, aper, band, magname, labs,zrange, outdir):
	import matplotlib.pyplot as plt
	import numpy as np
	import math
	import matplotlib.patches as mpatches
	pi=math.pi
	way1=1
	if not zrange:
		zrange='None Specified'
	else:
		zrange=str(zrange)
	nflag=len(dataflag)
	numflag=str(nflag) #string version
	nnot=len(datanot)
	numnot=str(nnot) #string version
	clumpflag=[]
	clumpnot=[]
	Naps=len(aper)
	
	if magname is 'mag_cmodel':
		tag='cmodel'
		print(tag)
	elif magname is'mag_aperture0':
		tag='aper'
		print(tag)
	else:
		tag=''
		print('no tag')
	
	stdvflag=[]
	stdvnot=[]
	for a in range(0, Naps):
		js=str(a)
		SB10flag=[]
		for n in range(0,nflag):
			magclump=dataflag[band+magname+js]
			sbclump=magclump+2.5*math.log10(4*pi*aper[a]**2)
			SB10=math.pow(10,sbclump[n]) #getting 10^SB
			SB10flag.append(SB10)
		SB10sum=np.sum(SB10flag)
		#print('SB10sum is ', SB10sum)
		SB10avg=SB10sum/nflag
		#print('SB10avg is ', SB10avg)
		sbavg=math.log10(SB10avg)
		clumpflag.append(sbavg)
	#making stdv
		
		if way1:
			magflag=dataflag[band+magname+js]
			SBflag=magflag+2.5*math.log10(4*pi*aper[a]**2)
			square10flag=(SBflag- sbavg)**2
			sum10sqfl=np.sum(square10flag)
			avg10fl=sum10sqfl/nflag
			sqrtfl=math.sqrt(avg10fl)
			stdvflag.append(sqrtfl)
		else:
			square10flag=(SB10flag- SB10avg)**2
			print('SB10flag is ', SB10flag)
			print('That squared with ',SB10avg,' is ', square10flag)
			sum10sqfl=np.sum(square10flag)
			print('sum of squares ', sum10sqfl)
			avg10fl=sum10sqfl/nflag
			print('average10 ', avg10fl)
			sqrtfl=math.sqrt(avg10fl)
			print('square root is ', sqrtfl)
			stdv1=math.log10(sqrtfl)
			print('standard deviation for first aperture? ', stdv1)
			#sbst=math.sqrt(np.mean(SBflag-sbavg)**2)
			stdvflag.append(stdv1)
		

		if nnot==0:
			clumpnot=np.zeros(Naps)
		else:
			SB10not=[]
			for n in range(0,nnot):
				magclump=datanot[band+magname+js]
				sbclump=magclump+2.5*math.log10(4*pi*aper[a]**2)
				SB10=math.pow(10,sbclump[n]) #getting 10^SB
				SB10not.append(SB10)
			SB10notsum=np.sum(SB10not)
			#print('SB10sum is ', SB10sum)
			SB10notavg=SB10notsum/nnot
			#print('SB10avg is ', SB10avg)
			sbnotavg=math.log10(SB10notavg)
			clumpnot.append(sbnotavg)
			if way1:
				magnot=datanot[band+magname+js]
				SBnot=magnot+2.5*math.log10(4*pi*aper[a]**2)
				square10not=(SBnot- sbnotavg)**2
				sum10sqnt=np.sum(square10not)
				avg10nt=sum10sqnt/nnot
				sqrtnt=math.sqrt(avg10nt)
				stdvnot.append(sqrtnt)
			else:
				square10not=(SB10not- SB10avg)**2	
				sum10sqnt=np.sum(square10not)
				print('sum of squares ', sum10sqnt)
				avg10nt=sum10sqnt/nnot
				print('average10 ', avg10nt)
				sqrtnt=math.sqrt(avg10nt)
				print('square root is ', sqrtnt)
				stdv2=math.log10(sqrtnt)
				print('standard deviation for first aperture? ', stdv2)
				#sbst=math.sqrt(np.mean(SBnot-sbavg)**2)
				stdvnot.append(stdv2)

	print('no flag stdv ',stdvflag, 'flagged stdv ', stdvnot)
#single plot
	fig=plt.figure()
	plt.plot(aper, clumpflag, c='r', marker='^')
	plt.errorbar(aper,clumpflag,stdvflag, c='m')
	if nnot==0:
		print('No Flagged Galaxies Here!')
	else:
		plt.plot(aper, clumpnot, c='b', marker='*')
		plt.errorbar(aper, clumpnot, stdvflag, c='c')
	plt.xlabel('Aperture Radius (arcseconds)', fontsize=10)
	plt.ylabel('Mean Surface Brightness (mag/arcsec^2)', fontsize=10)
	plt.suptitle("Mean Surface Brightness vs. Aperture Radius in "+band, fontsize=15)
	plt.title("Redshift range "+ zrange, fontsize=12)
	red_patch = mpatches.Patch(color='red', label=labs[1] +'('+numflag+')')
	blue_patch = mpatches.Patch(color='blue', label=labs[2] +'('+numnot+')')
	plt.legend(handles=[red_patch, blue_patch], loc=7,prop={'size':5})
	plt.gca().invert_yaxis()
	plt.show()
	fig.savefig(outdir+band+'_'+labs[3]+'+err_TFmeanSB.pdf')
	
	
	f=plt.figure()
	plt.plot(aper, clumpflag, c='r', marker='^')
	if nnot==0:
		print('No Flagged Galaxies Here!')
	else:
		plt.plot(aper, clumpnot, c='b', marker='*')
	plt.xlabel('Aperture Radius (arcseconds)', fontsize=10)
	plt.ylabel('Mean Surface Brightness (mag/arcsec^2)', fontsize=10)
	plt.suptitle("Mean Surface Brightness vs. Aperture Radius in "+band, fontsize=15)
	plt.title("Redshift range "+ zrange, fontsize=12)
	red_patch = mpatches.Patch(color='red', label=labs[1] +'('+numflag+')')
	blue_patch = mpatches.Patch(color='blue', label=labs[2] +'('+numnot+')')
	plt.legend(handles=[red_patch, blue_patch], loc=7,prop={'size':5})
	plt.gca().invert_yaxis()
	plt.show()
	f.savefig(outdir+band+'_'+labs[3]+'_TFmeanSB.pdf')