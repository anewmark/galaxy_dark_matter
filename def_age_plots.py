import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import math
outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'

def age_plot(x,y, start, end, tag=[]):
	#print('x spacing: ', xspace)
	print('x_center values: ', x)
	print('y values= ', y)
	start=np.array(start)
	end=np.array(end)
	arange=end-start
	arange=arange[1:]
	
	
	bwidth=arange
	print('bar width=', bwidth)
	print(len(bwidth))
	f=plt.figure()
	
	labels=np.append(start, end[len(end)-1])
	print('labels', labels)
	xspace=np.linspace(np.min(labels), np.max(labels), num=len(labels))
	#label=[str(n) for n in labels]
	
	#print(label)
	#plt.bar(x, y, width=bwidth, align='center', color='None')
	plt.bar(x, y, width=bwidth, align='center', color='None')
	plt.xlabel('Lookback Time (Gyr)')
	
	plt.xticks(labels, label=labels, fontsize=8, rotation='vertical') #smushed
	#plt.xticks(xspace, label=labels, fontsize=8, rotation='vertical') #ticks dont align with bar graph
	
	plt.ylabel('Stacked Mass Fractions')
	plt.title('Age vs. Mass Fractions')
	#plt.xlim(np.min(x)-bwidth[0], np.max(x)+bwidth[len(bwidth)-1]/2.0)
	plt.xlim(np.min(start), np.max(end))
	plt.ylim(0,1)
	plt.plot(0,0, label=tag[1], c='k')
	plt.plot(0,0, label=tag[2], c='k')
	plt.legend(loc=2,prop={'size':6.0})
	
	plt.show()
	outdirs=outdir+tag[0]+'agebin.pdf'
	f.savefig(outdirs)
	print(outdirs)
	
def age_plot1(x,y, start, end, yerr, tag=''):
	#print('x spacing: ', xspace)
	print('x_center values: ', x)
	print('y values= ', y)
	print('y errors= ', yerr)
	start=np.array(start)
	end=np.array(end)
	arange=end-start
	arange=arange

	bwidth=arange
	print('bar width=', bwidth)
	print(len(bwidth))
	f=plt.figure()
	
	labels=np.append(start, end[len(end)-1])
	#print('labels', labels)
	xspace=np.linspace(np.min(labels), np.max(labels), num=len(labels))
	label=[str(n) for n in labels]
	xsp=[math.log10(n) for n in labels]
	#print(xsp)
	
	#print(label)
	#plt.bar(x, y, width=bwidth, align='center', color='None')
	plt.bar(x, y, width=bwidth, align='center', color='None')
	plt.errorbar(x, y, yerr=yerr,label='Standard Error on the Mean per Age Bin', fmt='.', color='b')
	plt.xlabel('Lookback Time (Gyr)')
	plt.xscale('log')
	#plt.xticks(labels, label=labels, fontsize=8, rotation='vertical') #smushed
	plt.xticks(xsp, label=labels, fontsize=8, rotation='vertical')
	#plt.xticks(xspace, label=labels, fontsize=8, rotation='vertical') #ticks dont align with bar graph
	
	plt.ylabel('Stacked Mass Fractions')
	plt.title('Age vs. Mass Fractions')
	#plt.xlim(np.min(x)-bwidth[0], np.max(x)+bwidth[len(bwidth)-1]/2.0)
	plt.xlim(np.min(start), np.max(end))
	
	nmax=np.max(y)+yerr[len(yerr)-1]
	if nmax<1:
		nmax=1

	plt.ylim(0,nmax)
	plt.plot(0,0, label=tag[1], c='k')
	plt.plot(0,0, label=tag[2], c='k', marker='*')
	plt.plot(0,0, label=tag[3], c='k', marker=' ')
	plt.legend(loc=2,prop={'size':7.0})
	
	plt.show()
	outdirs=outdir+tag[0]+'agebin.pdf'
	f.savefig(outdirs)
	print(outdirs)

