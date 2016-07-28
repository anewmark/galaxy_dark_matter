import matplotlib as mpl
import matplotlib.pyplot as plt
def get_presentation():
	#will hopefully get presentation parameters for plt
	mpl.rcParams['axes.titlesize']=22
	mpl.rcParams['axes.labelsize']=18
	mpl.rcParams['lines.linewidth']=2
	mpl.rcParams['lines.markersize']=10
	mpl.rcParams['xtick.labelsize']=16
	mpl.rcParams['ytick.labelsize']=16
	mpl.rcParams['legend.borderpad']=0.15
	mpl.rcParams['legend.borderaxespad']=0.0
	mpl.rcParams['legend.labelspacing']=0.15