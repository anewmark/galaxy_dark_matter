import matplotlib
import numpy as np
import math
import matplotlib.pyplot as plt


##making functions for the transit document

def objarr(filename):		#just returns Flux and Time as different arrays
	data=np.genfromtxt('/Users/amandanewmark/repositories/usrp-sciprog/projects/transit/data/{}'.format(filename), dtype=None, names=True)
    return data['Time'],data['Flux']

Time, Flux = objarr("7016.01.txt")	assigns Time and Flux from objarr

def TFplot(x, y):		#plots Time v Flux
    plt.plot(x, y)
    plt.xlabel('Time(s)')
    plt.ylabel('Flux(W)')
    plt.show
    return plt.plot(x,y)


TF=TFplot(Time, Flux)
N=len(Time)


Pars=[]	#depth, duration, ingress duration, and center time

def trapfit(params, x,y): #depth, duration, ingress duration, and center time
    delta=params[0]
    T=params[1]
    tau=params[2]
    t0=params[3]
    y.where(x=False)