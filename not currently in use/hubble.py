
import matplotlib
import numpy as np
import math
import matplotlib.pyplot as plt

import sys
sys.path.append("/Users/amandanewmark/repositories/usrp-sciprog/projects/hubble/")

data=np.genfromtxt('/Users/amandanewmark/repositories/usrp-sciprog/projects/hubble/table1.txt', dtype=None, names=True)

N=len(data)
A=np.empty((N,2))
A[:,0]=data['R']
A[:,1]=1
params, _, _, _=np.linalg.lstsq(A, data['V'])
H0=params[0]

R=np.linspace(0,2.5,100)
fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(R, H0*R, 'k--')
ax.scatter(data["R"], data['V'],c='#ff9900')
ax.set_xlabel('Distance [MPC]')
ax.set_ylabel('Velocity[km/s]')
ax.set_xlim(xmin=0, xmax=2.5)
plt.show()

from coordinates import*
A=np.empty((N,4))
ra=Ra2Deg(data['RA'])
dec=Ra2Deg(data['DEC'])
A[:,0]=data['R']
A[:,1]=np.cos(ra*np.pi/180)*np.cos(dec*np.pi/180)
A[:,2]=np.sin(ra*np.pi/180)*np.cos(dec*np.pi/180)
A[:,3]=np.sin(dec*np.pi/180)