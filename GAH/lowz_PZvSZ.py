# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 11:09:02 2016

@author: elinor
"""

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('presentation')
from astropy.coordinates import SkyCoord
from astropy import units as u
import astropy.table as table 
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=100, Om0=0.27,)
import useful_el as usel
from SimpleTimer import TicToc
import astropy.coordinates as coord

#%%

#LOWZS = table.Table.read('/Users/elinor/Dropbox/Projects/CoolCore/DR12/galaxy_DR12v5_LOWZ_South.fits')
#hscmags = table.Table.read('//Users/elinor/Data/HSC_local/GAH/HSC_XMM_apmags_s15b.fits')

#matched within topcat
lowztable = table.Table.read('/Users/elinor/Data/HSC_local/GAH/LOWZ_HSCXMM_apmags+pz.fits')

f=plt.figure()
plt.subplots_adjust(hspace=0.001,bottom=0.1)
ax1=plt.subplot(221)
plt.plot(lowztable['Z'],lowztable['demp_photoz_median'],'k.',label='demp')
plt.plot([0,0.7],[0,0.7],'--k')
plt.axis([0,2,0,2])
plt.ylabel(r'$z_{demp}$')
ax2=plt.subplot(222)
plt.plot(lowztable['Z'],lowztable['mizuki_photoz_median'],'m.',label='mizuki')
plt.plot([0,2],[0,0.7],'--k')
plt.ylabel(r'$z_{mizuki}$')
plt.axis([0,2,0,2])
ax3=plt.subplot(223)
plt.plot(lowztable['Z'],lowztable['mlz_photoz_median'],'b.',label='mlz')
plt.plot([0,2],[0,0.7],'--k')
plt.xlabel(r'$z_{spec}$')
plt.ylabel(r'$z_{mlz}$')
plt.axis([0,2,0,2])
ax4=plt.subplot(224)
plt.plot(lowztable['Z'],lowztable['nnpz_photoz_median'],'r.',label='nnpz')
plt.plot([0,2],[0,0.7],'--k')
plt.xlabel(r'$z_{spec}$')
plt.ylabel(r'$z_{nnpz}$')
plt.axis([0,2,0,2])

f.savefig('LOWZXMM_PZmedvSZ.pdf')
