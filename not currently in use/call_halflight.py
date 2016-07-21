print('Finding Half Light Radii and Corresponding Slopes')

import astropy.table as table 
import numpy as np
import matplotlib.pyplot as plt
from defcuts import *
from def_get_mags import *
from def_lump_prof import *
from my_def_plots import *
from defflags import *

indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'

datatab = table.Table.read(indir+ 'LOWZ_HSCGAMA15_apmgs+cmodmag.fits')

bands=['g', 'r', 'i','z', 'y']

parm=['flags_pixel_saturated_center','flags_pixel_edge','flags_pixel_interpolated_center','flags_pixel_cr_center','flags_pixel_suspect_center', 'flags_pixel_clipped_any','flags_pixel_bad']

Flags=['flags_pixel_bright_object_center', 'brobj_cen_flag-', 'No Bright Ojbect Centers', 'Only Bright Object Centers', 'brobj_cen_flag']

daperture=[1.01,1.51,2.02,3.02,4.03,5.71,8.40,11.8,16.8,23.5]

aperture=[x*0.5 for x in daperture]

#get rid of cuts
#mag_cmodel
ne=[99.99, 199.99, 0.0]

mincut=0.1
maxcut=''
cutdata=not_cut(datatab, bands, 'mag_forced_cmodel', ne)

#get rid of flagged galaxies
for b in range(0, len(bands)-1):
	newdata=many_flags(cutdata, parm, bands[b])	#flags not in y?
	cutdata=newdata

bandi=['i']

outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/single_plot/'

Flagdat, Notdat,lab= TFflag(bandi,Flags, newdata)

halfrad, rad, lumdens=get_halflight(Flagdat, bands, aperture, scale='log')

objid=Flagdat['object_id_1']

m1=halflight_slopes(halfrad, rad, lumdens, objid, plots='yes', outdir=outdir)

#*****now for Not:

halfradN, radN, lumdensN=get_halflight(Notdat, bands, aperture, scale='log')

m2=halflight_slopes(halfradN, radN, lumdensN, objid, plots='yes', outdir=outdir)

figs=plt.figure()
bs=np.linspace(-2.0,-1.4,num=15, endpoint=False)
plt.hist(m1, bs, color='red', label='Not Flagged Galaxies', alpha=0.8)
plt.hist(m2,bs, color='blue', label='Flagged Galaxies', alpha=0.8)
plt.xlabel('Slopes', fontsize=10)
plt.legend(loc=0,prop={'size':7.0})
plt.ylabel('Frequency', fontsize=10)
plt.show()
