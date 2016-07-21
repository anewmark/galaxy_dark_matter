print('Will get absolute magnitudes from this')

import astropy.table as table 
import numpy as np
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

#Flags=['flags_pixel_bright_object_any', 'brobj_any_flag-', 'No Bright Ojbects', 'Only Bright Objects', 'brobj_any_flag']

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
Flagdat, Notdat,lab= TFflag(bandi,Flags, newdata)

#error='stdv'
error='bootstrap_stdv'
outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/flags/'


Li1, rad1, mean1, error1, radavg1= lump_prof2(Flagdat, bands, aperture, err=error)
Li2, rad2, mean2, error2, radavg2= lump_prof2(Notdat, bands, aperture,err=error)

#rads[n] reproduces each row with 10 apertures, same with Li1[n]

lum_overplot(Li1, Li2,rad1, rad2, mean1, mean2, error1, error2, radavg1, radavg2, lab, outdir, error)

print('error no flag= ', error1)
print('error only flag= ', error2)