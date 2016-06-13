print('This program should make plots of aperture size v aperture luminosity')

indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'

outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/'

import numpy as np
import math
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.table as table 

lowztable = table.Table.read(indir+ 'LOWZ_HSCGAMA15_apmgs.fits')

#lowztable=fits.open(indir+ 'LOWZ_HSCGAMA15_apmgs.fits') <-- probably better to use table
print('This worked')
hdulist=fits.open(indir+ 'LOWZ_HSCGAMA15_apmgs.fits')
hdulist.info()
#print(lowztable['Z'])
#lowztable['Z'] gives Z shift
#