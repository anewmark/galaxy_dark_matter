from defcuts import *
from defflags import *
from defdist import *
import astropy.table as table 


indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'

outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/single_plot/'

datatab = table.Table.read(indir+ 'LOWZ_HSCGAMA15_apmgs.fits')

band='i'

parm=['flags_pixel_saturated_center','flags_pixel_edge','flags_pixel_interpolated_center','flags_pixel_cr_center','flags_pixel_suspect_center', 'flags_pixel_clipped_any','flags_pixel_bad']	

Flags1=['flags_pixel_bright_object_center', 'No Flags', 'No Bright Ojbect Centers', 'Only Bright Object Centers', 'brobj_cen_flag']

Flags2=['flags_pixel_bright_object_any', 'No Flags', 'No Bright Ojbects', 'Only Bright Objects', 'brobj_any_flag']

mincut= 0.001
maxcut=''

zmin=0.2
zmax=0.3

#colname='mag_aperture00'
colname='mag_cmodel'
zcolname='Z'

cutdata=out_cut(datatab, band, colname,mincut, maxcut)


newdata=many_flags(cutdata, parm, band)	#this gets rid of multiple flags

#for bright_center
Flag, Not,lab= TFflag(band,Flags1, newdata)

hist_mag(Flag, Not, lab, band, colname, outdir)


#for bright_any

Flag, Not,lab= TFflag(band,Flags2, newdata)

hist_mag(Flag, Not, lab, band, colname, outdir)