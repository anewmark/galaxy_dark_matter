from defcuts import *
from defflags import *
from defclump import *
import astropy.table as table 


indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'

outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'

datatab = table.Table.read(indir+ 'LOWZ_HSCGAMA15_apmgs.fits')

band='i'

parm=['flags_pixel_saturated_center','flags_pixel_edge','flags_pixel_interpolated_center','flags_pixel_cr_center','flags_pixel_suspect_center', 'flags_pixel_clipped_any','flags_pixel_bad']	

Flags1=['flags_pixel_bright_object_center', 'No Flags', 'No Bright Ojbect Centers', 'Only Bright Object Centers', 'brobj_cen_flag']

Flags2=['flags_pixel_bright_object_any', 'No Flags', 'No Bright Ojbects', 'Only Bright Objects', 'brobj_any_flag']

aperture=[3.0*0.168,4.5*0.168,6.0*0.168,9.0*0.168,12.0*0.168,17.0*0.168,25.0*0.168,35.0*0.168,50.0*0.168,70.0*0.168]

mincut= 17.5
maxcut=18.5

minz=0.25
maxz=0.35

colname='mag_aperture0'
colnames='mag_cmodel'

datazcut, rangez=z_cut(datatab, minz, maxz)

cutdata=out_cut(datazcut, band, colnames,mincut, maxcut)
	#cutting the data according to cmodel, plotting according to apertures
newdata=many_flags(cutdata, parm, band)	#this gets rid of multiple flags

Flag, Not,lab= TFflag(band,Flags1, newdata)

TF_meanSB(Flag, Not, aperture, band, colname, lab,rangez, outdir)

