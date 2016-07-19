from defcuts import *
from defflags import *
from defclump import *
import astropy.table as table 


indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'

outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'

datatab = table.Table.read(indir+ 'LOWZ_HSCGAMA15_apmgs.fits')

band='i'

parm=['flags_pixel_saturated_center','flags_pixel_edge','flags_pixel_interpolated_center','flags_pixel_cr_center','flags_pixel_suspect_center', 'flags_pixel_clipped_any','flags_pixel_bad']	

Flags=['flags_pixel_bright_object_center', 'No Flags', 'No Bright Ojbect Centers', 'Only Bright Object Centers', 'brobj_cen_flag']

#Flags=['flags_pixel_bright_object_any', 'No Flags', 'No Bright Ojbects', 'Only Bright Objects', 'brobj_any_flag']

#Flags=['blendedness_flags', 'No Flags', 'Not Blended', 'Blended', 'blend']

daperture=[1.01,1.51,2.02,3.02,4.03,5.71,8.40,11.8,16.8,23.5]

aperture=[x*0.5 for x in daperture]

mincut= 17.5
maxcut=18.5

minz=0.25
maxz=0.35

colname='mag_aperture0'
cutcolname='mag_aperture05'
#cutcolname='mag_cmodel'

datazcut, rangez=z_cut(datatab, minz, maxz)

cutdata, crange=out_cut(datazcut, band, cutcolname,mincut, maxcut)

newdata=many_flags(cutdata, parm, band)	#this gets rid of multiple flags

Flag, Not,lab= TFflag(band,Flags, newdata)

TF_meanSB(Flag, Not, aperture, band, colname, lab,rangez, outdir)

