print('Please for the love of God work.')
import os, sys
import astropy.table as table
from defclump import *
from defcuts import *
from defflags import *

indir='/Users/amandanewmark/repositories/galaxy_dark_matter/GAH/'
filename='LOWZ_HSCGAMA15_apmgs.fits'

outdir='/Users/amandanewmark/repositories/galaxy_dark_matter/lumprofplots/clumps/'
if not os.path.exists(outdir):
    os.mkdir(outdir)

 
datatab = table.Table.read(indir+ filename)
objID=datatab['object_id']



aperture=[3.0*0.168,4.5*0.168,6.0*0.168,9.0*0.168,12.0*0.168,17.0*0.168,25.0*0.168,35.0*0.168,50.0*0.168,70.0*0.168]	#these are the apertures

bands=['r', 'i']
nbands=len(bands)
flags=['flags_pixel_saturated_center', 'flags_pixel_edge', 'flags_pixel_interpolated_center', 'flags_pixel_cr_center', 'flags_pixel_bad', 'flags_pixel_suspect_center', 'flags_pixel_clipped_any']
nflags=len(flags)

print(nbands)
for b in range(0,nbands):
	band=bands[b]
	#band='r'
	for f in range(0, nflags):	#This gets rid of single flags
		flag=flags[f]
		fs=str(f)
		print(flags[f], band)
		
		if flag==flags[0]:
			parms=['flags_pixel_saturated_center', 'No Flags', 'No Saturated Galaxies', 'Only Saturated Galaxies', 'sat_flag']	#these are the flag parameters
		elif flag==flags[1]:
			parms=['flags_pixel_edge', 'No Flags', 'No Edge Galaxies', 'Only Edge Galaxies', 'edge_flag']
		elif flag==flags[2] :
			parms=['flags_pixel_interpolated_center', 'No Flags', 'No Interpolated Centers', 'Only Interpolated Centers', 'interp_flag']
		elif flag==flags[3]:
			parms=['flags_pixel_cr_center', 'No Flags', 'No Cr Centers', 'Only Cr Centers', 'crcen_flag']
		elif flag==flags[4]:
			parms=['flags_pixel_bad', 'No Flags', 'No Bad Pixels', 'Only Bad Pixels', 'badpix_flag']
		elif flag==flags[5]:
			parms=['flags_pixel_suspect_center', 'No Flags', 'No Suspect Centers', 'Only Suspect Centers', 'sus_flag']
		elif flag==flags[6]:
			parms=['flags_pixel_clipped_any', 'No Flats', 'No Clipped Galaxies', 'Only Clipped Galaxies', 'clip_flag']
		else:
			print('woops')
	
		#parms=['flags_pixel_saturated_center', 'No Flags', 'No Saturated Galaxies', 'Only Saturated Galaxies', 'sat_flag']	#these are the flag parameters

		datazcut, rangez=z_cut(datatab, 0.15, 0.25) #give it data table, zmin, zmax

		Flag, Not,lab= TFflag(band,parms, datazcut)

		clump_plot(aperture, Flag, Not, band, 'mag_aperture0', lab, rangez, outdir)
			#note: for whatever reason, y axis not inversed in two plots
		#break
		#create new twoflags and three z cut functions