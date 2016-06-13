SELECT
  meas.object_id
, meas.parent_id
, meas.ira
, meas.idec

, meas.iflux_kron
, meas.iflux_kron_err
, meas.iflux_kron_flags
, meas.imag_kron
, meas.imag_kron_err

-- unforced cmodel is not loaded in meas table.
-- Get it from mosaic_measphoto__deepcoadd and mosaic_measflag_i__deepcoadd.
, measphoto.mag_cmodel         as imag_cmodel
, measphoto.mag_cmodel_err     as imag_cmodel_err
, measflag_i.cmodel_flux_flags as iflux_cmodel_flags
, measphoto.flux_cmodel        as iflux_cmodel
, measphoto.flux_cmodel_err    as iflux_cmodel_err


-- aperture photometry
,	phot.gmag_aperture00, phot.rmag_aperture00, phot.imag_aperture00, phot.zmag_aperture00, phot.ymag_aperture00, phot.gmag_aperture00_err, phot.rmag_aperture00_err, phot.imag_aperture00_err, phot.zmag_aperture00_err, phot.ymag_aperture00_err, 
	phot.gmag_aperture01, phot.rmag_aperture01, phot.imag_aperture01, phot.zmag_aperture01, phot.ymag_aperture01, phot.gmag_aperture01_err, phot.rmag_aperture01_err, phot.imag_aperture01_err, phot.zmag_aperture01_err, phot.ymag_aperture01_err, 
	phot.gmag_aperture02, phot.rmag_aperture02, phot.imag_aperture02, phot.zmag_aperture02, phot.ymag_aperture02, phot.gmag_aperture02_err, phot.rmag_aperture02_err, phot.imag_aperture02_err, phot.zmag_aperture02_err, phot.ymag_aperture02_err, 
	phot.gmag_aperture03, phot.rmag_aperture03, phot.imag_aperture03, phot.zmag_aperture03, phot.ymag_aperture03, phot.gmag_aperture03_err, phot.rmag_aperture03_err, phot.imag_aperture03_err, phot.zmag_aperture03_err, phot.ymag_aperture03_err, 
	phot.gmag_aperture04, phot.rmag_aperture04, phot.imag_aperture04, phot.zmag_aperture04, phot.ymag_aperture04, phot.gmag_aperture04_err, phot.rmag_aperture04_err, phot.imag_aperture04_err, phot.zmag_aperture04_err, phot.ymag_aperture04_err, 
	phot.gmag_aperture05, phot.rmag_aperture05, phot.imag_aperture05, phot.zmag_aperture05, phot.ymag_aperture05, phot.gmag_aperture05_err, phot.rmag_aperture05_err, phot.imag_aperture05_err, phot.zmag_aperture05_err, phot.ymag_aperture05_err, 
	phot.gmag_aperture06, phot.rmag_aperture06, phot.imag_aperture06, phot.zmag_aperture06, phot.ymag_aperture06, phot.gmag_aperture06_err, phot.rmag_aperture06_err, phot.imag_aperture06_err, phot.zmag_aperture06_err, phot.ymag_aperture06_err, 
	phot.gmag_aperture07, phot.rmag_aperture07, phot.imag_aperture07, phot.zmag_aperture07, phot.ymag_aperture07, phot.gmag_aperture07_err, phot.rmag_aperture07_err, phot.imag_aperture07_err, phot.zmag_aperture07_err, phot.ymag_aperture07_err, 
	phot.gmag_aperture08, phot.rmag_aperture08, phot.imag_aperture08, phot.zmag_aperture08, phot.ymag_aperture08, phot.gmag_aperture08_err, phot.rmag_aperture08_err, phot.imag_aperture08_err, phot.zmag_aperture08_err, phot.ymag_aperture08_err, 
	phot.gmag_aperture09, phot.rmag_aperture09, phot.imag_aperture09, phot.zmag_aperture09, phot.ymag_aperture09, phot.gmag_aperture09_err, phot.rmag_aperture09_err, phot.imag_aperture09_err, phot.zmag_aperture09_err, phot.ymag_aperture09_err 


-- columns which can be used for selection
, meas.tract
, meas.merge_peak_g
, meas.merge_peak_r
, meas.merge_peak_i
, meas.merge_peak_z
, meas.merge_peak_y
, meas.icountInputs
, meas.ideblend_has_stray_flux
, meas.iflags_pixel_bright_object_center
, meas.iflags_pixel_bright_object_any
, meas.iblendedness_flags
, meas.iblendedness_flags_noCentroid
, meas.iblendedness_flags_noShape
, meas.iblendedness_old
, meas.iblendedness_raw_flux
, meas.iblendedness_raw_flux_child
, meas.iblendedness_raw_flux_parent
, meas.iblendedness_abs_flux
, meas.iblendedness_abs_flux_child
, meas.iblendedness_abs_flux_parent
, meas.iblendedness_raw_shape_child[1]  as iblendedness_raw_shape_child_ixx
, meas.iblendedness_raw_shape_child[2]  as iblendedness_raw_shape_child_iyy
, meas.iblendedness_raw_shape_child[3]  as iblendedness_raw_shape_child_ixy
, meas.iblendedness_raw_shape_parent[1] as iblendedness_raw_shape_parent_ixx
, meas.iblendedness_raw_shape_parent[2] as iblendedness_raw_shape_parent_iyy
, meas.iblendedness_raw_shape_parent[3] as iblendedness_raw_shape_parent_ixy
, meas.iblendedness_abs_shape_child[1]  as iblendedness_abs_shape_child_ixx
, meas.iblendedness_abs_shape_child[2]  as iblendedness_abs_shape_child_iyy
, meas.iblendedness_abs_shape_child[3]  as iblendedness_abs_shape_child_ixy
, meas.iblendedness_abs_shape_parent[1] as iblendedness_abs_shape_parent_ixx
, meas.iblendedness_abs_shape_parent[2] as iblendedness_abs_shape_parent_iyy
, meas.iblendedness_abs_shape_parent[3] as iblendedness_abs_shape_parent_ixy

FROM

-- The meas table is split into "meas", "meas_centroid", and "meas_shape"
-- because PostgreSQL forbids a table more than ~1600 columns

s15b_wide.meas as meas
LEFT JOIN s15b_wide.meas as shape USING (object_id)
LEFT JOIN s15b_wide.meas as centroid USING (object_id)

LEFT JOIN s15b_wide.mosaic_measphoto__deepcoadd as measphoto on (meas.object_id = measphoto.id)
LEFT JOIN s15b_wide.mosaic_measflag_i__deepcoadd as measflag_i on (meas.object_id = measflag_i.id)
LEFT JOIN s15b_wide.photoobj_mosaic__deepcoadd__merged as phot on (meas.object_id = phot.id)  -- aperture photometry catalog
 
WHERE
-- Please uncomment to get a field you want

-- AEGIS
-- s15b_search_aegis(meas.patch_id)           AND

-- HECTOMAP
-- s15b_search_hectomap(meas.patch_id)        AND

-- GAMA09H
-- s15b_search_gama09h(meas.patch_id)         AND

-- WIDE12H
-- s15b_search_wide12h(meas.patch_id)         AND

-- GAMA15H
-- s15b_search_gama15h(meas.patch_id)         AND

-- VVDS
-- s15b_search_vvds(meas.patch_id)            AND

-- XMM
 s15b_search_xmm(meas.patch_id)             AND

 NOT meas.ideblend_skipped                  AND
 NOT meas.iflags_badcentroid                AND
 NOT centroid.icentroid_sdss_flags          AND
 NOT meas.iflags_pixel_edge                 AND
 NOT meas.iflags_pixel_interpolated_center  AND
 NOT meas.iflags_pixel_saturated_center     AND
 NOT meas.iflags_pixel_cr_center            AND
 NOT meas.iflags_pixel_bad                  AND
 NOT meas.iflags_pixel_suspect_center       AND
 NOT meas.iflags_pixel_clipped_any          AND
 meas.idetect_is_primary             	    AND
 NOT shape.ishape_hsm_regauss_flags    	    AND
 meas.iclassification_extendedness != 0     AND
 -- In postgres, all numbers are comparable including NaN
 shape.ishape_hsm_regauss_sigma != 'NaN'    AND
 measphoto.filter01 = 'HSC-I'
 ORDER BY meas.object_id
;
