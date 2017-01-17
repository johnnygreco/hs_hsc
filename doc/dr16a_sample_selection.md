# Selection of Massive Galaxies in DR1/S16A 

* Song Huang / 2016-12

---------

## Bright Galaxies

* SQL searches of galaxies with 
    - Either cModel or Kron magnitude in i-band brighter than 22.5 magnitude; 
    - `ncountinput_i >= 3`; do not consider other filters now;
    - Primary detection; no child in deblending; extended object;
    - Central pixel in g,r, and i-band is not bothered by saturation, cosmic-ray, 
      interpolation, or bad pixel;
    - No pixel in i-band is off the image 

* Will do quality check later

* Files: 
    - `s16a_wide_g09/g15/w12/vvd/xmm/hec_i22.5.fits`
    - Combine them into: `s16a_wide_galaxy_i22.5.fits`

---------

## The redMaPPer catalog 

* Switch to the `v6.3` now:
    - Slightly smaller sample size; some of the clusters selected in `v5.10` are excluded
    - Small differences in `Z_LAMBDA` and `LAMBDA` values

* Prepare the new central and member galaxy catalog:
    - Central: `redmapper_dr8_public_v6.3_catalog_flat.fits`
    - Member: `redmapper_dr8_public_v6.3_members_combined.fits`
    - The format is very similar to the earlier one 

-------- 

## The CAMIRA catalog 

* Using the `S16A_v2` version: `camira_mem_s16a_wide_v2`
    - Fits version: `camira_s16a_wide_v2.fits`
    - Only takes the clusters between redshift 0.2 and 0.5: `camira_s16a_wide_v2_z0.2-0.5.fits`

------- 

## Catalog matching 

* Match the `s16a_wide_galaxy_i22.5.fits` with the redMaPPer catalogs using 1.5 arcsec radius:
    - Central: `dr16a_redbcg_wide_i22.5_v6.3.fits`
        * 764 matched ones; 644 between $0.2\<z\<0.5$; 489 between $0.3\<z\<0.5$
        * For the $0.3\<z\<0.5$ ones: 
            - 1% affected by bad cModel;
            - 1-6% affected by EDGE pixel 
            - 2-11% affected by saturation (ANY)
            - 2-5% affected by interpolation (ANY)
            - 4-23% affected by clipped pixel (ANY)
            - 25-26% affected by bright object mask (CENTER); 33-37% (ANY)
    - Member: `dr16a_redmem_wide_i22.5_v6.3.fits`
        * 49844 matches

* Match the `s16a_wide_galaxy_i22.5.fits` with the CAMIRA catalog using 1 arcsec radius: 
    - Results: `dr16a_camira_wide_i22.5_v2.fits` 
    - 1361 matches 
