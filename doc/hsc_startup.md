## HSC Survey Related Start-up Project: 

### LSST and HSC Survey software pipeline

#### Installation

* HSC survey data pipeline is actually built upon the LSST data pipeline, which 
  is a very powerful, flexible photometry data pipeline; The current stable 
  LSST software stack is V9.0 (Winter2014).

* To install the LSST pipeline, follow the instructions here: 
  https://confluence.lsstcorp.org/display/LSWUG/LSST+Software+User+Guide 
  - This website only contains many useful information regarding the basic 
    concept of the software stack, and many, many more. **SHOULD PAY ATTENTION 
    TO IT** 

* To test the installation of LSST pipeline, follow: 
  https://confluence.lsstcorp.org/display/LSWUG/Testing+the+Installation
  - May encounter an error involves the GLIB-2.15 or such; The work-around is 
    remove the libm.so.6 soft-link in $ANACONDA_DIR/lib
  - Have passed the demo.sh tests on Thinkpad

* The procedures for installing the HSC survey software stack is somewhat 
  similar to the LSST pipeline; And, the instruction can be found at: 
  http://hsca.ipmu.jp/hscsphinx/pipeline/installation.html
  - Error about curl: put the ca_bundle.crt file to the desired location 
  - Error during building Python: tell EUPS to use the Python in the system 
    * put "Python system" in ~/.eups/manifest.remap 
  - Error when compiling PSFex, replace 'plwid' with 'plwidth' in psfex/src/cplot.c 
    or remove the plplot package from your computer 

* Should also read the FAQ page on the LSST trac page: 
  - https://dev.lsstcorp.org/trac/wiki/DM/Policy/UsingDMCode/FAQ 

#### How to use the pipeline

* Both of the LSST and HSC software stacks are maintained by EUPS: 
  - https://dev.lsstcorp.org/trac/wiki/EupsManual
  - Allow users to have different versions of the same package on one computer
  - eups distrib list <pkg>
  - eups distrib install <pkg> <version> 
  - eups declare <pkg> <version> -r none -m none -L dir/to/the/<pkg>.cfg
  - --nolocks parameters can sometimes become handy

* The batch processing is handeled by a system called TORQUE on Master 

* Term "RERUN" is borrowed from SDSS, and has the same meaning 
  - The same "RERUN" value indicates the data have been reduced in exactly the 
    same manner

* Important DATAID: 
  1: "VISIT"  : Specific exposure 
  2: "CCD"    : Specific CCD in that exposure
  3: "TRACT"  : Coordinate system used in the coadded images
  4: "PATCH"  : Coordinate system used in the coadded images
  5: "FIELD"  : The "OBJECT" entry from FITS header 
  6: "DATEOBS": The "DATE-OBS" entry from FITS header 
  7: "FILTER" : The "FILTER" entry from FITS header
  8: "TASK"   : A specific pipeline module
  - It’s impossible to have the same visit with a different filter or dateObs
  - Configuration parameters have a hierarchical form

* To reduce certain data on MASTER, you can: 
  1. `hscProcessCcd.py /data/ --id visit=1000 ccd=50 filter=HSC-I dateObs=2015-01-15`
  2. `hscProcessCoadd.py /data/ --id tract=0 patch=1,1 filter=HSC-I`

* The naming of HSC image: `HSC-VVVVVVV-CCC.fits`: 
  - VVVVVVV: 7-digit VISIT number 
  - CCC    : 3-digit CCD number

* `VISIT` and `CCD`   defines a single frame observation 
* `TRACT` and `PATCH` defines a co-add image
* `DateObs` and `FILTER` applies to both

* `_mapper` file indicates which instrument is used for this run 
  - Should be `lsst.obs.hsc.HscMapper`

* The directory tree for the data repository is: 
  - $SUBARU_DIR/OBJECT/DATE-OBS/POINTING/FILTER-NAME
  - e.g. `/data/Subaru/HSC/SSP_WIDE/2014-03-25/00814/HSC-Y/`
* The directory tree for pipeline output data: 
  - `.../rerun/POITING/FILTER/[corr] or [output] or [processCcd_metadata] or
    [qa] or [thumbs]`
  - `.../rerun/_parent`
  - `.../rerun/config/`
  - `.../rerun/schema/` 
* Among the output files, the useful ones are: 
  1. Corrected frames: `corr/CORR-VVVVVVV-CCC.fits`
  2. Background frame: `corr/BKGD-VVVVVVV-CCC.fits`
  3. Measurements on sources: `output/SRC-VVVVVVV-CCC.fits`
  4. Objects matched to catalog sources: `output/MATCH-VVVVVVV-CCC.fits`
  - Besides these, in the qa[quality assurance] and thumbs folders, there are 
    some useful PNG files

* The general procedure to run hscPipe on some data is like this: 
  http://hsca.ipmu.jp/hscsphinx/pipeline/quick.html

##### Making Detrends 

* After every step of detrend reduction, should run `genCalibRegistry.py` 
  to make the new calibration file available 
* If need to restart a job, or play around with some configurations, make sure 
  `--clobber-config` parameter is added to the command

1. Biases:  `reduceBias.py`
2. Darks:   `reduceDark.py`
3. Flats:   `reduceFlat.py`
4. Fringes: `reduceFringes.py`
   - Only necessary for Y-band

##### Single Frame Process 

* See help for `reduceFrames.py` and `hscProcessCcd.py` 
  - `reduceFrames.py $DATA_DIR --rerun RERUN_NAME --id DATAID 
       --queue NAME_OF_QUEUE --job NAME_OF_JOB --nodes NUMBER_OF_NODES 
       --procs NUMBER_OF_PROCESSES_ON_EACH_NODE `
  - `hscProcessCcd.py $DATA_DIR --rerun RERUN --id DATAID 
       --config CONFIG_PARA=VALUE --clobber-config 
       -C CONFIG_PARAMETER_FILE_NAME` 

##### Coadded Frame Process

* Steps: 
  1. SkyMap making --> Build the common coordinate system 
     - `TRACT` and `PATCH` is defined in this process, 
        the former one is at larger scale
     - `makeSkyMap.py $DATA_DIR --rerun=RERUN` (Whole sky)
     - `makeDiscreteSkyMap.py $DATA_DIR --rerun=RERUN --id DATAID` (partial) 
       * Only `TRACT=0` will be created in this step
     - `showVisitSkyMap.py` can be used to display visits on SkyMap
     - In DATAID, if all CCDs are needed, it should be defined as 
       `ccd=0..103`` (The rest CCDs are FOCUS)
  2. Mosaic making --> Uber-calibration 
     - `mosaic.py $DATA_DIR --rerun=RERUN --id DATAID` 
  3. Warp   making --> Resample the images from observed WCS to common SkyMap 
     (convert CCDs into `PATCHES`)
     - Each CCD will contribute to 4 patches
     - `makeCoaddTempExp.py DATA_DIR --rerun=RERUN --id tract=xx patch=y,z 
          filter=HSC-[G/R/I/Z/Y] --selectId visit=AAAAAA ccd=BBBBBB`
  4. Coadd  making --> Combine the wraped images
     - `assembleCoadd.py DATA_DIR --rerun=RERUN --id tract=xx patch=y,z 
          filter=HSC-[G/R/I/Z/Y] --selectId visit=AAAAAA ccd=BBBBBB`
  5. Process       --> Make detections and measurements on the co-add images
     - `hscProcessCoadd.py DATA_DIR --rerun RERUN --id tract=xxx patch=y,z 
          filter=HSC-[G/R/I/Z/Y]
     - The last three steps can be done with one command: 
       * `stack.oy DATA_DIR --rerun=IN_RERUN:OUT_RERUN --id DATA_ID[TRACT,FILTER] 
            --selectId SELECT_ID --queue QUEUE_NAME --nodes NUMBER_OF_NODES 
            --job JOB_NAME --procs NUMBER_OF_PROCESSES_ON_EACH_NODE`
       * If restacking is needed, then following parameters are necessary: 
         `--output=NEW_RERUN_LOCATION --clobber-config 
            --config doOverwriteOutput=True doOverwriteCoadd=True`

##### Pipeline Tool

* **Butler**: A generic mechanism for persisting and retrieving data using
              mappers.  
  - A Butler manages a collection of datasets known as a repository. Each
    dataset has a type representing its intended usage and a location
  - Finds datasets by scientifically-meaningful key/value pairs
  - Retrieves datasets as in-memory objects
  - Persists in-memory objects to datasets

* **Calexp**: A fully-qualified image, which includes a science pixel array,
              and concomitant data including a quality mask and a variance array. 

* The main pipeline tools you will like need to be concerned with are the
  following:
    1. The `butler` and `dataRef`: These are tools which can find and
       load various types of data for you.
    2. `Exposures`, `MaskedImages`, and `Images`: These are the
       various containers used to handle images.
    3. `SourceCatalogs`: These are containers for tabulated information
       about ‘sources’, including things like coordinates, fluxes

* The basic way to use butler is like this: 

`python 
   import lsst.daf.persistence as dafPersist 
   import hsc.pipe.base.butler as hscButler 

   dataDir = "/data3b/Subaru/HSC/rerun/yasuda/SSP1" 
   butler  = dafPersis.Butler( dataDir ) 

   # For single frame data 
   dataId = { 'visit':1332, 'ccd':39 } 
   # Get the CORR image 
   calexp_img = butler.get( 'calexp', dataId ) 
   # Get the location and filename for BIAS 
   bias_file  = butler.get( 'bias_filename', dataId ) 
   # Get the metadata for PSF 
   psf_mtdata = butler.get( 'psf_md', dataId )

   # For coadded data 
   dataId = { 'tract': 9374, 'patch':'1,7' } 
   # Coadded image 
   coadd_img  = butler.get( 'deepCoadd_calexp', dataId )

   # If need to work on the data from same dataId intensively, can :
   dataRef = hscButler.getDataRef( butler, dataId ) 
   coadd_psf  = butler.get( 'deepCoadd_psf' )

`

* Run `ls $OBS_SUBARU_DIR/policy/HscMapper.paf` can find the file that 
  includes the descriptions for all Butler targets; And the most commonly used 
  ones are:
* For single frame: 
  - `BIAS`; `DARK`; `FLAT`; `FRINGE`
  - `CALEXP`, `PSF`, `SRC`
  - `postISRCCD`: after instrument signature removal (NOT WRITTEN TO DISK 
       by default, need to define `isr.doWrite=True`
* For coadd image: 
  - `deepCoadd_calexp`; `deepCoadd_psf`; `deepCoadd_src`


--------------------------------------------------------------------------------

### How to query the HSC Survey Database 

* The HSC Survey database at IPMU can be accessed through PostgreSQL queries: 
  - `psql -h hscdb.ipmu.jp -U kensaku -d dr_early -p 5432 -f xxx.sql` 
  - The password is "104chips", which is the default password for HSC survey
  - The current database for early data is called "dr_early"

* The current two schemas are: 
  1. ssp_s14a0_wide_20140523a 
  2. ssp_s14a0_udeep_20140523a


### BCGs, HSC Survey, and MaNGA Survey

#### MaNGA BCG Proposal: Brightest Cluster Galaxies and their Dark Matter Halos

* By Claire Lackner, Alexie Leauthaud, et. al. 
* To study the stellar and halo properties of BCG by: 
  1. Using the MaNGA instrument to observe a more complete sample of nearby 
     BCGs ($$z \le 0.15$$) 
  2. Combined the photometry from HSCS and the kinematic information from MaNGA 

--------------------------------------------------------------------------------

* Need to find out how many BCGs are located within the overlap region of MaNGA 
  footprints and HSC survey wide field 
  - At this point, using the "BCG" sample selected from the Yang+2007 DR7 Group 
    Catalog: The central galaxies of any group whose total halo mass is larger 
    than $$10^{13.75} M_{sun}$$ is considered as BCGs 
  - The RA/Dec definition of all wide fields can be found at: 
    http://hscsurvey.pbworks.com/w/page/47058337/Wide%20survey%20fields
  - There 3469 BCG-like galaxies in the catalog, 615 of them have 
    in_manga==1 (from a random selected test sample) 
  - A quick search using the catalog returns: 
    1. 28 BCGs in the HectoMap field: (4 have in_manga==1; 
       5 have $$grp\_Mhalo\_Lest \ge 14.0$$; 1 has $$grp\_Mhalo\_Lest \ge 14.25$$);
    2. 65 BCGs in the Stripe82 field: (23 have in\_manga==1; 
       25 have grp\_Mhalo\_Lest >= 14.0; 6 has grp\_Mhalo\_Lest>=14.25 ); 

* Make a list of BCGs in really massive halos (grp_Mhalo_Lest>14.00) that are 
  located within the HectoMap and Stripe82 region; And try to identify a few
  interesting objects at intermediate redshift
  - By checking the SDSS color images, most of the objects are indeed massive, 
    early-type, old galaxies, but a lot of them appear to be quite "isolated" 
    galaxy 
  - Alexie has selected a few interesting objects for the MaNGA project: 
  |#  | ra           | dec         | z         |grp\_Mhalo\_Lest| M\_r       |       
  |---|:------------:|:-----------:|:---------:|:--------------:|:---------:|
  |h1 | 245.362216361| 42.761290060| 0.1355748 | 14.4627        |-22.9543381|
  |h2 | 214.322165473| 43.521017585| 0.1371584 | 14.1359        |-22.2173538|       
  |h4 | 246.426331360| 43.931768058| 0.1331387 | 14.1077        |-22.1327534|        
  |s1 |  14.437667261| -0.419432998| 0.0438390 | 14.5595        |-22.2236042|
  |s2 |  13.302346094| -0.800121073| 0.1363631 | 14.5463        |-22.5665398|
  |s3 |  18.739971721|  0.430799451| 0.0448190 | 14.4605        |-22.1312408|
  |s4 |   7.368434375| -0.212632542| 0.0598650 | 14.3447        |-22.0776558|
  |s5 |  11.920603743| -0.857469044| 0.1145691 | 14.3126        |-22.4481277|  

* By querying the `ssp_s14a0_wide_20140523a.photoobj_mosaic__deepcoadd__iselect`
  catalog using `imag_kron < 22.0 mag` as constraint, the database returns 
  1185043 detections. 
  - Wtihin the existing HSC/GAMA early field, there are 15 "BCGs" from the 
    `Yang_MaNGA.fits` catalog 
  - Extract the basic information (ra, dec, redshit, MHalo), and SDSS images 
    (both JPEG and FITS ones) from Skyserver
  | id                | ra          | dec          |i_kron |y_kron| z     |Mhalo_L| Mr    |
  |:-----------------:|:-----------:|:------------:|:-----:|:----:|:-----:|:-----:|:-----:| 
  |2638574486412078503|219.432661025|-0.31603409809|16.621 |16.170|0.13793|14.481 |-22.689|  
  |2775106067913246734|218.537679870| 1.61692103095|15.925 |15.471|0.13786|14.357 |-22.600|  
  |2706409680921441283|216.568121070| 0.83762930374|15.981 |15.504|0.12499|14.191 |-22.158|  
  |2706445415049333151|215.696902169| 1.22992602790|15.504 |15.048|0.11726|14.153 |-22.600|  
  |2638600599813239067|218.946991044|-0.47711472421|15.765 |15.329|0.10568|13.996 |-22.140|
  - We will call these 5 BCGs `bcgs_sdss_1/2/3/4/5`: 
    1. bcgs_sdss_1: sdss_i=14.82; logM_pass_port=11.85; in_maxbcg; 
    2. bcgs_sdss_2: sdss_i=15.21; logM_pass_port=11.90; in_maxbcg;
    3. bcgs_sdss_3: sdss_i=15.21; logM_pass_port=11.56; in_maxbcg;
    4. bcgs_sdss_4: sdss_i=14.80; logM_pass_port=11.74; in_maxbcg;
    5. bcgs_sdss_5: sdss_i=15.01; logM_pass_port=11.66; not_in_maxbcg;
  - Besides of the SDSS single frame i-band image; Also generate the 0.1 deg 
    mosaic images in i- and z-band using the following services: 
      1. http://sdss.physics.nyu.edu/sdss3/mosaic/index.php
      2. http://dr10.sdss3.org./mosaics

* Extract snapshot HSC-I band images for these five BCGs: 
  - There DataID are (deeociadd_iselect): 
    1. bcgs_sdss_1: tract=9374, patch=3,6
    2. bcgs_sdss_2: tract=9859, patch=5,1
    3. bcgs_sdss_3: tract=9615, patch=3,5
    4. bcgs_sdss_4: tract=9615, patch=7,7
    5. bcgs_sdss_5: tract=9374, patch=6,5

--------------------------------------------------------------------------------

* We also need to find out how many BCGs are there within the current HSC survey 
  data ? 
  - At this point, about 20 square degree field within the GAMA field has been 
    observed using r and i-band filters; The data should already reach 26 mag 
    per square arcsec in i-band. 
  - Now we have two BCG catalog to cross-match: 
    1. The Yang+2007 Group Catalog: (Again, z<0.15; log(M_halo)>13.75) 
       * According to Claire, and a sutdent from Priceton, about 80% of these 
         brightest galaxy in a massive group can be seen as genuine BCGs
    2. RedMapper catalog for SDSS DR8: 
       * RedMapper is a very sophisticated, "red-sequence"-based, cluster finder 
         - http://risa.stanford.edu/redmapper/ 
         - http://adsabs.harvard.edu/abs/2014ApJ...785..104R 
         - http://adsabs.harvard.edu/abs/2014ApJ...783...80R 

--------------------------------------------------------------------------------

#### Basic Information about MaNGA Survey

* Should always refer to the SDSS-IV Trac website 
  - Has already required the approval to log-in

--------------------------------------------------------------------------------

### HSC SSP Tables 

#### Summary of useful catalogs 

1.  `calibframe__deepcoadd`
2.  `frame`
3.  `frame_forcelist__deepcoadd__iselect`
4.  `frame_forcephoto__deepcoadd__iselect`
5.  `frame_matchlist`
6.  `frame_matchphoto`
7.  `frame_sourcelist`
8.  `frame_sourcephoto`
9.  `mosaic__deepcoadd`
10. `mosaic_forcelist__deepcoadd__iselect`
11. `mosaic_forcephoto__deepcoadd__iselect`
12. `mosaic_sourcelist__deepcoadd`
13. `mosaic_sourcephoto__deepcoadd`

#### Also, these catalog may also be helpful 

1.  `exposure`
2.  `warped__deepcoadd`
3.  `wcs`

#### Image Metadata: 

##### Frame & Frame_Mng

* FRAME: contains information for CORR (reduced CCD images); Useful information 
         includes: 
         1. Frame_id, Frame_Num....All kinds of ID and Index information 
         2. RA, Dec for the frame cetner, four corners, WCS information 
         3. Information about the observation (Date, Time, Airmass...)
         4. Information about the CCD (GAIN, Zeropt ...)
  - `CTYPE1, CTYPE2, CUNIT1, CUNIT2, CRPIX1, CRPIX2, CRVAL1, CRVAL2` and 
    `CD1_1, CD1_2, CD2_1, CD2_2`: Standard WCS information 
  - `LLC\ULC\URC\LRC[RA/DECL]: RA, Dec for four corners
  - `PA, INSROT`: Angle of the CCD during observation
  - `DATE_OBS, MJD, TAIOBS, UT, HST, LST`: Date and time
  - `Azimuth, Elevation, Airmass`: Configuration during observation
  - `FILTER, FILTER01`: Filter information
  - `EXPTIME`: Exposure time 
  - `DATA_TYPE`: OBJECT or STANDARD_STAR
  - `GAIN1/2/3/4, OSLEVEL1/2/3/4, OSSIGMA1/2/3/4`: GAIN and OverScan information
  - `SKYLEVEL, SIGMA_SKY`: median and sigma of sky value
  - `SEEING, ELLIPT, ELL_PA`: Seeing and PSF information
  - `ZEROPT, MAGZERO_RMS, MAGZERO_NOBJ`: Magnitude zeropoint information
  - `NOBJ_BRIGHT, NOBJ_MATCHED`: Number of bright objects, and objects for 
    astrometric match

* FRAME_MNG: Management table; file locations
  - `IMG_LOC, CAT_LOC, MAT_LOC, PSF_LOC`

##### MOSAIC, CALIBFRAME, WARPED
  - MOSAIC --> Information for the CALEXP coadded images 
    * `TRACT, PATCH, PATCH_NUM, POINTING, RERUN, MOS_RERUN`
    * `NAXIS1, NAXIS2`
    * `RA2000, DECL2000`: Frame center 
    * The other information are very similar to FRAME
  - MOSAIC_MNG --> Management catalog
    * `IMG_LOC, CAT_LOC, MAT_LOC`: Location of files

#### Table for catalog data: 

* At the current stage, the following two tables are the most useful ones: 
  1. **PHOTOOBJ_MOSAIC**: Coadd forced sources 
  2. **PHOTOOBJ_FRAME** : Reduced single frame image 
  - In the future there will be `deepcoadd` and `bestseeing` release for 
    the Mosaic catalog;  So far, only deepcoadd is available
  - For mosaic forced measurements, there are `iselect` and `rselect`
    information (for HSC_Wide right now, only `iselect` is available)

* The full name of a table should be: 
  `schema_name.table_root_name__(mos_rerun)__(cat_rerun)`
  - e.g. `ssp_s14a0_udeep_20140523a.photoobj_mosaic__deepcoadd__iselect`

* For single frame data, the following catalog are useful: 
  - FRAME: Use frame_ID as index, provide information about the frame 
  - FRAME_SOURCELIST: use frame_ID as index; 
      * `DEBLEND_NCHILD`: Number of children
      * `CENTROID_[NAIVE|SDSS]_[X|Y]`: Center of the object using NAIVE and 
        SDSS algorithm 
      * `SHAPE_SDSS`: Shape measured with SDSS adaptive moment algorithm 
      * `SHAPE_SDSS_CENTROID_[X|Y]`: Centroid 
      * `SHAPE_SDSS_PSF`: Adaptive moments of the PSF model
      * `FLUX_[NAIVE|GAUSSIAN|PSF|SINC|KRON]`: Different photometry
      * `CLASSIFICATION_EXTENDEDNESS`: Probability of being extended object
      * `FLAG_FLAGS_PIXEL_[EDGE/INTERPOLATE_ANY/INTERPOLATED_CENTER/
           SATURATED_ANY/SATURATED_CENTER/CR_ANY/CR_CENTER/BAD/SUSPECT_ANY/
           SUSPECT_CENTER]: Set if certain issue is related to the object 
      * `FLAG_CALIB_PSF_[CANDIDATE|USED]: if used for PSF   
      * `FLAG_CLASSIFICATION_PHOTOMETRIC`: If used for photometric calibration
      * `FLAG_DEBLEND_[DEBLENDED_AS_PSF|TOO_MANY_PEAKS|FAILED|PARENT_TOO_BIG| 
           SKIPPED|RAMPED_TEMPLATE|PATCHED_TEMPLATE|HAS_STRAY_FLUX]: Set if such 
           issue appears in the DEBLEND process
      * `FLAG_FLUX_KRON_FLAGS`: Set if Flux.Kron failed 
      * `FLAG_FLUX_KRON_RADIUS`: Set if Kron radius is bad 
      * `FLAG_FLUX_KRON_SMALLRADIUS`: Set if Kron radius is smaller than PSF
  - FRAME_SOURCEPHOTO: 
      * `FILTER01`: Filter name 
      * `MAG_[GAUSSIAN|NAIVE|PSF|KRON|SINC]` and `_ERR`
      * `MAG_[CMODEL|CMODEL_EXP|CMODEL_DEV]` **NOT AVAILABLE YET**
  - FRAME_FORCELIST, FRAME_FORCEPHOTO will become available in the future

* For mosaic images, `MOSAIC_SOURCELIST__DEEPCOADD` and
  `MOSAIC_SOURCEPHOTO__DEEPCOADD` are the most useful ones
  - `TRACT`, `PATCH`, `POINTING`, and `ID` should be used to 
    join these two table (In the future, only `ID` would be enough)

  - MOSAIC_SOURCELIST: 
    * `RERUN, MOS_RERUN, FILTER01`
    * All the flux and flags from FRAME_SOURCELIST
    * `MULTISHAPELET_PSF_[INNER|OUTER|ELLIPSE|CHISQ|INTEGRAL]`: Shapelet PSF 
    * `CMODEL_[FLUX|FLUX_ERR|CENTER|EXP_FLUX|EXP_FLUX_ERR|EXP_ELLIPSE|
                 DEV_FLUX|DEV_FLUX_ERR|DEV_ELLIPSE]`: Cmodel Flux measurements
    * There are a bunch of FLAGS for CMODEL photometry, which are not available 
      right now.  Probably the CMODEL photometry is not very good
  - MOSAIC_SOURCEPHOTO: 
    * All the magnitude from FRAME_SOURCEPHOTO 
    * `MAG_CMODEL[_ERR|_EXP|_EXP_ERR|_DEV|_DEV_ERR]: Cmodel magnitude
      **NOT AVAILABLE YET**
  - MOSAIC_FORCEFLAG_[G|R|I|Z|Y} are not available yet

#### Useful functions 

##### Aggregate 

* `qmedian( expression )`
* `quantile( expression, fraction )`
  - e.g. `SELECT quantile( seeing, 0.3)`
  - e.g. `SELECT quantile( seeing, array[ 0.3, 0.5, 0.7 ] )`
* `weighted_mean( values, weight_values )`
* `mean, variance, stddev, skewness, kurtosis`
  - e.g. `SELECT mean( mag_kron, '<', 99.99 )`

#### Spatial Searches 

* `f_getobj_circle( ra, dec, radius, table_name )`
  - e.g., `SELECT * FROM f_getobj_circle(150.403189, 1.485288, 2.0,
    'ssp_s14a0_udeep_20140523a.frame_forcelist__deepcoadd__iselect')`
* `f_getobj_rectangle(ra, dec, delta_ra, delta_dec, table_name)`

#### Utils 

* `frameid2visitccd` and `visitccd2frameid`
  - e.g. `frameid2visitccd('HSCA00000301')`
  - e.g. `visitccd2frameid(2, 27)`
* `hms2deg` and `dms2deg`; `deg2hms` and `deg2dms`

#### WCS related 

* `sky2pix` and `pix2sky` 
  - `sky2pix( ra, dec, schema, tract, frame-id )`
    * e.g., `sky2pix(150.5, 1.5,'ssp_s14a0_udeep_20140523a', 0, 'HSCA00188753')`
  - `pix2sky( x, y, schema, tract, frame-id )`
    * e.g., `pix2sky(1750.325,359.630,'ssp_s14a0_udeep_20140523a', 0, 'HSCA00188753')`

* `shape_sky2pix` and `shape_pix2sky`
  - e.g., `shape_pix2sky(shape_sdss, centroid_sdss_x, centroid_sdss_y,
    'ssp_s14a0_wide_20140523a', tract, frame_id)`

--------------------------------------------------------------------------------

### Extract Useful Information from the current Databse 

* ssp_s14a0_wide_20140523a:  HSC-I and HSC-Y 
* ssp_s14a0_udeep_20140523a: HSC-G, R, I, Z, Y (but g-band is bad)

#### Frame information

1. hsc_wide_20140523a_frame_mng.csv             10399 
2. hsc_wide_i_20140523a_frameinfo.csv            4887 
3. hsc_wide_y_20140523a_frameinfo.csv            5512 

4. hsc_udeep_20140523a_frame_mng.csv            11905  
5. hsc_udeep_g_20140523a_frameinfo.csv            311
6. hsc_udeep_r_20140523a_frameinfo.csv           1248 
7. hsc_udeep_i_20140523a_frameinfo.csv           1144
8. hsc_udeep_z_20140523a_frameinfo.csv           1664
9. hsc_udeep_y_20140523a_frameinfo.csv           7538 

#### Mosaic information 

1. hsc_wide_i_20140523a_mosaic.csv                839
2. hsc_wide_y_20140523a_mosaic.csv                972

3. hsc_udeep_g_20140523a_mosaic.csv                75
4. hsc_udeep_r_20140523a_mosaic.csv                83
5. hsc_udeep_i_20140523a_mosaic.csv                84
6. hsc_udeep_z_20140523a_mosaic.csv                84
7. hsc_udeep_y_20140523a_mosaic.csv                85

--------------------------------------------------------------------------------

#### Basic Sample Selection

* DEBLEND_NCHILD and CLASSIFICATION_EXTENDEDNESS are not available for 
  photoobj_mosaic__deepcoadd__iselect

* The photometry is calibrated against SDSS photometry, so it is basically on 
  AB system 
* Aperture photometry: in the config file, there are two groups of parameters 
  for aperture photometry, and the radius definitions are different: 
  1. flux.aperture: 3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0
  2. flux.aperture.elliptical: 1.0, 1.5625, 2.44140625, 3.814697265625,
     5.9604644775390625, 9.313225746154785, 14.551915228366852, 22.737367544323206,
     35.52713678800501, 55.51115123125783 
  - The latter ones should be the ones used for forced photometry

##### WIDE (GAMA) 

* Select from photoobj_mosaic__deepcoadd__iselect & 
  mosaic_forceflag__deepcoadd__iselect
  - FLUX_KRON_FLAGS, FLUX_KRON_FLAGS_RADIUS, FLUX_SINC_FLAGS & 
    CENTROID_SDSS_FLAGS are required to be NOT True
  1. hsc_wide_20140523a_deepcoadd_iselect_21.0.csv    1047431  
  2. hsc_wide_20140523a_deepcoadd_iselect_22.0.csv     917793
  3. hsc_wide_20140523a_deepcoadd_iselect_23.0.csv    1708780
  4. hsc_wide_20140523a_deepcoadd_iselect_23.5.csv    1307332 
  5. hsc_wide_20140523a_deepcoadd_iselect_24.0.csv    1634468
  6. hsc_wide_20140523a_deepcoadd_iselect_24.5.csv    1830765

* Updated ones
  - The earlier ones forgot to set 
    `mosaic_forceflag__deepcoadd__iselect.filter01 = 'HSC-I'` during the search, 
    so the number of objects is incorrect, and over-estimated
  1. hsc_wide_20140523a_deepcoadd_iselect_13.0_21.0.csv    563158
  2. hsc_wide_20140523a_deepcoadd_iselect_21.0_22.0.csv    476372
  3. hsc_wide_20140523a_deepcoadd_iselect_22.0_23.0.csv    907567
  4. hsc_wide_20140523a_deepcoadd_iselect_23.0_24.0.csv   1636039
  5. hsc_wide_20140523a_deepcoadd_iselect_24.0_25.0.csv   2148431
  6. hsc_wide_20140523a_deepcoadd_iselect_25.0_26.0.csv   1528992  


##### UDEEP (COSMOS) 

* Select from photoobj_mosaic__deepcoadd__iselect & 
  mosaic_forceflag__deepcoadd__iselect
  - FLUX_KRON_FLAGS, FLUX_KRON_FLAGS_RADIUS, FLUX_SINC_FLAGS & 
    CENTROID_SDSS_FLAGS are required to be NOT True
  1. hsc_udeep_20140523a_deepcoadd_iselect_23.5.csv    1025301  
  2. hsc_udeep_20140523a_deepcoadd_iselect_24.5.csv     710942  
  3. hsc_udeep_20140523a_deepcoadd_iselect_25.5.csv     584215  
  4. hsc_udeep_20140523a_deepcoadd_iselect_26.5.csv     230620  
  5. hsc_udeep_20140523a_deepcoadd_iselect_27.5.csv      37889  

* Updated ones: 
  1. hsc_udeep_20140523a_deepcoadd_iselect_21.0.csv    41993 
  2. hsc_udeep_20140523a_deepcoadd_iselect_23.0.csv   162423 
  3. hsc_udeep_20140523a_deepcoadd_iselect_26.0.csv   583368 
  4. hsc_udeep_20140523a_deepcoadd_iselect_27.0.csv   612105 
  5. hsc_udeep_20140523a_deepcoadd_iselect_28.0.csv   617016 
  6. hsc_udeep_20140523a_deepcoadd_iselect_28.5_dirty.csv   679171 
  7. hsc_udeep_20140523a_deepcoadd_iselect_29.5_dirty.csv   680964 
  8. hsc_udeep_20140523a_deepcoadd_iselect_29.5_clean.fits  601848 
  9. hsc_udeep_20140523a_deepcoadd_iselect_29.5_cleaner.fits  538982 
     - Include the SHAPE_SDSS_FLAGS

* Join photoobj_mosaic__deepcoadd__iselect with:   
  - mosaic_forcelist__deepcoadd__iselect on `photoobj_mosaic.id =
      forcelist.object_id`
  - mosaic_forcephoto__deepcoadd__iselect on `photoobj_mosaic.id =
      forcephoto.id`
  - mosaic_forceflag__deepcoadd__iselect on `photoobj_mosaic.id =
      forceflag.object_id`
  * **VERY CONFUSING** 
  * hsc_udeep_20140523a_combine_iselect_28.5_dirty.csv 565845 
  * hsc_udeep_20140523a_combine_iselect_28.5_clean.csv 509803 
    - When combined the FORCEPHOTO catalog, the number of objects are much
      smaller. **TODO: Why? What if only select from Forcephoto?**
  * hsc_udeep_20140523a_combine_iselect_29.5_dirty.csv     567429
  * hsc_udeep_20140523a_combine_iselect_29.5_clean.fits    511330
  * hsc_udeep_20140523a_combine_iselect_29.5_cleaner.fits  466446
    - GALAXY field = classification_extendedness 
    - DEBLEND_NCHILD

--------------------------------------------------------------------------------

### Match with COSMOS 

* Use the v1.8 photo-z catalog provided by Alexie.  
* Select all the useful galaxies: `TYPE == 0 && STAR_ACS == 1 && AUTO_FLAG !=-1`
  - `TYPE == 0`: Galaxy 
  - `STAR_ACS == 1`: Galaxy; Based on star-galaxy separation from Alexie's catalog 
  - `AUTO_FLAG != -1`: Defined by Capak+2007; Has useful total magnitude
  - 652848/2017800
  - If further request `i_mag_err >= 0`:  641947

* **WARNING**: The following FLAGS are not actually working!!
  1. flags_pixel_cr_any
  2. flux\_kron\_flags / flux\_sinc\_flags / flux\_kron\_flags\_radius 


* Use hsc_udeep_20140523a_deepcoadd_iselect_26.0.csv (583368) 
  - flags_pixel_edge                8647 
  - flags_pixel_interpolate_any     8324
  - flags_pixel_interpolate_center  4442
  - flags_pixel_saturated_any      16477 
  - flags_pixel_saturated_center    6204 
  - flags_pixel_suspect_any          252
  - flags_pixel_suspect_center       203
  - The flags_pixel_cr_[any/center] flags do not seem to work --> 0
  - `! flags_pixel_edge && ! flags_pixel_interpolated_any && !
      flags_pixel_suspect_any&& ! flags_pixel_saturated_any`
      * Define a CLEAN sample: 565254 (97%) 

* Use hsc_udeep_20140523a_deepcoadd_iselect_27.0.csv (612105) 
  - flags_pixel_edge                8661 
  - flags_pixel_interpolate_any     8338
  - flags_pixel_interpolate_center  4447
  - flags_pixel_saturated_any      16643 
  - flags_pixel_saturated_center    6273 
  - flags_pixel_suspect_any          252
  - flags_pixel_suspect_center       203
  - The flags_pixel_cr_[any/center] flags do not seem to work --> 0
  - `! flags_pixel_edge && ! flags_pixel_interpolated_any && !
      flags_pixel_suspect_any&& ! flags_pixel_saturated_any`
      * Define a CLEAN sample: 593816 (97%) 

* Use hsc_udeep_20140523a_deepcoadd_iselect_28.0.csv (617015) 
  - flags_pixel_edge                8666 
  - flags_pixel_interpolate_any     8343
  - flags_pixel_interpolate_center  4448
  - flags_pixel_saturated_any      16671 
  - flags_pixel_saturated_center    6282 
  - flags_pixel_suspect_any          252
  - flags_pixel_suspect_center       203
  - The flags_pixel_cr_[any/center] flags do not seem to work --> 0
  - `! flags_pixel_edge && ! flags_pixel_interpolated_any && !
      flags_pixel_suspect_any&& ! flags_pixel_saturated_any`
      * Define a CLEAN sample: 598694 (97%) 

* Use hsc_udeep_20140523a_deepcoadd_iselect_28.5_dirty.csv (679170)
  - Can define a clean sample from it: 600137

* Cross-Match using 1.5 arcsec separation 
  - 349351 Pairs (54%); 317773 clean ones
* Cross-Match using 1.0 arcsec separation 
  - 346808 Pairs (53%); 315648 clean ones
* Cross-Match using 0.5 arcsec separation 
  - 343963 Pairs      ; 313556 clean ones 
* Cross-Match using 0.3 arcsec separation 
  - 317298 Pairs      ; 290594 clean ones

#### New XMatch (2014/08/13) 


* For HSC survey: 
  1. ...deepCoadd_iselect_29.5_dirty    680963 
  2. ...deepCoadd_iselect_29.5_cleaner  538982 
  - **For above HSC catalog, there is no Star/Galaxy Separation right now!!**
  3. hsc_udeep_20140523a_combine_iselect_29.5_dirty.fits         567429 
     - GALAXY: `galaxy == 1`:                                    499684
     - GALAXY_NON\_DEBLEND: `galaxy == 1 && deblend_nchild < 2`: 397444
  4. hsc_udeep_20140523a_combine_iselect_29.5_cleaner.fits       466446 
     - GALAXY: `galaxy == 1`:                                    406711
     - GALAXY_NON\_DEBLEND: `galaxy == 1 && deblend_nchild < 2`: 309276

* For COSMOS survey: 
  1. cosmos_galaxy_photoz_v1.8          652848 
     - Select galaxies with useful imag: 
     - `imag_Subaru >= 10.0 && imag_Subaru <= 32.0`  641947
     - Correct the Subaru r, i, and z-band magnitude to AUTO by adding the 
       auto_offset back in
  2. cosmos_ultravista_rselect          773712
     - Select a clean galaxy sample using: 
     - `star == 0 && contamination == 0 && nan_contam < 3 && eip >= 0.0` 551573

* Two cross-match
  1. "Detect": sanity checks for source detections; Using the "Dirty" sample; 
     - Separation = 1.0 arcsec 
  2. "Photom": for the comparison of photometry; Using the "Clean" or "Useful" 
     sample 
     - Separation = 0.5 arcsec

##### Group1 

- Session: **hsc_cosmos_compare_2.fits**

1. hsc_cosmos_xmatch_detect_1: HSC Coadd Dirty v.s. COSMOS Photoz
   - Sep = 1.0 arcsec 
   - Return 347893 objects

2. hsc_cosmos_xmatch_photom_1: HSC Coadd Cleaner v.s. COSMOS Photoz (imag)
   - Sep = 0.5 arcsec 
   - Return 276573 objects 
   
3. hsc_cosmos_xmatch_photom_2: HSC Coadd Cleaner v.s. Ultravista r-select (clean)
   - Sep = 0.5 arcsec 
   - Return 267986 objects

##### Group2 

- Session: **hsc_cosmos_compare_3.fits**

1. hsc_cosmos_xmatch_detect_2: HSC Combine Dirty (galaxy) v.s. COSMOS Photoz 
   - Sep = 1.0 arcsec 
   - Return 274656

2. hsc_cosmos_xmatch_photom\_3: HSC Combine Cleaner (galaxy_non_blend) v.s.
   COSMOS Photoz (imag) 
   - Sep = 0.5 arcsec 
   - Return 203300

3. hsc_cosmos_xmatch_photom\_4: HSC Combine Cleaner (galaxy) v.s. COSMOS Photoz (imag) 
   - Sep = 0.5 arcsec 
   - Return 226595

4. hsc_cosmos_xmatch_photom\_5: HSC Combine Cleaner (galaxy) v.s. Ultravista 
   r-select (clean)
   - Sep = 0.5 arcsec 
   - Return 216407

* F814W from Subaru i and z; According to Capak+2007: 
  - `( zmag_kron - a_z + 0.032 ) + 0.632*( imag_kron - a_i + 0.093 - zmag_kron - a_z
- 0.032) - 0.116*pow( ( imag_kron - a_i + 0.093 - zmag_kron + a_z - 0.032), 2.0
) - 0.001`

* [2014/08/14] The F814W Mag_Auto has been extracted from the Capak+2008 catalog 
  and joined with the cosmos_galaxy_photoz_v1.8 catalog. 
  - cosmos_galaxy_photoz_v1.8_acs.fits 

* [2014/08/14] The total flux for the SERSICFIT and BULGEFIT models in 
  Claire's catelog have been calculated
  - The total flux has been merged into the original table 
  - The zeropoint to convert the flux into F814W magnitude is not very clear 
    right now; probably close to 26.0

* [2014/08/14] New Cross Match to study the detection of objects 
  - cosmos_galaxy_photoz_v1.8_acs [652,663] 
  - hsc...combine_iselect_29.5_dirty_galaxy [499684] 
  - the session is saved as **hsc_cosmos_xmatch_compare_4.fits**
  - 1.0 arcsec max distance: Matched:       274655
  - 1.5 arcsec max distance: HSC_No_COSMOS: 223158
  - 1.5 arcsec max distance: COSMOS_No_HSC: 376137
  - Also define a common area for the last two subsamples
    * xmatch_hsc_no_cosmos_common.fits  [54934]
    * xmatch_cosmos_no_hsc_common.fits [153210]

* [2014/08/15] Use new COSMOS_Photoz_v1.8_ACS and COSMOS_Lackner_Model_Totflux
  tables for Xmatch 
  - Session saved as **hsc_cosmos_compare_5.fits**

1. **hsc_cosmos_xmatch_photom\_6**: 
   -HSC Combine Cleaner (galaxy_non_blend) v.s. COSMOS Photoz ACS (imag) 
   - Sep = 0.5 arcsec 
   - Return 203607 

2. **hsc_cosmos_xmatch_photom\_7**: 
   -HSC Combine Cleaner (dirty) v.s. COSMOS Photoz ACS (imag) 
   - Sep = 0.5 arcsec 
   - Return 272543 

3. **hsc_cosmos_xmatch_photom\_8**: 
   -HSC Combine Cleaner (cleaner) v.s. COSMOS Lackner's Modelfit (totflux) 
   - Sep = 0.5 arcsec 
   - Return 37280 

---- 

## Weekly Report 

### Sample Selection for HSC UDeep Early-Data 

* The HSC sample is made from the i-band selected objects from the 
  `deepCoadd` mosaic images.
  - 1. All valid detections with `imag_kron < 29.5 mag` (i-band Kron magnitude)
       are extracted from the `Photoobj_Mosaic` table of the IPMU database; 
       There are about **680963** detections. This sample is labeled as **DIRTY**
  - 2. These detections are further joined with the `Mosaic_Forcelist`, 
       `Mosaic_Forcephoto`, and `Mosaic_Forceflag` tables to get other necessary 
       information.  And the `classification_extendedness == 1` criteria is 
       applied to exclude stellar objects. This reduced the sample size to 
       **499684**. This sample is labeled as **CLEAN**
  - 3. Available quality-control flags are used to select a sample of galaxies  
       with (relative) reliable photometry. These galaxies have well-defined 
       center, are not too close to the edge of the image, are not affected by 
       interpolated, saturated pixels. Furthermore, they have 
       `deblend_nchild == 1`.  This further reduce the sample size to
       **309276**.  The sample is labeled as **DECENT**. 
  - 4. Galactic extinction has been corrected using the value provided in the 
       `Photoobj_Mosaic` catalog. 

### Sample Selection from Existed Photometric Data in COSMOS Field  

* The reference sample is selected from the COSMOS multi-band photometry 
  catalog (`Photoz_v1.8`).  This catalog includes photometric measurements from 
  the Subaru Suprime-Camera observations of the COSMOS field (The g+, r+, i+, 
  and z+ data are relevant).  Star-Galaxy separation parameter, along with other 
  quality control flags are applied to select a sample of galaxies with useful 
  Subaru i+ band photometry (641947).  Then, the r+, i+, and z+ magnitude are 
  corrected back to AUTO_MAG using the **auto_offset** parameter. 
    
### Number of Detected Object in i-band 
       
* According to Appendix A, the current HSC survey UDEEP footprint is slightly 
  larger than the original HST/ACS COSMOS field.  However, it is easy to notice 
  that the number of detections using the i-band `deepCoadd` images is 
  much lower than the galaxies found on HST/ACS F814W images (Also see Fig.1).
  It is likely due to the limited depth of the currently available i-band
  images.

Fig.1: The i-band magnitude distribution for galaxies from the Suprime-Camera 
       observations (Red line), and for galaxies detected in the HSC survey 
       early data (Filled region with different color). 

* For Suprime-Camera data, the 5-$\sigma$ depth (95% complete) in a 2.1 arcsec
  aperture is around 25.9 to 26.2 magnitude; Such value for the current HSC 
  survey i-band data is not accurately estimated.
