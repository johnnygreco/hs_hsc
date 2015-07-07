#!/usr/bin/env python
# encoding: utf-8

# Examples:

################################################################################
## Get the cutout image
#import coaddImageCutout as cic
#
## Example galaxies in the WIDE field
#wideRoot = '/lustre/Subaru/SSP/rerun/yasuda/SSP3.4.1_20141224/'
## From redmapper ID=16474; z=0.38
#wideRa  = 133.02153
#wideDec = 1.3656108
#wideSize = 200
#
#cic.coaddImageCutout(wideRoot, wideRa, wideDec, wideSize, saveMsk=True,
#                     saveSrc=True, filt='HSC-I', prefix='hsc_wide1_cutout',
#                     circleMatch=True)
#
## Example of very large galaxy cutout
#largeRa   = 133.65253
#largeDec  = 0.64257544
#largeSize = 1000
##
#cic.coaddImageCutFull(wideRoot, largeRa, largeDec, largeSize, savePsf=True,
                     #saveSrc=True, filt='HSC-I', prefix='hsc_wide2_cutout',
                     #verbose=True, visual=True)
#
## Exampe galaxies in the UDEEP field
## From COSMOS campact galaxies: ID=6026402
#udeepRa  = 149.75626
#udeepDec = 2.794748
#
#################################################################################
#import coaddBatchCutout as cbc

#wideRoot = '/lustre/Subaru/SSP/rerun/yasuda/SSP3.4.1_20141224/'
#inCat = '/home/song/work/hs_hsc/py/data/hsc_redmapper_xmm.fits'
#cbc.coaddBatchCutout(wideRoot, inCat, size=400, filter='HSC-I', prefix='redmapper_xmm')

#################################################################################
## Get the common patches shape for surveys at different bands
#import coaddPatchShape as cps
#
#wkbList = cps.batchPatchShape('/home/song/work/early/cosmos_341_corners/',
#                               'ssp341_cosmos*corners.fits')
#
#################################################################################
## Get the corner (Ra,Dec) for a list of images
#import coaddImgCornerRaDec as cicr
#
#multi = cicr.coaddAllBandsCorners('/lustre/Subaru/SSP/rerun/yasuda/ \
#                                  SSP3.4.1_20141224/',
#                                  prefix='ssp341_cosmos')
#

#################################################################################
## Get the NO_DATA mask for a list of images
#import coaddPatchNoData as cpnd
#
## Example of single image
#rootDir = '/lustre/Subaru/SSP/rerun/yasuda/SSP3.4.1_cosmos_setWeight'
#tract   = 0
#patch   = '0,6'
#filter  = 'HSC-I'
#
#noDataConvex = cpnd.coaddPatchNoData(rootDir, tract, patch, filter,
#                                     prefix='ssp341_cosmos', savePNG=True)
#
## Example of batch mode
#rootDir = '/lustre/Subaru/SSP/rerun/yasuda/SSP3.4.1_cosmos_setWeight'
#prefix  = 'ssp341_cosmos'
#filter  = 'HSC-I'
#
#noDataUse, noDataComb = cpnd.batchPatchNoData(rootDir, filter=filter,
#                                              prefix=prefix)

#################################################################################
# Match external catalogs

#""" For Thinkpad """
##inCat = "/home/hs/Dropbox/work/project/hs_hsc/py/data/test_extra_cat.fits"
##acpMask = "/home/hs/Dropbox/work/project/hs_hsc/py/ipynb/ssp341_cosmos_HSC-I_corners.wkb"
##rejMask = "/home/hs/Dropbox/work/project/hs_hsc/py/ipynb/ssp341_cosmos_0_HSC-I_nodata_all.wkb"

#""" For Master """
#inCat = "/home/song/work/hs_hsc/py/data/test_extra_cat.fits"
#acpMask = "/home/song/work/hs_hsc/py/ipynb/ssp341_cosmos_HSC-I_corners.wkb"
#rejMask = "/home/song/work/hs_hsc/py/ipynb/ssp341_cosmos_0_HSC-I_nodata_all.wkb"

#import coaddMaskMatch as cmm
#cmm.coaddMaskMatch(inCat, acpMask, raField='RA', decField='DEC',
                   #showMatch=True, infoField1=None, infoText1=None,
                   #infoField2=None, infoText2=None, outCat=None,
                   #rejMask=rejMask)

#################################################################################

# Test Cutout Prepare

import coaddCutoutPrepare as coaddPre

#dataDir = '/Users/songhuang/Downloads/cutout/redmapper/example'
dataDir = '/home/hs/Downloads/redmapper_example'
prefix = 'hsc_wide2_cutout_HSC-I_full'

coaddPre.coaddCutoutPrepare(prefix, root=dataDir, srcCat=None, verbose=True,
                           bSizeH=10, bSizeC=400, thrH=3.5, thrC=1.5, mask=2,
                           growC=5.0, growW=2.0, growH=1.5, kernel=6, central=1,
                           galX=None, galY=None, galR1=None, galR2=None, galR3=None,
                           galQ=None, galPA=None, visual=True, suffix='',
                           combBad=True)

