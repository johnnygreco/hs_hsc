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
#
#cic.coaddImageCutout(wideRoot, largeRa, largeDec, largeSize, saveMsk=True,
#                     saveSrc=True, filt='HSC-I', prefix='hsc_wide2_cutout',
#                     circleMatch=True)
#
## Exampe galaxies in the UDEEP field
## From COSMOS campact galaxies: ID=6026402
#udeepRa  = 149.75626
#udeepDec = 2.794748
#
#################################################################################
## Get the common patches shape for surveys at different bands
#import coaddPatchShape as cps
#
#wkbList = cps.batchPatchNoData('/home/song/work/early/cosmos_341_corners/',
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
