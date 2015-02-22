#!/usr/bin/env python

import argparse
import numpy                  as np
import matplotlib.pyplot      as pyplot
import lsst.daf.persistence   as dafPersist
import lsst.afw.image         as afwImage
from astropy.io import fits

def getPixelScale(wcs):

    pixelScale = wcs.pixelScale().asArcseconds()
    return pixelScale

def getImgArr(mImg):

    # Get the numpy ndarray for the image
    imgArr = mImg.getImage().getArray()
    # Get the numpy ndarray for the mask image
    mskArr = mImg.getMask().getArray()
    # Get the numpy ndarray for the variance
    varArr = mImg.getVariance().getArray()
    # Turns the variance array into sigma array
    sigArr = np.sqrt(varArr)

    return imgArr, mskArr, sigArr

def getBadArr(calExp):

    # Get the mask image
    mskImgBad = calExp.getMaskedImage().getMask()
    # Remove the "DETECTED" Mask from the mask image
    mskImgBad.clearMaskPlane(5)
    # Remove the "EDGE" Mask from the mask image XXX TODO: Check if this is good
    mskImgBad.clearMaskPlane(4)

    return mskImgBad.getArray()

def getDetArr(calExp):

    # Get the mask image
    mskImgDet = calExp.getMaskedImage().getMask()
    # Remove all the mask planes except for the "DETECTED" one
    # TODO: This is very stupid way of doing this !
    mskImgDet.clearMaskPlane(0)   # BAD
    mskImgDet.clearMaskPlane(1)   # SAT
    mskImgDet.clearMaskPlane(2)   # INTRP
    mskImgDet.clearMaskPlane(3)   # CR
    mskImgDet.clearMaskPlane(4)   # EDGE
    mskImgDet.clearMaskPlane(6)   # DETECTED_NEGATIVE
    mskImgDet.clearMaskPlane(7)   # SUSPECT
    mskImgDet.clearMaskPlane(8)   # NO_DATA
    mskImgDet.clearMaskPlane(9)   # CROSS_TALK
    mskImgDet.clearMaskPlane(10)  # UNMASKEDNAN

    return mskImgDet.getArray()

def getBkgArr(bImg):

    bkgArr = bImg.getImage().getArray()
    return bkgArr

def savePsfArr(calExp, psfFile):

    psfImg = calExp.getPsf().computeImage()
    psfImg.writeFits(psfFile)
    return psfImg.getArray()

def getZeroPoint(butler, dataId):

    calMeta = butler.get('deepCoadd_md', dataId, immediate=True)
    zeropoint = (2.5 * np.log10(calMeta.get("FLUXMAG0")))
    return zeropoint

def coaddImagePrepare(rootDir, tract, patch, filt, prefix):

    # Make a butler and specify the dataId
    butler = dafPersist.Butler(rootDir)
    dataId = {'tract': tract, 'patch': patch, 'filter': filt}

    # Get the prefix of the output files
    prefix = prefix + '-' + str(tract) + '-' + patch + '-' + filt

    # Get the exposure from the butler
    calExp = butler.get('deepCoadd', dataId, immediate=True)
    # Get the WCS information
    wcs = calExp.getWcs()
    # Get the pixel scale of the image
    pixelScale = getPixelScale(wcs)

    # Get the flux zeropoint
    zeropoint = getZeroPoint(butler, dataId)

    # Get the maskedImageObject
    mImg = calExp.getMaskedImage()
    imgArr, mskArr, sigArr = getImgArr(mImg)

    # Get the "Bad" mask array
    badArr = getBadArr(calExp)
    # Get the "Detected" mask array
    detArr = getDetArr(calExp)

    # Get the numpy ndarray for the background
    bImg = butler.get('deepCoadd_calexpBackground', dataId, immediate=True)
    bkgArr = getBkgArr(bImg)

    # Add the background back on to the image
    #oriArr = (imgArr + bkgArr)

    # Get the PSF array
    psfFile = prefix + '_psf.fits'
    psfArr = savePsfArr(calExp, psfFile)

    # Get the header of the images
    fitsName = butler.get('deepCoadd_filename', dataId)
    print fitsName
    hduList = fits.open(fitsName[0])
    imgHeader = hduList[1].header
    mskHeader = hduList[2].header
    varHeader = hduList[3].header

    # Update the header with more information
    imgHeader.set('PIXSCALE', pixelScale)
    imgHeader.set('ZP_PHOT',  zeropoint)

    # Write the Image array
    imgFile = prefix + '_img.fits'
    fits.writeto(imgFile, imgArr, imgHeader)
    # Write the Mask array
    mskFile = prefix + '_msk.fits'
    fits.writeto(mskFile, mskArr, mskHeader)
    badFile = prefix + '_bad.fits'
    fits.writeto(badFile, badArr, mskHeader)
    detFile = prefix + '_det.fits'
    fits.writeto(detFile, detArr, mskHeader)
    # Write the Sigma array
    sigFile = prefix + '_sig.fits'
    fits.writeto(sigFile, sigArr, varHeader)
    # Write the Image + Background array
    bkgFile = prefix + '_bkg.fits'
    fits.writeto(bkgFile, bkgArr, imgHeader)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("tract", type=int, help="Tract to show")
    parser.add_argument("patch",  help="Patch to show")
    parser.add_argument('-f', '--filter', dest='filt', help="Filter",
                       default='HSC-I')
    parser.add_argument('-p', '--prefix', dest='outfile',
                        help='Prefix of the output file', default='hsc_coadd')
    args = parser.parse_args()

    coaddImagePrepare(args.root, args.tract, args.patch, args.filt, args.outfile)
