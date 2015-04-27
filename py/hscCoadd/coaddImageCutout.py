#!/usr/bin/env python
# encoding: utf-8

import os
import copy
import argparse
import warnings
import numpy as np

# HSC Pipeline
import lsst.daf.persistence   as dafPersist
import lsst.afw.coord         as afwCoord
import lsst.afw.image         as afwImage
import lsst.afw.geom          as afwGeom
import lsst.afw.table         as afwTable

# Astropy
from astropy import wcs
from astropy.io import fits

# Personal
import coaddColourImage as cdColor


def getCoaddPsfImage(calExp, coord):

    # Get the WCS information
    wcs = calExp.getWcs()
    # The X,Y coordinate of the image center
    coordXY = wcs.skyToPixel(coord)
    # Get the PSF object for the exposure
    psf = calExp.getPsf()
    try:
        psfImg = psf.computeImage(coordXY)
        return psfImg
    except Exception:
        warnings.warn("### Can not compute PSF Image !!!")
        return None


def getCoaddMskPlane(calExp, bitmask):

    # Get the mask image
    mskImg = calExp.getMaskedImage().getMask()
    newMsk = copy.deepcopy(mskImg)
    # Extract specific plane from it
    newMsk &= newMsk.getPlaneBitMask(bitmask)

    return newMsk

def getCoaddBadMsk(calExp):

    # Get the mask image
    mskImg = calExp.getMaskedImage().getMask()

    badMsk = copy.deepcopy(mskImg)
    # Clear the "EDGE" plane
    badMsk.clearMaskPlane(4)
    # Clear the "DETECTED" plane
    badMsk.clearMaskPlane(5)
    # Clear the "DETECTED_NEGATIVE" plane
    badMsk.clearMaskPlane(6)
    # Clear the "CROSSTALK" plane XXX TODO: Check later to see if it is still
    # appropriate
    badMsk.clearMaskPlane(9)

    return badMsk

def getCircleRaDec(ra, dec, size):

    # Get a small but representative set of (RA, DEC) that describe a circle
    # region around the central input coordinate

    # Convert the size from pixel unit to degress
    sizeDegree = (size * 0.168) / 3600.0
    # representative set of polar angles
    angles = np.array([0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0])
    phi = np.array(angles * np.pi / 180.0)

    raList  = ra  + sizeDegree * np.cos(phi)
    decList = dec + sizeDegree * np.sin(phi)

    return raList, decList


def coaddImageCutout(root, ra, dec, size, saveMsk=True, saveSrc=True,
                     filt='HSC-I', prefix='hsc_coadd_cutout',
                     circleMatch=True, verbose=True, extraField1=None, extraValue1=None):

    # Get the SkyMap of the database
    butler = dafPersist.Butler(root)
    skyMap = butler.get("deepCoadd_skyMap", immediate=True)

    # Get the expected cutout size
    dimExpect  = (2 * size + 1)
    cenExpect  = (dimExpect/2.0, dimExpect/2.0)
    sizeExpect = dimExpect ** 2
    # Cutout size in unit of degree
    sizeDeg = size * 0.168 / 3600.0

    # Verbose
    if verbose:
        print '####################################################################'
        print " Input Ra, Dec: %10.5f, %10.5f" % (ra, dec)
        print " Cutout size is expected to be %d x %d" % (dimExpect, dimExpect)

    ############################################################################
    # First, search for the central (Ra, Dec)
    # Define the Ra, Dec pair
    point = afwGeom.Point2D(ra, dec)
    coord = afwCoord.IcrsCoord(point)

    # Search for overlapped tract, patch pairs
    matches = skyMap.findClosestTractPatchList([coord])
    # Number of matched tracts
    nTract = len(matches)
    # Number of matched (patches)
    nPatch = 0
    for tt in range(nTract):
        nPatch += len(matches[tt][1])
    if verbose:
        print "### Find %d possible matches !" % nPatch

    matchCen = []
    for tract, patch in matches:

        # Get the (tract, patch) ID
        tractId = tract.getId()
        patchId = "%d,%d" % patch[0].getIndex()
        if verbose:
            print "### Choose (Tract, Patch) for center: %d, %s !" % (tractId, patchId)
        matchCen.append((tractId, patchId))

        # Get the coadd images
        # Try to load the coadd Exposure; the skymap covers larger area than the
        # available data, which will cause Butler to fail sometime
        try:
            coadd = butler.get("deepCoadd", tract=tractId,
                               patch=patchId, filter=filt, immediate=True)

        except Exception, errMsg:

            print "#############################################"
            print " The desired coordinate is not available !!! "
            print "#############################################"
            print errMsg

            """ TODO """
            coaddFound = False
            noData = True
            partialCut = True
            continue

        else:

            """
            It's still possible that the matched location actually has no useful data
            (e.g. have been interpolated, or in the NO_DATA part of the patch)
            For this situation, no psf image can be generated !!
            """
            coaddFound = True
            # Get the Coadded PSF image
            psfImg = getCoaddPsfImage(coadd, coord)
            if psfImg is None:
                noData = True
                partialCut = True
                continue
            else:
                noData = False

            # Get the WCS information
            wcs = coadd.getWcs()

            # Convert the central coordinate from Ra,Dec to pixel unit
            pixel = wcs.skyToPixel(coord)
            pixel = afwGeom.Point2I(pixel)

            # Define the bounding box for the central pixel
            bbox = afwGeom.Box2I(pixel, pixel)

            # Grow the bounding box to the desired size
            bbox.grow(int(size))

            # Compare to the coadd image, and clip
            bbox.clip(coadd.getBBox(afwImage.PARENT))

            if bbox.isEmpty():
                noData = True
                partialCut = True
                continue
            else:
                if bbox.getArea() < sizeExpect:
                    partialCut = True
                    if verbose:
                        print "### Cut out image dimension " + \
                                "is : %d x %d " % (bbox.getWidth(), bbox.getHeight())
                else:
                    partialCut = False

            # Make a new ExposureF object for the cutout region
            subImage = afwImage.ExposureF(coadd, bbox, afwImage.PARENT)
            # Get the WCS
            subWcs = subImage.getWcs()
            # Get the central pixel coordinates on new subImage WCS
            newX, newY = subWcs.skyToPixel(raDec)
            # Get the new origin pixel
            newOriX, newOriY = subImage.getImage().getXY0()
            newX = newX - newOriX
            newY = newY - newOriY

            # Get the header of the new subimage
            subHead = subImage.getMetadata()
            subHead.set('RA_CUT', ra)
            subHead.set('DEC_CUT', dec)
            subHead.set('NEW_X', newX)
            subHead.set('NEW_Y', newY)
            if partialCut:
                subHead.set('PARTIAL', 1)
            else:
                subHead.set('PARTIAL', 0)
            if (extraField1 is not None) and (extraValue1 is not None):
                subHead.set(extraField1, extraValue1)

            # To see if data are available for all the cut-out region
            if partialCut:
                warnings.warn("### Only part of the cutout region is available" + \
                        " : %d, %s" % (tractId, patchId))
                outPre = prefix + '_' + str(tractId) + '_' + patchId + '_' + \
                        filt + '_cent'
            else:
                outPre = prefix + '_' + str(tractId) + '_' + patchId + '_' + \
                        filt + '_full'

            # Define the output file name
            outImg = outPre + '.fits'
            outPsf = outPre + '_psf.fits'

            # Save the cutout image to a new FITS file
            subImage.writeFits(outImg)

            # Save the PSF image
            psfImg.writeFits(outPsf)

            if saveMsk is True:
                # Get different mask planes
                mskDetec = getCoaddMskPlane(subImage, 'DETECTED')
                mskIntrp = getCoaddMskPlane(subImage, 'INTRP')
                mskSatur = getCoaddMskPlane(subImage, 'SAT')
                mskDetec.writeFits(outPre + '_detec.fits')
                mskIntrp.writeFits(outPre + '_intrp.fits')
                mskSatur.writeFits(outPre + '_satur.fits')

                # Get the "Bad" mask plane
                mskBad = getCoaddBadMsk(subImage)
                mskBad.writeFits(outPre + '_bad.fits')


            if saveSrc is True:

                # Get the source catalog
                """ Sometimes the forced photometry catalog might not be available """
                try:
                    srcCat = butler.get('deepCoadd_forced_src', tract=tractId,
                                        patch=patchId, filter=filt, immediate=True,
                                        flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)
                    # Get the pixel coordinates for all objects
                    srcRa  = np.array(map(lambda x: x.get('coord').getRa().asDegrees(),
                                          srcCat))
                    srcDec = np.array(map(lambda x: x.get('coord').getDec().asDegrees(),
                                          srcCat))
                    # Simple Box match
                    indMatch = ((srcRa > (ra - sizeDeg)) & (srcRa < (ra + sizeDeg)) &
                                (srcDec > (dec - sizeDeg)) & (srcDec < (dec + sizeDeg)))
                    # Extract the matched subset
                    srcMatch = srcCat.subset(indMatch)
                    # Save the src catalog to a FITS file
                    outSrc = outPre + '_src.fits'
                    srcMatch.writeFits(outSrc)
                except:
                    print "### Tract: %d  Patch: %s" % (tractId, patchId)
                    warnings.warn("### Can not find the forced photometry catalog !")
                    if not os.path.isfile('no_src.lis'):
                        noSrc = open('no_src.lis', 'w')
                        noSrc.write("%d  %s \n" % (tractId, patchId))
                        noSrc.close()
                    else:
                        noSrc = open('no_src.lis', 'a+')
                        noSrc.write("%d  %s \n" % (tractId, patchId))
                        noSrc.close()



    # If only part of the desired cutout region is covered, and the circleMatch
    # Flag is set, find all the Patches that overlap with a circle region around
    # the input Ra, Dec
    if partialCut and circleMatch:

        if verbose:
            print "####### Search for other overlapped patches #######"

        # Return the list of RA, DEC that described a circle region around the
        # input (RA, DEC).  The radius is the input size in unit of arcsec
        raList, decList = getCircleRaDec(ra, dec, (size * 0.8))
        points = map(lambda x, y: afwGeom.Point2D(x, y), raList, decList)
        coords = map(lambda x: afwCoord.IcrsCoord(x), points)

        # Search for overlapped tract, patch pairs
        matches = skyMap.findTractPatchList(coords)
        # Number of matched tracts
        nTract = len(matches)
        # Number of matched (patches)
        nPatch = 0
        for tt in range(nTract):
            nPatch += len(matches[tt][1])
        if verbose:
            print "### Find %d possible overlap patches" % nPatch

        for tract, patch in matches:

            # Example all possible matches
            tractId = tract.getId()
            for ii in range(len(patch)):
                patchId = "%d,%d" % patch[ii].getIndex()

                # Skip the image that has been used
                if (tractId, patchId) in matchCen:
                    if verbose:
                        print "### %d - %s is the patch used for CENT cutout" % (tractId, patchId)
                    continue

                # Get the coadd images
                # Try to load the coadd Exposure; the skymap covers larger area than the
                # available data, which will cause Butler to fail sometime
                try:
                    coadd = butler.get("deepCoadd", tract=tractId,
                                       patch=patchId, filter=filt, immediate=True)

                except Exception, errMsg:

                    print "#############################################"
                    print "              PARTIAL OVERLAP                "
                    print " The desired coordinate is not available !!! "
                    print "#############################################"
                    print errMsg

                else:

                    # Get the WCS information
                    wcs = coadd.getWcs()

                    # Convert the central coordinate from Ra,Dec to pixel unit
                    pixel = wcs.skyToPixel(coord)
                    pixel = afwGeom.Point2I(pixel)

                    # Define the bounding box for the central pixel
                    bbox = afwGeom.Box2I(pixel, pixel)

                    # Grow the bounding box to the desired size
                    bbox.grow(int(size))

                    # Compare to the coadd image, and clip
                    bbox.clip(coadd.getBBox(afwImage.PARENT))

                    if bbox.isEmpty():
                        continue
                    elif bbox.getArea() < int(sizeExpect * 0.1):
                        # Ignore small overlapped image
                        if verbose:
                            print "### %d - %s has very small overlapped region" % (tractId,
                                                                                    patchId)
                        continue
                    else:
                        if verbose:
                            print "### Find one useful overlap: %d, %s" % (tractId,
                                                                           patchId)

                    # Make a new ExposureF object for the cutout region
                    subImage = afwImage.ExposureF(coadd, bbox, afwImage.PARENT)
                    # Get the WCS
                    subWcs = subImage.getWcs()
                    # Get the central pixel coordinates on new subImage WCS
                    newX, newY = subWcs.skyToPixel(raDec)
                    # Get the new origin pixel
                    newOriX, newOriY = subImage.getImage().getXY0()
                    newX = newX - newOriX
                    newY = newY - newOriY

                    # Get the header of the new subimage
                    subHead = subImage.getMetadata()
                    subHead.set('RA_CUT', ra)
                    subHead.set('DEC_CUT', dec)
                    subHead.set('NEW_X', newX)
                    subHead.set('NEW_Y', newY)
                    if (extraField1 is not None) and (extraValue1 is not None):
                        subHead.set(extraField1, extraValue1)

                    # To see if data are available for all the cut-out region
                    outPre = prefix + '_' + str(tractId) + '_' + patchId + '_' + \
                             filt + '_part'
                    # Define the output file name
                    outImg = outPre + '.fits'
                    # Save the cutout image to a new FITS file
                    subImage.writeFits(outImg)

                    if saveMsk is True:
                        # Get different mask planes
                        mskDetec = getCoaddMskPlane(subImage, 'DETECTED')
                        mskIntrp = getCoaddMskPlane(subImage, 'INTRP')
                        mskSatur = getCoaddMskPlane(subImage, 'SAT')
                        mskDetec.writeFits(outPre + '_detec.fits')
                        mskIntrp.writeFits(outPre + '_intrp.fits')
                        mskSatur.writeFits(outPre + '_satur.fits')

                        # Get the "Bad" mask plane
                        mskBad = getCoaddBadMsk(subImage)
                        mskBad.writeFits(outPre + '_bad.fits')

    return coaddFound, noData, partialCut


def coaddImageCutFull(root, ra, dec, size, saveMsk=True, saveSrc=True,
                      filt='HSC-I', prefix='hsc_coadd_cutout', verbose=True,
                      extraField1=None, extraValue1=None, skyMap=None):

    # Get the SkyMap of the database
    if skyMap is None:
        try:
            butler = dafPersist.Butler(root)
            skyMap = butler.get("deepCoadd_skyMap", immediate=True)
        except Exception:
            print '### Can not load the correct SkyMap!'

    # Check the filter
    if not cdColor.isHscFilter(filt, short=False):
        raise Exception("### Wrong Filter for HSC Data!")

    # (Ra, Dec) Pair for the center
    raDec = afwCoord.Coord(ra*afwGeom.degrees, dec*afwGeom.degrees)
    # [Ra, Dec] list
    raList, decList = cdColor.getCircleRaDec(ra, dec, size)
    points = map(lambda x, y: afwGeom.Point2D(x, y), raList, decList)
    raDecList = map(lambda x: afwCoord.IcrsCoord(x), points)

    # Expected size and center position
    dimExpect = (2 * size +1)
    cenExpect = (dimExpect/2.0, dimExpect/2.0)
    sizeExpect = dimExpect ** 2
    # Get the half size of the image in degree
    sizeDegree = size * 0.168 / 3600.0

    # Create empty arrays
    imgEmpty = np.zeros((dimExpect, dimExpect), dtype="float")
    mskEmpty = np.zeros((dimExpect, dimExpect), dtype="uint8")
    varEmpty = np.zeros((dimExpect, dimExpect), dtype="float")
    sigEmpty = np.zeros((dimExpect, dimExpect), dtype="float")
    # Fill it with NaN
    imgEmpty.fill(np.nan)
    mskEmpty.fill(np.nan)
    varEmpty.fill(np.nan)
    sigEmpty.fill(np.nan)

    # Figure out the area we want, and read the data.
    # For coadds the WCS is the same in all bands, but the code handles the general case
    # Start by finding the tract and patch
    matches = skyMap.findTractPatchList(raDecList)
    tractList, patchList = cdColor.getTractPatchList(matches)
    nPatch = len(patchList)
    # Prefix of the output file
    outPre = prefix + '_' + filt + '_full'

    # Verbose
    if verbose:
        print '####################################################################'
        print " Input Ra, Dec: %10.5f, %10.5f" % (ra, dec)
        print " Cutout size is expected to be %d x %d" % (dimExpect, dimExpect)

    newX = []
    newY = []
    boxX = []
    boxY = []
    boxSize = []
    #
    trList = []
    paList = []
    zpList = []
    #
    imgArr = []
    mskArr = []
    varArr = []
    sigArr = []
    srcArr = []

    # Go through all these images
    for j in range(nPatch):
        # Tract, patch
        tract, patch = tractList[j], patchList[j]
        print "### Dealing with %d - %s" % (tract, patch)
        # Check if the coordinate is available in all three bands.
        try:
            # Get the coadded exposure
            coadd = butler.get("deepCoadd", tract=tractId,
                               patch=patchId, filter=filt, immediate=True)
        except Exception, errMsg:
            print "#########################################################"
            print " No data is available in %d - %s" % (tract, patch)
            print "#########################################################"
            print errMsg
        else:
            # Get the WCS information
            wcs = coadd.getWcs()
            if j is 0:
                # Get the CD Matrix of the WCS
                cdMatrix = wcs.getCDMatrix()
                # Get the pixel size in arcsec
                pixScale = wcs.pixelScale().asDegrees() * 3600.0
            # Convert the central coordinate from Ra,Dec to pixel unit
            pixel = wcs.skyToPixel(raDec)
            pixel = afwGeom.Point2I(pixel)
            # Define the bounding box for the central pixel
            bbox = afwGeom.Box2I(pixel, pixel)
            # Grow the bounding box to the desired size
            bbox.grow(int(size))
            xOri, yOri = bbox.getBegin()
            # Compare to the coadd image, and clip
            bbox.clip(coadd.getBBox(afwImage.PARENT))
            # Get the masked image
            try:
                subImage  = afwImage.ExposureF(coadd, bbox,
                                               afwImage.PARENT)
                # Extract the image array
                imgArr.append(subImage.getMaskedImage().getImage().getArray())
                # Extract the mask array
                mskBad = getCoaddBadMsk(subImage)
                mskArr.append(mskBad.getArray())
                # Extract the variance array
                imgVar = subImage.getMaskedImage().getVariance().getArray()
                # Convert it into sigma array
                imgSig = np.sqrt(imgVar)
                varArr.append(imgVar)
                sigArr.append(imgSig)
                # Get the source catalog
                if saveSrc:
                    """ Sometimes the forced photometry catalog might not be available """
                    try:
                        srcCat = butler.get('deepCoadd_forced_src', tract=tract,
                                            patch=patch, filter=filt,
                                            immediate=True,
                                            flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)
                        # Get the pixel coordinates for all objects
                        srcRa  = np.array(map(
                            lambda x: x.get('coord').getRa().asDegrees(), srcCat))
                        srcDec = np.array(map(
                            lambda x: x.get('coord').getDec().asDegrees(), srcCat))
                        # Simple Box match
                        indMatch = ((srcRa > (ra - sizeDegree)) &
                                    (srcRa < (ra + sizeDegree)) &
                                    (srcDec > (dec - sizeDegree)) &
                                    (srcDec < (dec + sizeDegree)))
                        # Extract the matched subset
                        srcArr.append(srcCat.subset(indMatch))
                    except:
                        print "### Tract: %d  Patch: %s" % (tractId, patchId)
                        warnings.warn("### Can not find the forced photometry catalog !")
                        srcArr.append(None)
                        if not os.path.isfile('no_src.lis'):
                            noSrc = open('no_src.lis', 'w')
                            noSrc.write("%d  %s \n" % (tractId, patchId))
                            noSrc.close()
                        else:
                            noSrc = open('no_src.lis', 'a+')
                            noSrc.write("%d  %s \n" % (tractId, patchId))
                            noSrc.close()
                # Save the width of the BBox
                boxX.append(bbox.getWidth())
                # Save the heigth of the BBox
                boxY.append(bbox.getHeight())
                # Save the size of the BBox in unit of pixels
                boxSize.append(bbox.getWidth() * bbox.getHeight())
                # New X, Y origin coordinates
                newX.append(bbox.getBeginX() - xOri)
                newY.append(bbox.getBeginY() - yOri)
                # Tract, Patch
                trList.append(tract)
                paList.append(patch)
                # Photometric zeropoint
                zpList.append(2.5 * np.log10(coadd.getFluxMag0()[0]))
            except:
                print '### SOMETHING IS WRONG WITH THIS BOUNDING BOX !!'
                print "    %d -- %s -- %s " % (tract, patch, filtArr[i])
                print "    Bounding Box Size: %d" % (bbox.getWidth() * bbox.getHeight())


    # Number of returned RGB image
    nReturn = len(newX)
    print "### Return %d Useful Images" % nReturn
    indSize = np.argsort(boxSize)
    # Go through the returned images, put them in the cutout region
    for n in range(nReturn):
        ind = indSize[n]
        # Put in the image array
        imgUse = imgArr[ind]
        imgEmpty[newY[ind]:(newY[ind] + boxY[ind]),
                 newX[ind]:(newX[ind] + boxX[ind])] = imgUse[:, :]
        # Put in the mask array
        mskUse = mskArr[ind]
        mskEmpty[newY[ind]:(newY[ind] + boxY[ind]),
                 newX[ind]:(newX[ind] + boxX[ind])] = mskUse[:, :]
        # Put in the variance array
        varUse = varArr[ind]
        varEmpty[newY[ind]:(newY[ind] + boxY[ind]),
                 newX[ind]:(newX[ind] + boxX[ind])] = varUse[:, :]
        # Put in the sigma array
        sigUse = sigArr[ind]
        sigEmpty[newY[ind]:(newY[ind] + boxY[ind]),
                 newX[ind]:(newX[ind] + boxX[ind])] = sigUse[:, :]
        if n is (nReturn - 1):
            # This is the largest available sub-image
            phoZp = zpList[ind]
            trLarge, paLarge = trList[ind], paList[ind]
            """ TODO: Get the total exposure time at the image center """

    # Save the source catalog
    if saveSrc:
        srcCount = 0
        for m in range(nReturn):
            if srcArr[m] is not None:
                if srcCount is 0:
                    srcUse = srcArr[m]
                    srcCount += 1
                else:
                    for item in srcArr[m]:
                        srcUse.append(item)
                    srcCount += 1

    # Create a WCS for the combined image
    outWcs = wcs.WCS(naxis=2)
    outWcs.wcs.crpix = [dimExpect/2.0, dimExpect/2.0]
    outWcs.wcs.crval = [ra, dec]
    outWcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    outWcs.wcs.cdelt = cdMatrix
    # Output to header
    outHead = outWcs.to_header()
    outHead.set("PIXEL", pixScale, "Pixel Scale [arcsec/pix]")
    outHead.set("PHOZP", phoZp, "Photometric Zeropoint")
    outHead.set("EXPTIME", 1.0, "Set exposure time to 1 sec")
    outHead.set("GAIN", 3.0, "Average GAIN for HSC CCDs")

    # Define the output file name
    # Save the image array
    outImg = outPre + '_img.fits'
    hduImg = fits.PrimaryHDU(imgEmpty, header=outHead)
    hduList = fits.HDUList([hduImg])
    hduList.writeto(outImg)
    # Save the mask array
    outMsk = outPre + '_bad.fits'
    hduMsk = fits.PrimaryHDU(mskEmpty, header=outHead)
    hduList = fits.HDUList([hduMsk])
    hduList.writeto(outMsk)
    # Save the variance array
    outVar = outPre + '_var.fits'
    hduVar = fits.PrimaryHDU(varEmpty, header=outHead)
    hduList = fits.HDUList([hduVar])
    hduList.writeto(outVar)
    # Save the sigma array
    outSig = outPre + '_sig.fits'
    hduSig = fits.PrimaryHDU(sigEmpty, header=outHead)
    hduList = fits.HDUList([hduSig])
    hduList.writeto(outSig)
    # If necessary, save the source catalog
    if saveSrc:
        outSrc = outPre + '_src.fits'
        srcUse.writeFits(outSrc)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("ra",   type=float, help="RA  to search")
    parser.add_argument("dec",  type=float, help="Dec to search")
    parser.add_argument("size", type=float, help="Half size of the cutout box")
    parser.add_argument('-f', '--filter', dest='filt', help="Filter",
                       default='HSC-I')
    parser.add_argument('-p', '--prefix', dest='outfile',
                        help='Prefix of the output file',
                        default='hsc_coadd_cutout')
    args = parser.parse_args()

    coaddImageCutout(args.root, args.ra, args.dec, args.size,
                     filt=args.filt, prefix=args.outfile)
