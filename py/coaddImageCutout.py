#!/usr/bin/env python

import copy
import argparse
import numpy as np
import lsst.daf.persistence   as dafPersist
import lsst.afw.coord         as afwCoord
import lsst.afw.image         as afwImage
import lsst.afw.geom          as afwGeom
import lsst.afw.table         as afwTable

def getCoaddPsfImage(calExp, coord):

    # Get the WCS information
    wcs = calExp.getWcs()
    # The X,Y coordinate of the image center
    coordXY = wcs.skyToPixel(coord)
    # Get the PSF object for the exposure
    psf = calExp.getPsf()
    psfImg = psf.computeImage(coordXY)

    return psfImg

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
                     circleMatch=True):

    # Get the SkyMap of the database
    butler = dafPersist.Butler(root)
    skyMap = butler.get("deepCoadd_skyMap", immediate=True)

    # Get the expected cutout size
    dimExpect  = (2 * size + 1)
    sizeExpect = dimExpect ** 2
    # Cutout size in unit of degree
    sizeDeg = size * 0.168 / 3600.0

    # Verbose
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
    print "### Find %d possible matches !" % nPatch

    matchCen = []
    for tract, patch in matches:

        # Get the (tract, patch) ID
        tractId = tract.getId()
        patchId = "%d,%d" % patch[0].getIndex()
        print "### Choose (Tract, Patch) for center: %d, %s !" % (tractId, patchId)
        matchCen.append((tractId, patchId))

        # Get the coadd images
        # Try to load the coadd Exposure; the skymap covers larger area than the
        # available data, which will cause Butler to fail sometime
        try:
            coadd = butler.get("deepCoadd", tract=tractId,
                               patch=patchId, filter=filt, immediate=True)

        except Exception, errMsg:

            coaddFound = False
            print "#############################################"
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
            bbox.clip(coadd.getBBox(afwImage.PARENT))

            if bbox.isEmpty():
                continue
            else:
                if bbox.getArea() < sizeExpect:
                    partialCut = True
                else:
                    partialCut = False

            # Make a new ExposureF object for the cutout region
            subImage = afwImage.ExposureF(coadd, bbox, afwImage.PARENT)

            # To see if data are available for all the cut-out region
            if partialCut:
                print "### Only part of the desired cutout-region is returned !"
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

            # Get the Coadded PSF image
            psfImg = getCoaddPsfImage(coadd, coord)
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

            coaddFound = True

            if saveSrc is True:

                # Get the source catalog
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

    # If only part of the desired cutout region is covered, and the circleMatch
    # Flag is set, find all the Patches that overlap with a circle region around
    # the input Ra, Dec
    if partialCut and circleMatch:

        print "####### Search for other overlapped patches #######"

        # Return the list of RA, DEC that described a circle region around the
        # input (RA, DEC).  The radius is the input size in unit of arcsec
        raList, decList = getCircleRaDec(ra, dec, (size * 0.7))
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
        print "### Find %d possible overlap patches" % nPatch

        for tract, patch in matches:

            # Example all possible matches
            tractId = tract.getId()
            for ii in range(len(patch)):
                patchId = "%d,%d" % patch[ii].getIndex()

                # Skip the image that has been used
                if (tractId, patchId) in matchCen:
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
                    bbox.grow(int(size * 1.1))
                    bbox.clip(coadd.getBBox(afwImage.PARENT))

                    if bbox.isEmpty():
                        continue
                    elif bbox.getArea() < int(sizeExpect * 0.1):
                        # Ignore small overlapped image
                        print "### %d - %s has very small overlapped region" % (tractId, patchId)
                        continue
                    else:
                        print "### Find one useful overlap: %d, %s" % (tractId,
                                                                       patchId)

                    # Make a new ExposureF object for the cutout region
                    subImage = afwImage.ExposureF(coadd, bbox, afwImage.PARENT)

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

    return coaddFound


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
