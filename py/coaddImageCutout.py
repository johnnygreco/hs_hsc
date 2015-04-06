#!/usr/bin/env python

import argparse
import numpy                  as np
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
    # Extract specific plane from it
    mskImg &= mskImg.getPlaneBitMask(bitmask)

    return mskImg

def getCoaddBadMsk(calExp):

    # Get the mask image
    badMsk = calExp.getMaskedImage().getMask()
    # Clear the "EDGE" plane XXX TODO
    # badMsk.clearMaskPlane(4)
    # Clear the "DETECTED" plane
    badMsk.clearMaskPlane(5)
    # Clear the "DETECTED_NEGATIVE" plane
    badMsk.clearMaskPlane(6)
    # Clear the "CROSSTALK" plane XXX TODO: Check later to see if it is still
    # appropriate
    badMsk.clearMaskPlane(9)

    return badMsk

def coaddImageCutout(root, ra, dec, size, saveMsk=True, saveSrc=True,
                     filt='HSC-I', prefix='hsc_coadd_cutout'):

    # Get the SkyMap of the database
    butler = dafPersist.Butler(root)
    skyMap = butler.get("deepCoadd_skyMap", immediate=True)

    # Define the Ra, Dec pair
    point = afwGeom.Point2D(ra, dec)
    coord = afwCoord.IcrsCoord(point)

    # Search for overlapped tract, patch pairs
    # It is possible that certain coordinate is included in two patches!
    for tract, patch in skyMap.findClosestTractPatchList([coord]):

        # Get the (tract, patch) ID
        tractId = tract.getId()
        patchId = "%d,%d" % patch[0].getIndex()
        print "Find (Tract, Patch): %d, %s !" % (tractId, patchId)

        # Get the coadd images
        try:
            coadd = butler.get("deepCoadd", tract=tractId,
                               patch=patchId, filter=filt, immediate=True)

        # XXX TODO: Better handle of exception: LsstCppExceptions
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
            bbox.grow(size)
            bbox.clip(coadd.getBBox(afwImage.PARENT))
            if bbox.isEmpty():
                continue

            # Make a new ExposureF object for the cutout region
            subImage = afwImage.ExposureF(coadd, bbox, afwImage.PARENT)

            # Define the output file name
            outPre = prefix + '_' + str(tractId) + '_' + patchId + '_' + filt
            outImg = outPre + '.fits'
            outPsf = outPre + '_psf.fits'

            # Save the cutout image to a new FITS file
            subImage.writeFits(outImg)

            # Get the Coadded PSF image
            psfImg = getCoaddPsfImage(coadd, coord)
            psfImg.writeFits(outPsf)

            if saveMsk is True:
                # Get different mask planes
                mskDetec = getCoaddMskPlane(coadd, 'DETECTED')
                mskIntrp = getCoaddMskPlane(coadd, 'INTRP')
                mskSatur = getCoaddMskPlane(coadd, 'SAT')
                mskDetec.writeFits(outPre + '_detec.fits')
                mskIntrp.GwriteFits(outPre + '_intrp.fits')
                mskSatur.writeFits(outPre + '_satur.fits')

                # Get the "Bad" mask plane
                mskBad = getCoaddBadMsk(coadd)
                mskBad.writeFits(outPre + '_bad.fits')

            coaddFound = True

            if saveSrc is True:

                # Get the source catalog
                srcCat = butler.get('deepCoadd_src', tract=tractId,
                                    patch=patchId, filter=filt, immediate=True,
                                    flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)
                # Get the pixel coordinates for all objects
                srcX, srcY = srcCat.getX(), srcCat.getY()
                # Get the pixel coordinate of the cutout center
                cenX, cenY = pixel.getX(),  pixel.getY()
                # Simple Box match
                indMatch = ((srcX < (cenX - size)) & (srcX > (cenX + size)) &
                            (srcY < (cenY - size)) & (srcY < (cenY - size)))
                # Extract the matched subset
                srcMatch = srcCat.subset(indMatch)
                # Save the src catalog to a FITS file
                outSrc = outPre + '_src.fits'
                srcMatch.writeFits(outSrc)

        finally:
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
