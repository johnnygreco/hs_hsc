#!/usr/bin/env python

import os
import argparse
import numpy   as np
import lsst.daf.persistence as dafPersist
import lsst.afw.display.rgb as afwRgb
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.image as afwImage

def getCircleRaDec(ra, dec, size):

    # Get a small but representative set of (RA, DEC) that describe a circle
    # region around the central input coordinate

    # Convert the size from pixel unit to degress
    sizeDegree = (size * 0.168) / 3600.0
    # representative set of polar angles
    angles = np.array([0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0])
    phi = np.array(angles * np.pi / 180.0)

    # Get the (Ra, Dec) coordinates of these points
    raList  = ra  + sizeDegree * np.cos(phi)
    decList = dec + sizeDegree * np.sin(phi)

    # Also include the center
    raList  = np.append(raList,  ra)
    decList = np.append(decList, dec)

    return raList, decList

def saveRgbPng(outRgb, imgRgb):

    import matplotlib.pyplot as plt

    fig = plt.figure(dpi=120)
    fig.imshow(imgRgb, interpolation='none')
    fig.savefig(outRgb)

def isHscFilter(filter, full=True):

    if full:
        hscFilters = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y']
    else:
        hscFilters = ['g', 'r', 'i', 'z', 'y']

    return (filter in hscFilters)

def coaddColourImage(root, ra, dec, size, filt='gri',
                     prefix='hsc_coadd_cutout', info=None,
                     min=-0.0, max=0.6, Q=4):

    # Get the SkyMap of the database
    butler = dafPersist.Butler(root)
    skyMap = butler.get("deepCoadd_skyMap", immediate=True)

    # [Ra, Dec] pair
    raDec = afwCoord.Coord(ra*afwGeom.degrees, dec*afwGeom.degrees)

    # Check the choice of filters
    if len(filt) is not 3:
        raise Exception("Have to be three filters!")
    elif not (isHscFilter(filt[0]) & isHscFilter(filt[1]) &
              isHscFilter(filt[2])):
        raise Exception("Not all filters are valid !")

    #Figure out the area we want, and read the data.
    #For coadds the WCS is the same in all bands, but the code handles the general case
    #Start by finding the tract and patch
    tractInfo, patchInfo = skyMap.findTractPatchList([raDec])[0]
    tract = tractInfo.getId()
    patch = "%d,%d" % patchInfo[0].getIndex()

    # Check if the coordinate is available in all three bands.
    try:

        # Get the correct HSC filter name
        filter1 = "HSC-%s" % filt[0].upper()
        filter2 = "HSC-%s" % filt[1].upper()
        filter3 = "HSC-%s" % filt[2].upper()

        # Get the metadata
        md1 = butler.get("deepCoadd_md", immediate=True,
                        tract=tract, patch=patch, filter=filter1)
        md2 = butler.get("deepCoadd_md", immediate=True,
                        tract=tract, patch=patch, filter=filter2)
        md3 = butler.get("deepCoadd_md", immediate=True,
                        tract=tract, patch=patch, filter=filter3)

        filtArr = [filter1, filter2, filter3]

    except Exception, errMsg:

        print "#############################################"
        print " The desired coordinate is not available !!! "
        print "#############################################"
        print errMsg

    else:

        #Then we can read the desired pixels
        images = {}
        wcs    = {}
        cutoutSize = int(size)

        for i in range(3):

            # Find the file of the coadd image
            coadd = butler.get("deepCoadd", immediate=True,
                                  tract=tract, patch=patch, filter=filtArr[i])

            # Get the WCS information
            wcs = coadd.getWcs()

            # Convert the central coordinate from Ra,Dec to pixel unit
            pixel = wcs.skyToPixel(raDec)
            pixel = afwGeom.Point2I(pixel)

            # Define the bounding box for the central pixel
            bbox = afwGeom.Box2I(pixel, pixel)

            # Grow the bounding box to the desired size
            bbox.grow(int(cutoutSize))
            bbox.clip(coadd.getBBox(afwImage.PARENT))

            # Get the masked image
            subImage  = afwImage.ExposureF(coadd, bbox, afwImage.PARENT)
            # Get the WCS
            subWcs = subImage.getWcs()

            images[i] = subImage.getMaskedImage()
            wcs[i] = subWcs

        # Define the Blue, Green, and Red channels
        # TODO: Need to be tested
        B, G, R = images[0].getImage(), images[1].getImage(), images[2].getImage()
        cenB, cenG, cenR = wcs[0].skyToPixel(raDec), wcs[1].skyToPixel(raDec), \
                wcs[2].skyToPixel(raDec)

        # Get the image coordinate for the input (RA, DEC)
        # TODO

        # To see if data are available for all the cut-out region
        dimExpect = (2 * size +1)
        if (B.getHeight < dimExpect) or (B.getWidth() < dimExpect):
            print " ### Only part of the desired cutout-region is returned !"
            # Define the name of the output file
            outRgb = prefix + '_' + filt + '_part_color.png'
        else:
            outRgb = prefix + '_' + filt + '_color.png'

        # Generate the RGB image
        imgRgb = afwRgb.makeRGB(R, G, B, min=min, range=(max - min), Q=Q,
                               saturatedPixelValue=None)

        # Better way to show the image
        # TODO
        saveRgbPng(outRgb, imgRgb)
        # afwRgb.writeRGB(outRgb, imgRgb)

    return imgRgb

def getTractPatchList(matches):

    tract = []
    patch = []

    for match in matches:
        tractInfo, patchInfo = match
        tractId = tractInfo.getId()
        for patchItem in patchInfo:
            tract.append(tractId)
            patch.append("%d,%d" % patchItem.getIndex())

    return tract, patch

def getFitsImgName(root, tract, patch, filter, imgType='deepCoadd'):

    if root[-1] is not '/':
        root += '/'

    imgName = root + imgType + '/' + filter + '/' + str(tract) + '/' + patch + '.fits'

    return imgName


def coaddColourImageFull(root, ra, dec, size, filt='gri',
                        prefix='hsc_coadd_cutout', info=None,
                        min=-0.0, max=0.6, Q=4):

    # Get the SkyMap of the database
    butler = dafPersist.Butler(root)
    skyMap = butler.get("deepCoadd_skyMap", immediate=True)

    # [Ra, Dec] pair
    raDec = afwCoord.Coord(ra*afwGeom.degrees, dec*afwGeom.degrees)
    # [Ra, Dec] List
    raList, decList = getCircleRaDec(ra, dec, size*0.8)
    points = map(lambda x, y: afwGeom.Point2D(x, y), raList, decList)
    raDec  = map(lambda x: afwCoord.IcrsCoord(x), points)

    # Check the choice of filters
    if len(filt) is not 3:
        raise Exception("Have to be three filters!")
    elif not (isHscFilter(filt[0]) & isHscFilter(filt[1]) &
              isHscFilter(filt[2])):
        raise Exception("Not all filters are valid !")
    # Get the correct HSC filter name
    filter1 = "HSC-%s" % filt[0].upper()
    filter2 = "HSC-%s" % filt[1].upper()
    filter3 = "HSC-%s" % filt[2].upper()
    filtArr = [filter1, filter2, filter3]

    cutoutSize = int(size)

    #Figure out the area we want, and read the data.
    #For coadds the WCS is the same in all bands, but the code handles the general case
    #Start by finding the tract and patch
    #tractInfo, patchInfo = skyMap.findClosestTractPatchList(raDec)
    matches = skyMap.findClosestTractPatchList(raDec)
    # Return a list of tracts and patches
    tracts, patches = getTractPatchList(matches)
    # Number of unique tracts (?)
    nTract = len(matches)
    nPatch = len(patches)
    print "### Match returns %d tract(s) and %d patch(es) !" % (nTract, nPatch)

    for tract, patch in zip(tracts, patches):

        for filter in filtArr:
            fitsName = getFitsImgName(root, tract, patch, filter,
                                      imgType='deepCoadd')
            if not os.path.isfile():
                print " XXX Can not find image: %s" % fitsName
            else:
                coadd = butler.get("deepCoadd", immediate=True,
                                   tract=tract, patch=patch, filter=filter)
                # Get the WCS information
                wcs = coadd.getWcs()
                # Convert the central coordinate from Ra,Dec to pixel unit
                pixel = wcs.skyToPixel(raDec)
                pixel = afwGeom.Point2I(pixel)
                # Define the bounding box for the central pixel
                bbox = afwGeom.Box2I(pixel, pixel)
                # Grow the bounding box to the desired size
                bbox.grow(int(cutoutSize))
                bbox.clip(coadd.getBBox(afwImage.PARENT))
                # Get the masked image
                subImage  = afwImage.ExposureF(coadd, bbox, afwImage.PARENT)
                # Get the WCS
                subWcs = subImage.getWcs()

        ## TODO: To be continued ...
            images[i] = subImage.getMaskedImage()
            wcs[i] = subWcs

        # Define the Blue, Green, and Red channels
        # TODO: Need to be tested
        B, G, R = images[0].getImage(), images[1].getImage(), images[2].getImage()
        cenB, cenG, cenR = wcs[0].skyToPixel(raDec), wcs[1].skyToPixel(raDec), \
                wcs[2].skyToPixel(raDec)

        # Get the image coordinate for the input (RA, DEC)
        # TODO

        # To see if data are available for all the cut-out region
        dimExpect = (2 * size +1)
        if (B.getHeight < dimExpect) or (B.getWidth() < dimExpect):
            print " ### Only part of the desired cutout-region is returned !"
            # Define the name of the output file
            outRgb = prefix + '_' + filt + '_part_color.png'
        else:
            outRgb = prefix + '_' + filt + '_color.png'

        # Generate the RGB image
        imgRgb = afwRgb.makeRGB(R, G, B, min=min, range=(max - min), Q=Q,
                               saturatedPixelValue=None)

        # Better way to show the image
        # TODO
        saveRgbPng(outRgb, imgRgb)
        # afwRgb.writeRGB(outRgb, imgRgb)

    return imgRgb


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("ra",   type=float, help="RA  to search")
    parser.add_argument("dec",  type=float, help="Dec to search")
    parser.add_argument("size", type=float, help="Half size of the cutout box")
    parser.add_argument('-f', '--filters', dest='filt',
                        help="Combination of three filters for the colour image",
                        default='gri')
    parser.add_argument('-p', '--prefix', dest='outfile',
                        help='Prefix of the output file',
                        default='hsc_coadd_cutout')
    parser.add_argument('-i', '--info', dest='info',
                        help='Information to show on the image',
                        default=None)
    args = parser.parse_args()

    coaddColourImage(args.root, args.ra, args.dec, args.size,
                     filt=args.filt, prefix=args.outfile,
                     info=args.info)
