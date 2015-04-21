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


def saveRgbPng(outRgb, imgRgb, xCen=None, yCen=None, name=None, info=None):

    """
    Save the RGB image as a PNG figure
    TODO: Need more works!
    """

    import matplotlib.pyplot as plt
    import matplotlib.axes   as axes

    fig = plt.figure(dpi=120, frameon=False)

    # Show the image
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(imgRgb, interpolation='none')
    ax.set_aspect('equal')

    # Suppress all the ticks and tick labels
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    # Highlight the new center: TODO
    if (xCen is not None) and (yCen is not None):
        ax.scatter(xCen, yCen, s=60, lw=2.0, marker='o', edgecolors='r',
                   facecolors='none')

    # Add some information on the image
    if name is not None:
        ax.text(0.5, 0.10, name, fontsize=20, fontweight='bold',
                ha='center', va='center', color='w',
                transform=ax.transAxes)
    if info is not None:
        ax.text(0.8, 0.88, info, fontsize=18, fontweight='bold',
                ha='center', va='center', color='w',
                transform=ax.transAxes)

    ax.margins(0.00, 0.00, tight=True)

    fig.savefig(outRgb, bbox_inches='tight', pad_inches=0, ad_inches=0)
    plt.close(fig)


def isHscFilter(filter, short=True):

    if not short:
        hscFilters = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y']
    else:
        hscFilters = ['g', 'r', 'i', 'z', 'y']

    return (filter in hscFilters)


def coaddColourImage(root, ra, dec, size, filt='gri',
                     prefix='hsc_coadd_cutout', info=None,
                     min=-0.0, max=0.70, Q=10, name=None):

    # Get the SkyMap of the database
    butler = dafPersist.Butler(root)
    skyMap = butler.get("deepCoadd_skyMap", immediate=True)

    # [Ra, Dec] pair
    raDec = afwCoord.Coord(ra*afwGeom.degrees, dec*afwGeom.degrees)

    # Expected size and center position
    dimExpect = (2 * size +1)
    cenExpect = (dimExpect/2.0, dimExpect/2.0)

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
        newcen = {}
        height = {}
        width  = {}
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

            # Compare to the coadd image, and clip
            bbox.clip(coadd.getBBox(afwImage.PARENT))

            # Get the masked image
            subImage  = afwImage.ExposureF(coadd, bbox, afwImage.PARENT)
            # Get the WCS
            subWcs = subImage.getWcs()

            # Get the central pixel coordinates on new subImage WCS
            newcen[i] = subWcs.skyToPixel(raDec)

            height[i] = bbox.getHeight()
            width[i]  = bbox.getWidth()

            images[i] = subImage.getMaskedImage()

        # Define the Blue, Green, and Red channels
        B, G, R = images[0].getImage(), images[1].getImage(), images[2].getImage()

        for j in range(3):
            newCenX, newCenY = newcen[j]
            newOriX, newOriY = images[j].getImage().getXY0()
            newCenX = newCenX - newOriX
            newCenY = height[j] - (newCenY - newOriY)
        newX = newCenX
        newY = newCenY

        # To see if data are available for all the cut-out region
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
        saveRgbPng(outRgb, imgRgb, xCen=newX, yCen=newY, name=name, info=info)

    #return imgRgb

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
