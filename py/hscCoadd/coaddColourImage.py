#!/usr/bin/env python

from __future__ import division

import os
import argparse
import numpy   as np

import lsst.daf.persistence as dafPersist
import lsst.afw.display.rgb as afwRgb
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.image as afwImage

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import matplotlib.axes   as axes


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


def saveRgbPng(outRgb, imgRgb, cenMark=False, xCen=None, yCen=None, name=None,
               info1=None, info2=None, info3=None, sLength=None,
               sString=None):

    """
    Save the RGB image as a PNG figure
    """

    # Decide the image size
    sizeX, sizeY, dim = imgRgb.shape
    sizeX = int(sizeX / 100) if (sizeX / 100) < 15 else 15
    sizeY = int(sizeY / 100) if (sizeY / 100) < 15 else 15
    sizeX = sizeX if sizeX > 6 else 6
    sizeY = sizeY if sizeY > 6 else 6

    fig = plt.figure(figsize=(sizeX, sizeY),
                     dpi=100, frameon=False)

    # Show the image
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(imgRgb, interpolation='none')
    ax.set_aspect('equal')

    # Suppress all the ticks and tick labels
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    # Highlight the new center: TODO
    if cenMark and (xCen is not None) and (yCen is not None):
        ax.scatter(xCen, yCen, s=80, lw=0.5, marker='o', edgecolors='r',
                   facecolors='none')

    # Add some information on the image
    if name is not None:
        ax.text(0.5, 0.10, name, fontsize=20, fontweight='bold',
                ha='center', va='center', color='w',
                transform=ax.transAxes)
    if info1 is not None:
        ax.text(0.8, 0.90, info1, fontsize=18, fontweight='bold',
                ha='center', va='center', color='w',
                transform=ax.transAxes)
    if info2 is not None:
        ax.text(0.8, 0.82, info2, fontsize=18, fontweight='bold',
                ha='center', va='center', color='w',
                transform=ax.transAxes)
    if info3 is not None:
        ax.text(0.8, 0.74, info3, fontsize=18, fontweight='bold',
                ha='center', va='center', color='w',
                transform=ax.transAxes)
        if sLength is not None:
            ax.plot([0.14, 0.14+sLength], [0.88, 0.88], 'w-', lw=2.5,
                    transform=ax.transAxes)
            if sString is not None:
                ax.text((0.28 + sLength)/2.0, 0.85, sString, fontsize=15,
                        ha='center', va='center', color='w',
                        fontweight='bold', transform=ax.transAxes)

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
                     prefix='hsc_coadd_cutout', info1=None, info2=None, info3=None,
                     min=-0.0, max=0.70, Q=10, name=None, localMax=True,
                     scaleBar=10, butler=None):

    # Get the SkyMap of the database
    if butler is None:
        try:
            butler = dafPersist.Butler(root)
        except Exception:
            print '### Can not load the correct Butler!'
    skyMap = butler.get("deepCoadd_skyMap", immediate=True)

    # [Ra, Dec] pair
    raDec = afwCoord.Coord(ra*afwGeom.degrees, dec*afwGeom.degrees)

    # Expected size and center position
    dimExpect = (2 * size +1)
    cenExpect = (dimExpect/2.0, dimExpect/2.0)
    # Create a empty array
    # For RGB image, the data type should be uint8
    rgbEmpty = np.zeros((dimExpect, dimExpect, 3), dtype="uint8")

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
        images  = {}
        #newcen = {}
        newX    = {}
        newY    = {}
        cutoutSize = int(size)

        for i in range(3):

            # Find the file of the coadd image
            coadd = butler.get("deepCoadd", tract=tract, patch=patch,
                                filter=filtArr[i])

            # Get the WCS information
            wcs = coadd.getWcs()

            # Convert the central coordinate from Ra,Dec to pixel unit
            pixel = wcs.skyToPixel(raDec)
            pixel = afwGeom.Point2I(pixel)

            # Define the bounding box for the central pixel
            bbox = afwGeom.Box2I(pixel, pixel)

            # Grow the bounding box to the desired size
            bbox.grow(int(cutoutSize))
            xOri, yOri = bbox.getBegin()

            # Compare to the coadd image, and clip
            bbox.clip(coadd.getBBox(afwImage.PARENT))
            newX[i] = bbox.getBeginX() - xOri
            newY[i] = bbox.getBeginY() - yOri

            # Get the masked image
            subImage  = afwImage.ExposureF(coadd, bbox, afwImage.PARENT)
            # Get the WCS
            #subWcs = subImage.getWcs()
            # Get the central pixel coordinates on new subImage WCS
            #newcen[i] = subWcs.skyToPixel(raDec)

            # Extract the image array
            images[i] = subImage.getMaskedImage().getImage()

        # Define the Blue, Green, and Red channels
        # These cutouts are still HSC ImageF object, not numpy array
        bCut, gCut, rCut = images[0], images[1], images[2]

        if localMax:
            maxArr = []
            for m in range(3):
                imgPad = np.zeros((dimExpect, dimExpect), dtype=float)
                imgPad[newY[m]:(newY[m] + images[m].getHeight()),
                       newX[m]:(newX[m] + images[m].getWidth())] = images[m].getArray()
                globalMax = np.max(images[m].getArray())
                localMax = np.max(imgPad[cenExpect[0]-10:cenExpect[0]+10,
                                  cenExpect[1]-10:cenExpect[1]+10])
                maxArr.append(localMax / globalMax)
                #print "### %d : %6.3f" % (m+1, localMax/globalMax)
            maxShow = np.max(np.asarray(maxArr))
        else:
            maxShow = max

        #for j in range(3):
            #newCenX, newCenY = newcen[j]
            #newOriX, newOriY = images[j].getXY0()
            #newCenX = newCenX - newOriX
            #newCenY = images[j].getHeight() - (newCenY - newOriY)

        # To see if data are available for all the cut-out region
        if (bCut.getHeight() < dimExpect) or (bCut.getWidth() < dimExpect):
            print " ### Only part of the desired cutout-region is returned !"
            # Define the name of the output file
            outRgb = prefix + '_' + filt + '_part_color.png'
            partial = True
        else:
            outRgb = prefix + '_' + filt + '_color.png'
            partial = False

        # Generate the RGB image
        # 15/04/22: min ==> minimum
        imgRgb = afwRgb.makeRGB(rCut, gCut, bCut, minimum=min,
                               range=(maxShow - min), Q=Q,
                               saturatedPixelValue=None)
        if partial:
            for k in range(3):
                rgbEmpty[newY[k]:(newY[k] + images[k].getHeight()),
                        newX[k]:(newX[k] + images[k].getWidth()), k] = imgRgb[:, :, k]
            imgRgb = rgbEmpty

        # Add a scale bar
        if scaleBar is not None:
            sLength = (scaleBar * 1.0) / 0.170 / (dimExpect * 1.0)
            sString = "%d\"" % int(scaleBar)
        else:
            sLength = None
            sString = None

        # Better way to show the image
        if partial:
            saveRgbPng(outRgb, imgRgb, cenMark=True, xCen=cenExpect[0],
                       yCen=cenExpect[1], sLength=sLength, sString=sString,
                       name=name, info1=info1, info2=info2, info3=info3)
        else:
            saveRgbPng(outRgb, imgRgb,  name=name,
                       info1=info1, info2=info2, info3=info3,
                       sLength=sLength, sString=sString)

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
                         prefix='hsc_coadd_cutout',
                         info1=None, info2=None, info3=None,
                         min=-0.0, max=0.70, Q=10, name=None, localMax=True,
                         scaleBar=10, butler=None, verbose=False):

    # Get the SkyMap of the database
    if butler is None:
        try:
            butler = dafPersist.Butler(root)
            if verbose:
                print "### Load in the Butler"
        except Exception:
            print '### Can not load the correct Butler!'
    skyMap = butler.get("deepCoadd_skyMap", immediate=True)

    # [Ra, Dec] list
    raDec = afwCoord.Coord(ra*afwGeom.degrees, dec*afwGeom.degrees)
    raList, decList = getCircleRaDec(ra, dec, size)
    points = map(lambda x, y: afwGeom.Point2D(x, y), raList, decList)
    raDecList = map(lambda x: afwCoord.IcrsCoord(x), points)

    # Expected size and center position
    dimExpect = (2 * size +1)
    cenExpect = (dimExpect/2.0, dimExpect/2.0)
    # Create a empty array
    # For RGB image, the data type should be uint8
    rgbEmpty = np.zeros((dimExpect, dimExpect, 3), dtype="uint8")

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
    # Cutout size
    cutoutSize = int(size)

    #Figure out the area we want, and read the data.
    #For coadds the WCS is the same in all bands, but the code handles the general case
    #Start by finding the tract and patch
    matches = skyMap.findTractPatchList(raDecList)
    tractList, patchList = getTractPatchList(matches)
    nPatch = len(patchList)
    # Output RGB image
    if verbose:
        print "### WILL DEAL WITH %d (TRACT, PATCH)" % nPatch
    outRgb = prefix + '_' + filt + '_color.png'

    newX = []
    newY = []
    boxX = []
    boxY = []
    boxSize = []
    rgbArr = []
    # Go through all these images
    for j in range(nPatch):
        # Tract, patch
        tract, patch = tractList[j], patchList[j]
        if verbose:
            print "### Dealing with %d - %s" % (tract, patch)
        # Check if the coordinate is available in all three bands.
        try:
            # Get the metadata
            md1 = butler.get("deepCoadd_md", immediate=True,
                            tract=tract, patch=patch, filter=filter1)
            md2 = butler.get("deepCoadd_md", immediate=True,
                            tract=tract, patch=patch, filter=filter2)
            md3 = butler.get("deepCoadd_md", immediate=True,
                            tract=tract, patch=patch, filter=filter3)
        except Exception, errMsg:
            print "#########################################################"
            print " The galaxy is not available in %d - %s" % (tract, patch)
            print "#########################################################"
            print errMsg
        else:
            #Then we can read the desired pixels
            images  = {}
            # Go through the three bands
            for i in range(3):
                # Find the file of the coadd image
                coadd = butler.get("deepCoadd", tract=tract, patch=patch,
                                    filter=filtArr[i], immediate=True)
                # Get the WCS information
                wcs = coadd.getWcs()
                # Convert the central coordinate from Ra,Dec to pixel unit
                pixel = wcs.skyToPixel(raDec)
                pixel = afwGeom.Point2I(pixel)
                # Define the bounding box for the central pixel
                bbox = afwGeom.Box2I(pixel, pixel)
                # Grow the bounding box to the desired size
                bbox.grow(int(cutoutSize))
                xOri, yOri = bbox.getBegin()
                # Compare to the coadd image, and clip
                bbox.clip(coadd.getBBox(afwImage.PARENT))
                # Get the masked image
                try:
                    subImage  = afwImage.ExposureF(coadd, bbox,
                                                   afwImage.PARENT)
                    # Extract the image array
                    images[i] = subImage.getMaskedImage().getImage()
                    bboxGood = True
                    if i == 1:
                        boxX.append(bbox.getWidth())
                        boxY.append(bbox.getHeight())
                        boxSize.append(bbox.getWidth() * bbox.getHeight())
                        newX.append(bbox.getBeginX() - xOri)
                        newY.append(bbox.getBeginY() - yOri)
                except:
                    print '### SOMETHING IS WRONG WITH THIS BOUNDING BOX !!'
                    print "    %d -- %s -- %s " % (tract, patch, filtArr[i])
                    print "    Bounding Box Size: %d" % (bbox.getWidth() * bbox.getHeight())
                    bboxGood  = False
                    continue

            if bboxGood:
                # Define the Blue, Green, and Red channels
                # These cutouts are still HSC ImageF object, not numpy array
                bCut, gCut, rCut = images[0], images[1], images[2]
                # Generate the RGB image
                # 15/04/22: min ==> minimum
                imgRgb = afwRgb.makeRGB(rCut, gCut, bCut, minimum=min,
                                       range=(max - min), Q=Q,
                                       saturatedPixelValue=None)
                rgbArr.append(imgRgb)
            else:
                continue

    # Number of returned RGB image
    nReturn = len(rgbArr)
    if verbose:
        print "### Return %d Useful Images" % nReturn
    if nReturn > 0:
        indSize = np.argsort(boxSize)
        # Go through the returned images, put them in the cutout region
        for n in range(nReturn):
            ind = indSize[n]
            # This could lead to problem FIXME
            rgbUse = rgbArr[ind]
            for k in range(3):
                rgbEmpty[newY[ind]:(newY[ind] + boxY[ind]),
                         newX[ind]:(newX[ind] + boxX[ind]), k] = rgbUse[:, :, k]

        imgRgb = rgbEmpty
        # Add a scale bar
        if scaleBar is not None:
            sLength = (scaleBar * 1.0) / 0.170 / (dimExpect * 1.0)
            sString = "%d\"" % int(scaleBar)
        else:
            sLength = None
            sString = None
        # Better way to show the image
        saveRgbPng(outRgb, imgRgb, name=name,
                   info1=info1, info2=info2, info3=info3,
                   sLength=sLength, sString=sString)

    else:
        print "### NO COLOR IMAGE IS GENERATED FOR THIS OBJECT !!"


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

#def coaddColourImageFull(root, ra, dec, size, filt='gri',
                         #prefix='hsc_coadd_cutout',
                         #info1=None, info2=None, info3=None,
                         #min=-0.0, max=0.70, Q=10, name=None, localMax=True,
                         #scaleBar=10, butler=None, verbose=False):
