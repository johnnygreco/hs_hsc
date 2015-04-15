#!/usr/bin/env python

from __future__ import division

import copy
import os
import argparse
import numpy as np

from astropy.io import fits
from astropy    import wcs

import lsst.daf.persistence   as dafPersist

# Matplotlib default settings
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.figsize'] = 12, 10
mpl.rcParams['xtick.major.size'] = 8.0
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.size'] = 4.0
mpl.rcParams['xtick.minor.width'] = 1.5
mpl.rcParams['ytick.major.size'] = 8.0
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.size'] = 4.0
mpl.rcParams['ytick.minor.width'] = 1.5
mpl.rc('axes', linewidth=2)

# Shapely related imports
from shapely.geometry import Polygon, LineString
from shapely          import wkb
from shapely.ops      import cascaded_union

from scipy import ndimage
from skimage.measure import find_contours, approximate_polygon

def showHscMask(coords, large=None, title='No Data Mask Plane',
                pngName=None):

    fig = plt.figure(figsize=(10, 10), dpi=120)

    ax = fig.add_subplot(111)
    fontsize = 14
    ax.minorticks_on()

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    # Set title
    ax.set_title(title, fontsize=25, fontweight='bold')
    ax.title.set_position((0.5,1.01))

    # Outline all the mask regions
    for raDec in coords:
        ax.plot(raDec[:, 1], raDec[:, 0], '-r', linewidth=1.5)

    # Using polygon to highlight all the large ones
    if large is not None:
        for raDec in large:
            ax.plot(raDec[:, 1], raDec[:, 0], '-b', linewidth=2.0)

    fig.subplots_adjust(hspace=0.1, wspace=0.1,
                        top=0.95, right=0.95)

    if pngName is not None:
        fig.savefig(pngName)
        plt.close(fig)

    return fig


def imgAddNoise(im, gaussian, factor):

    im = ndimage.gaussian_filter(im, gaussian)
    im += factor * np.random.random(im.shape)

    return im

def getPixelRaDec(wcs, xx, yy):

    coord = wcs.pixelToSky(xx, yy).toIcrs()

    ra  = coord.getRa().asDegrees()
    dec = coord.getDec().asDegrees()

    return ra, dec

def polySaveWkb(poly, wkbName):

    polyWkb = wkb.dumps(poly)

    wkbFile = open(wkbName, 'w')
    wkbFile.write(polyWkb.encode('hex'))
    wkbFile.close()


# Save the Polygon region into a DS9 .reg file
def getPolyLine(polyCoords):

    coordShow = map(lambda x: str(x[0]) + ' ' + str(x[1]) + ' ', polyCoords)

    # The format for Polygon in DS9 is:
    # Usage: polygon x1 y1 x2 y2 x3 y3 ...
    polyLine = 'polygon '
    for p in coordShow:
        polyLine += p

    polyLine += '\n'

    return polyLine

# !!! Make sure that the Polygon is simple
def polySaveReg(poly, regName, listPoly=False, color='blue'):

    # DS9 region file header
    head1 = '# Region file format: DS9 version 4.1\n'
    head2 = 'global color=%s width=2\n' % color
    head3 = 'icrs\n'
    # Open the output file and write the header
    regFile = open(regName, 'w')
    regFile.write(head1)
    regFile.write(head2)
    regFile.write(head3)

    if listPoly:
        for pp in poly:
            # Get the coordinates for every point in the polygon
            polyCoords = pp.boundary.coords[:]
            polyLine = getPolyLine(polyCoords)
            regFile.write(polyLine)
    else:
        polyCoords = poly.boundary.coords[:]
        polyLine = getPolyLine(polyCoords)
        regFile.write(polyLine)

    regFile.close()


def listAllImages(rootDir, filter):

    import glob

    if rootDir[-1] is '/':
        searchDir = rootDir + 'deepCoadd/' + filter.upper() + '/*/*.fits'
    else:
        searchDir = rootDir + '/deepCoadd/' + filter.upper() + '/*/*.fits'

    return map(lambda x: x, glob.glob(searchDir))


def coaddPatchNoData(rootDir, tract, patch, filter, prefix='hsc_coadd',
                     savePNG=True, verbose=True, tolerence=4,
                     minArea=10000):

    # Make a butler and specify the dataID
    butler = dafPersist.Butler(rootDir)
    dataId = {'tract':tract, 'patch':patch, 'filter':filter}

    # Get the name of the input fits image
    if rootDir[-1] is '/':
        fitsName = rootDir + 'deepCoadd/' + filter + '/' + str(tract).strip() \
                   + '/' + patch + '.fits'
    else:
        fitsName = rootDir + '/deepCoadd/' + filter + '/' + str(tract).strip() \
                   + '/' + patch + '.fits'
    if not os.path.isfile(fitsName):
        raise Exception('Can not find the input fits image: %s' % fitsName)
    # Read in the original image and get the wcs
    # TODO: This is not perfect
    hduList = fits.open(fitsName)
    header = hduList[1].header
    imgWcs = wcs.WCS(header)

    # Get the name of the wkb and deg file
    strTractPatch = (str(tract).strip() + '_' + patch + '_' + filter)
    ## For all the accepted regions
    noDataAllWkb = prefix + '_' + strTractPatch + '_nodata_all.wkb'
    noDataAllReg = prefix + '_' + strTractPatch + '_nodata_all.reg'
    ## For all the big mask regions
    noDataBigWkb = prefix + '_' + strTractPatch + '_nodata_big.wkb'
    noDataBigReg = prefix + '_' + strTractPatch + '_nodata_big.reg'

    # Get the name of the png file
    titlePng = prefix + strTractPatch + '_NODATA'
    noDataPng = prefix + '_' + strTractPatch + '_nodata.png'

    if verbose:
        print "## Reading Fits Image: %s" % fitsName

    # Get the exposure from the butler
    calExp = butler.get('deepCoadd', dataId, immediate=True)

    # Get the object for mask plane
    mskImg = calExp.getMaskedImage().getMask()

    # Extract the NO_DATA plane
    # TODO: NO_DATA is not a system mask, maybe should use INTRP later
    noData = copy.deepcopy(mskImg)
    noData &= noData.getPlaneBitMask('NO_DATA')
    # Return the mask image array
    noDataArr = noData.getArray()

    # Set all masked pixels to be 1
    noDataArr /= 256
    # Pad the 2-D array by a little
    noDataArr = np.lib.pad(noDataArr, ((1, 1), (1, 1)), 'constant',
                           constant_values=0)

    # Try a very different approach: Using the find_contours and
    # approximate_polygon methods from scikit-images package
    maskShapes = []  # For all the accepted mask regions
    maskCoords = []  # For the "corner" coordinates of these regions
    maskAreas  = []  # The sizes of all regions

    # Only find the 0-level contour
    contoursAll = find_contours(noDataArr, 0)
    if verbose:
        print "### %d contours have been detected" % len(contoursAll)
    for maskContour in contoursAll:
        # Approximate one extracted contour into a polygon
        # tolerance decides the accuracy of the polygon, hence
        # the number of coords for each polygon.
        # Using large tolerance also means smaller number of final
        # polygons
        contourCoords = approximate_polygon(maskContour, tolerance=tolerence)
        # Convert these coordinates into (RA, DEC) using the WCS information
        contourSkyCoords = map(lambda x: [x[1], x[0]], contourCoords)
        contourRaDec     = imgWcs.wcs_pix2world(contourSkyCoords, 1)
        # Require that any useful region must be at least an triangular
        if len(contourCoords) >= 3:
            # Form a lineString using these coordinates
            maskLine = LineString(contourRaDec)
            # Check if the lineString is valid and simple, so can be used
            # to form a closed and simple polygon
            if maskLine.is_valid and maskLine.is_simple:
                contourPoly = Polygon(contourRaDec)
                maskShapes.append(contourPoly)
                maskCoords.append(contourRaDec)
                maskAreas.append(Polygon(contourCoords).area)

    if verbose:
        print "### %d regions are useful" % len(maskAreas)
    # Save all the masked regions to a .reg file
    polySaveReg(maskShapes, noDataAllReg, listPoly=True, color='red')
    # Also create a MultiPolygon object, and save a .wkb file
    maskAll = cascaded_union(maskShapes)
    polySaveWkb(maskAll, noDataAllWkb)

    # Isolate the large ones
    maskBigList = np.array(maskShapes)[np.where(np.array(maskAreas) > minArea)]
    maskBigList = map(lambda x: x, maskBigList)
    coordBigList = np.array(maskCoords)[np.where(np.array(maskAreas) > minArea)]
    nBig = len(maskBigList)
    if nBig > 0:
        if verbose:
            print "### %d regions are larger than the minimum mask sizes" % nBig
        # Save all the masked regions to a .reg file
        polySaveReg(maskBigList, noDataBigReg, listPoly=True, color='blue')
        # Also create a MultiPolygon object, and save a .wkb file
        maskBig = cascaded_union(maskBigList)
        polySaveWkb(maskBig, noDataBigWkb)
    else:
        maskBig = None
        if verbose:
            print "### No region is larger than the minimum mask sizes"

    if savePNG:
        showHscMask(maskCoords, large=coordBigList, title=titlePng,
                    pngName=noDataPng)

    return maskAll, maskBig


def batchPatchNoData(rootDir, filter='HSC-I', prefix='hsc_coadd'):

    # Get the list of coadded images in the direction
    imgList = listAllImages(rootDir, filter)

    # Get the list of tract and patch for these images
    tract = map(lambda x: int(x.split('/')[-2]), imgList)
    patch = map(lambda x: x.split('/')[-1].split('.')[0], imgList)

    results = map(lambda x, y: coaddPatchNoData(rootDir, x, y, filter,
                                prefix=prefix, savePNG=True, verbose=True, tolerence=4,
                                minArea=10000), tract, patch)
    allList = map(lambda x: x[0], results)
    bigList = map(lambda x: x[1], results)

    allUse = []
    for ss in allList:
        if ss is not None:
            allUse.append(ss)
    bigUse = []
    for tt in bigList:
        if tt is not None:
            bigUse.append(tt)

    # Make a cascaded union of them
    allComb = cascaded_union(allUse)
    bigComb = cascaded_union(bigUse)

    # Break them down into list
    ## ALL
    if allComb.geom_type is "MultiPolygon":
        allPolys = map(lambda x: x, allComb.geoms)
        polySaveReg(allPolys, prefix + '_' + filter + '_nodata_all_combined.reg',
                    listPoly=True, color='red')
    elif allComb.geom_type is "Polygon":
        polySaveReg(allPolys, prefix + '_' + filter + '_nodata_all_combined.reg',
                    color='red')
    ## BIG
    if bigComb.geom_type is "MultiPolygon":
        bigPolys = map(lambda x: x, bigComb.geoms)
        polySaveReg(bigPolys, prefix + '_' + filter + '_nodata_big_combined.reg',
                    listPoly=True, color='blue')
    elif bigComb.geom_type is "Polygon":
        polySaveReg(bigPolys, prefix + '_' + filter + '_nodata_big_combined.reg',
                    color='blue')

    # Save these polygons as a .wkb file
    polySaveWkb(allComb, prefix + '_' + filter + '_nodata_all_combined.wkb')
    polySaveWkb(bigComb, prefix + '_' + filter + '_nodata_big_combined.wkb')

    return allComb, bigComb

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("root",   help="Root directory of data repository")
    parser.add_argument("filter", help="HSC filter")
    parser.add_argument('-p', '--prefix', dest='prefix',
                        help='Prefix of the output file',
                        default='hsc_coadd')
    args = parser.parse_args()

    batchPatchNoData(args.root, args.filter, prefix=args.prefix)
