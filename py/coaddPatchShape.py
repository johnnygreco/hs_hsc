#!/usr/bin/env python

from __future__ import division

import numpy as np
import os

# Shapely related imports
from shapely.geometry import Polygon
from shapely.geometry import MultiPoint
from shapely.geometry import Point
from shapely.ops      import cascaded_union
from shapely          import wkb

from descartes        import PolygonPatch

# Astropy related imports
from astropy.io    import fits
from astropy.table import Table

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

BLUE = '#6699cc'
GRAY = '#999999'

# Read a fits catalog and return its data
# TODO: Moved to coaddUtils.py
def readFitsCat(catFile, hdu=1):

    ## Open and read the catalog file

    hduList = fits.open(catFile)
    catData = hduList[hdu].data

    return catData

# Construct a list of Polygons using the coordinates of the four corners
def getPolyList(catData):

    # Get the polygon shape for each image

    polyList = map(lambda x: MultiPoint(((x.field('raLL'), x.field('decLL')),
                                     (x.field('raUL'), x.field('decUL')),
                                     (x.field('raUR'), x.field('decUR')),
                                     (x.field('raLR'), x.field('decLR')))).convex_hull, catData)
    return polyList


# Get the cascaded union of a list of Polygongs
def getPolyUnion(polyList):

    # Get the combined poly region

    polyUnion = cascaded_union(polyList)
    polyArea  = polyUnion.area

    return polyUnion, polyArea

# Return the min/max Ra/DEC for input catalog
def getMinMaxRaDec(catData):

    # Get the min, max of the Ra, Dec
    minRa  = np.min([catData.field('raLL').min(), catData.field('raLR').min(),
                    catData.field('raUL').min(), catData.field('raUR').min()])
    minDec = np.min([catData.field('decLL').min(), catData.field('decLR').min(),
                    catData.field('decUL').min(), catData.field('decUR').min()])
    maxRa  = np.max([catData.field('raLL').max(), catData.field('raLR').max(),
                    catData.field('raUL').max(), catData.field('raUR').max()])
    maxDec = np.max([catData.field('decLL').max(), catData.field('decLR').max(),
                    catData.field('decUL').max(), catData.field('decUR').max()])

    return minRa, maxRa, minDec, maxDec

# Separate the UnionPoly into a multiple parts polygons
def unionBreakDown(unionPoly):

    # Break the multi-part polygon down to separated regions

    nUnionPart = len(unionPoly.geoms)
    if nUnionPart is 1:
        return unionPoly, nUnionPart
    else:
        unionParts = map(lambda x: x, unionPoly.geoms)
        return unionParts, nUnionPart

# Whether certain RA/DEC is included in certain Polygon
# Should work for UnionPoly even when there are separated parts
def isRaDecInside(poly, ra, dec):

    # See whether certain coordinate is inside the polygon

    point = Point((ra, dec))
    isInside = point.within(poly)

    return isInside

# Match a list of RA/DEC to a .wkb file
# Return a list of True/False
def wkbListMatch(wkbFile, raList, decList):

    if not os.path.isfile(wkbFile):
        raise Exception('Can not find the .wkb file: %s') % wkbFile
    else:
        polyWkb = polyReadWkb(wkbFile)
        matches = np.array(map(lambda x,y : isRaDecInside(polyWkb, x, y),
                               raList, decList))
        return matches

# Plot a polygon region
def plotSinglePoly(poly, outPNG='poly.png'):

    # Only plot a single Polygon region

    fig = plt.figure(1, figsize=(14, 10), dpi=100)
    ax = fig.add_subplot(111)

    partShow = PolygonPatch(poly, fc=BLUE, ec=GRAY,
                            alpha=0.8, zorder=1)
    ax.add_patch(partShow)

    polyBounds = poly.bounds
    raRange  = [(polyBounds[0]-0.2), (polyBounds[2]+0.2)]
    decRange = [(polyBounds[1]-0.2), (polyBounds[3]+0.2)]

    ax.set_xlim(*raRange)
    ax.set_ylim(*decRange)

    ax.set_xlabel(r'RA\ (deg)',  fontsize=22, fontweight='bold')
    ax.set_ylabel(r'DEC\ (deg)', fontsize=22, fontweight='bold')

    fontsize = 16
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    ax.minorticks_on()
    plt.tick_params(which='major', width=2.0, length=8.0, labelsize=20)
    plt.tick_params(which='minor', width=1.8, length=6.0)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2.5)

    ax.grid(alpha=0.6, color='k', linewidth=1.5)

    fig.savefig(outPNG, dpi=100)


# WKB File related
# Save a Polygon region to a .wkb file
def polySaveWkb(poly, wkbName):

    polyWkb = wkb.dumps(poly)

    wkbFile = open(wkbName, 'w')
    wkbFile.write(polyWkb.encode('hex'))
    wkbFile.close()

# Read a .wkb file into a Polygon shape
def polyReadWkb(wkbName, load=True):

    wkbFile = open(wkbName, 'r')
    polyWkb = wkbFile.read().decode('hex')
    wkbFile.close()

    if load is True:
        return wkb.loads(polyWkb)
    else:
        return polyWkb

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
def polySaveReg(poly, regName, listPoly=False):

    # DS9 region file header
    head1 = '# Region file format: DS9 version 4.1\n'
    head2 = 'global color=blue width=2\n'
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


# Given a list of Polygon shapes, find the common region among them
def polyMultiCommon(polyMulti):

    nPoly = len(polyMulti)
    outerMulti = cascaded_union(polyMulti).bounds

    # Multiple intersections
    if nPoly > 2:
        interMulti = polyMulti[0].intersection(polyMulti[1])
        for pp in polyMulti[2:]:
            interMulti = interMulti.intersection(pp)
    elif nPoly == 2:
        interMulti = polyMulti[0].intersection(polyMulti[1])
    else:
        print("WARNING: only one polygon found!")
        interMulti = polyMulti[0]

    return interMulti, outerMulti

# Main function
# Read a Corner coordinates catalog, and return the UnionPoly region for all the
# images contained.  If the UnionPoly has separated parts, break them down, and
# save individual regions.
def coaddPatchShape(catFile, savePng=True):

    # Define the name of the output wkb file
    catPrefix = os.path.splitext(catFile)[0]

    # Check the cat file
    if not os.path.isfile(catFile):
        raise Exception('Can not find the input catlog: %s' % catFile)
    # Read in the data
    catData = readFitsCat(catFile)

    # Get the list of polygons
    polyList = getPolyList(catData)
    # Save a DS9 .reg file for each individual images
    polySaveReg(polyList, catPrefix + '_list.wkb', listPoly=True)

    # Get the Union of all polygons and its total area
    polyUnion, areaUnion = getPolyUnion(polyList)
    # Save a region file
    wkbUnion = polySaveWkb(polyUnion, catPrefix + '.wkb')

    # Break the Union into separated regions if necessary
    polyParts, nParts = unionBreakDown(polyUnion)
    if nParts > 1:
        for ii in range(nParts):
            # Save a separate .wkb file
            partWkb = catPrefix + '_' + str(ii+1).strip() + '.wkb'
            polySaveWkb(polyParts[ii], partWkb)
            # Save a separated .reg file
            partReg = catPrefix + '_' + str(ii+1).strip() + '.reg'
            polySaveReg(polyParts[ii], partReg)
            # Save a PNG figure if necessary
            if savePng:
                partPng = catPrefix + '_' + str(ii+1).strip() + '.png'
                plotSinglePoly(polyParts[ii], outPNG=partPng)
    else:
        # Save a separated .reg file
        polyReg = catPrefix + '.reg'
        polySaveReg(polyUnion, polyReg)
        # Save a PNG file if necessary
        if savePng:
            polyPng = catPrefix + '.png'
            plotSinglePoly(polyUnion, outPNG=polyPng)

    return polyList, polyUnion


# Read in a list of .wkb files and store the common region
def coaddCommonArea(wkbList, clobber=True, prefix='poly_common',
                    savePng=True):

    # Number of Wkb files
    nWkb = len(wkbList)
    if nWkb == 0:
        raise Exception("Can not find any matched WKB file!")
    else:
        print "### %d wkb files to be dealt with!" % nWkb

    # Get a list of Polygons from these wkb files
    polyMulti = map(lambda x: polyReadWkb(x), wkbList)
    # Get the overlapped region and the outer boundary of all the Polygons
    interMulti, outerMulti = polyMultiCommon(polyMulti)
    # Save them to .wkb files
    polySaveWkb(interMulti, prefix + '_inter.wkb')
    polySaveWkb(outerMulti, prefix + '_outer.wkb')
    # If the interMulti has separated regions, break it down too:
    commonParts, nCommon = unionBreakDown(interMulti)
    if nCommon > 1:
        for ii in range(nCommon):
            # Save a separate .wkb file
            commonWkb = prefix + '_' + str(ii+1).strip() + '_inter.wkb'
            polySaveWkb(commonParts[ii], commonWkb)
            # Save a separated .reg file
            commonReg = prefix + '_' + str(ii+1).strip() + '_inter.reg'
            polySaveReg(commonParts[ii], commonReg)
            # Save a PNG figure if necessary
            if savePng:
                commonPng = prefix + '_' + str(ii+1).strip() + '_inter.png'
                plotSinglePoly(commonParts[ii], outPNG=commonPng)
    else:
        # Save a separated .reg file
        commonReg = prefix + '_inter.reg'
        polySaveReg(interMulti, commonReg)
        # Save a PNG file if necessary
        if savePng:
            commonPng = prefix + '_inter.png'
            plotSinglePoly(interMulti, outPNG=commonPng)

    return commonParts, nCommon

# Search for corner catalog using certain pattern, and run coaddPatchShape on
# every one of them
def batchPatchShape(location, pattern, clobber=False):

    #location = '.'
    #pattern  = 'ssp*_corners.fits'
    # Check and update the location and pattern
    if location[-1] is not '/':
        location += '/'
    pattern = location + pattern

    import glob
    # Get the list of corner catalogs
    cornerList = glob.glob(pattern)
    nCorner = len(cornerList)
    if nCorner == 0:
        raise Exception("Can not find any matched catalog!")
    else:
        print "### %d corner catalogs to be dealt with!" % nCorner
    # Get the list of Wkb files
    wkbList = map(lambda x: x.replace('.fits', '.wkb'), cornerList)

    for ii in range(nCorner):

        print "### Deal with catalog: %s" % cornerList[ii]
        if (not os.path.isfile(wkbList[ii])) or clobber:
            print "### %s is not found! Generate it!" % wkbList[ii]
            polyList, polyUnion = coaddPatchShape(cornerList[ii])
        else:
            print "### %s is available" % wkbList[ii]

        if not os.path.isfile(wkbList[ii]):
            raise Exception("Can not generate .wkb file %s" % wkbList[ii])

    # Get the regions that are covered by all five bands
    import re
    commonPrefix = wkbList[0].replace('.fits', '')
    commonPrefix = re.sub(r"HSC-._", "", commonPrefix)
    print "### Prefix for the common region is %s" % commonPrefix
    commonParts, nCommon = coaddCommonArea(wkbList, clobber=True,
                                           prefix=commonPrefix)

    return wkbList
