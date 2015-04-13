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
from matplotlib import pyplot
from matplotlib import rc, font_manager

BLUE = '#6699cc'
GRAY = '#999999'

def readFitsCat(catFile, hdu=1):

    ## Open and read the catalog file

    hduList = fits.open(catFile)
    catData = hduList[hdu].data

    return catData


def getPolyList(catData):

    # Get the polygon shape for each image

    polyList = map(lambda x: MultiPoint(((x.field('raLL'), x.field('decLL')),
                                     (x.field('raUL'), x.field('decUL')),
                                     (x.field('raUR'), x.field('decUR')),
                                     (x.field('raLR'), x.field('decLR')))).convex_hull, catData)
    return polyList


def getPolyUnion(polyList):

    # Get the combined poly region

    polyUnion = cascaded_union(polyList)
    polyArea  = polyUnion.area

    return polyUnion, polyArea


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


def unionBreakDown(unionPoly):

    # Break the multi-part polygon down to separated regions

    nUnionPart = len(unionPoly.geoms)
    if nUnionPart is 1:
        return unionPoly, nUnionPart
    else:
        unionParts = map(lambda x: x, unionPoly.geoms)
        return unionParts, nUnionPart

# Match related
def isRaDecInside(poly, ra, dec):

    # See whether certain coordinate is inside the polygon

    point = Point((ra, dec))
    isInside = point.within(poly)

    return isInside

def wkbListMatch(wkbFile, raList, decList):

    if not os.path.isfile(wkbFile):
        raise Exception('Can not find the .wkb file: %s') % wkbFile
    else:
        polyWkb = polyReadWkb(wkbFile)
        matches = np.array(map(lambda x,y : isRaDecInside(polyWkb, x, y),
                               raList, decList))
        return matches

# Plotting related
def plotPolyPart(poly, index, outPNG='poly.png'):

    # Only plot a part of the union

    fig = pyplot.figure(1, figsize=(14, 10), dpi=100)
    ax = fig.add_subplot(111)

    polyParts = poly.geoms
    partShow = PolygonPatch(polyParts[index], fc=BLUE, ec=GRAY,
                            alpha=0.8, zorder=1)
    ax.add_patch(partShow)

    polyBounds = polyParts[index].bounds
    raRange  = [(polyBounds[0]-0.2), (polyBounds[2]+0.2)]
    decRange = [(polyBounds[1]-0.2), (polyBounds[3]+0.2)]

    ax.set_xlim(*raRange)
    ax.set_ylim(*decRange)

    #plt.rc('text', usetex=True)
    ax.set_xlabel(r'RA\ (deg)',  fontsize=22, fontweight='bold')
    ax.set_ylabel(r'DEC\ (deg)', fontsize=22, fontweight='bold')

    fontsize = 16
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    ax.minorticks_on()
    pyplot.tick_params(which='major', width=2.0, length=8.0, labelsize=20)
    pyplot.tick_params(which='minor', width=1.8, length=6.0)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2.5)

    ax.grid(alpha=0.6, color='k', linewidth=1.5)


    fig.savefig(outPNG, dpi=100)


# WKB File related
def polySaveWkb(poly, wkbName):

    polyWkb = wkb.dumps(poly)

    wkbFile = open(wkbName, 'w')
    wkbFile.write(polyWkb.encode('hex'))
    wkbFile.close()

def polyReadWkb(wkbName, load=True):

    wkbFile = open(wkbName, 'r')
    polyWkb = wkbFile.read().decode('hex')
    wkbFile.close()

    if load is True:
        return wkb.loads(polyWkb)
    else:
        return polyWkb

# Multiple catalog related
def polyMultiCommon(polyMulti):

    nPoly = len(polyMulti)
    boundsMulti = cascaded_union(polyMulti).bounds

    # Multiple intersections TODO: Don't know if there is better way to do this?
    if nPoly > 2:
        interMulti = polyMulti[0].intersection(polyMulti[1])
        for pp in polyMulti[2:]:
            interMulti = interMulti.intersection(pp)
    elif nPoly == 2:
        interMulti = polyMulti[0].intersection(polyMulti[1])
    else:
        print("WARNING: only one polygon found!")
        interMulti = polyMulti[0]

    return interMulti, boundsMulti

# TODO
def coaddCommonArea(catList):


# Main function
def coaddPatchShape(catFile):

    # Define the name of the output wkb file
    catPrefix = os.path.splitext(catFile)[0]

    # Check the cat file
    if not os.path.isfile(catFile):
        raise Exception('Can not find the input catlog: %s' % catFile)
    # Read in the data
    catData = readFitsCat(catFile)

    # Get the list of polygons
    polyList = getPolyList(catData)

    # Get the Union of all polygons and its total area
    polyUnion, areaUnion = getPolyUnion(polyList)
    # Save a region file
    wkbUnion = polySaveWkb(polyUnion, catPrefix + '.wkb')

    # Break the Union into separated regions if necessary
    polyParts, nParts = unionBreakDown(polyUnion)
    if nParts > 1:
        for ii in range(nParts):
            polySaveWkb(polyParts[ii], catPrefix + '_' + str(ii+1).strip() + '.wkb')

    return polyList, polyUnion


def batchPatchShape(location, pattern):

    #location = '.'
    #pattern  = 'ssp*_corners.fits'

    for fileName in os.listdir(location):
        if fnmatch.fnmatch(fileName, pattern):
            print ('## ' + fileName)
            polyList, polyUnion = coaddPatchShape(fileName)


