#!/usr/bin/env python

from __future__ import division

import os
import copy
import argparse
import numpy as np
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
from shapely.geometry import Polygon
from shapely.ops      import cascaded_union
from descartes        import PolygonPatch

import coaddPatchShape  as coaddPS
#getPolyUnion, polySaveWkb, polyReadWkb
import coaddPatchNoData as coaddND
#polySaveReg

""" Get a (RA, DEC) pair from WCS [NOT USED] """
def wcsGetRaDecPair(wcsInfo, xpos, ypos):

    raRad, decRad = wcsInfo.pixelToSky(xpos, ypos)
    return raRad.asDegrees(), decRad.asDegrees()

""" Get the dimension of a tract [NOT USED] """
def getTractXYDim(tractInfo):

    tractExtent = tractInfo.getBBox().getDimensions()
    return tractExtent.getX(), tractExtent.getY()

""" Get the Corner RA, DEC based on tract's WCS [NOT WORKING] """
def fourCornerRaDec(tractWcs, xDim, yDim):

    xCorners = [1, 1, xDim, xDim]
    yCorners = [1, yDim, yDim, 1]

    corners = []
    for x, y in zip(xCorners, yCorners):
        corners.append(wcsGetRaDecPair(tractWcs, x, y))

    return corners

def getTractList(rootDir, filter, imgType='deepCoadd', toInt=True,
                 prefix='hsc_coadd', toFile=True):
    """
    Get all the tractID from certain rootDir
    """

    if rootDir[-1] is not '/':
        rootDir += '/'
    tractDir = rootDir + imgType + '/' + filter + '/'

    if not os.path.exists(tractDir):
        raise Exception("Can not find the directory: %s" % tractDir)
    else:
        tractStrs = [d for d in os.listdir(tractDir) if
                     os.path.isdir(os.path.join(tractDir, d))]
        """ Save to a list file """
        if toFile:
            tractName = prefix + '_' + filter + '_tract.lis'
            tractFile = open(tractName, 'w')
            for tt in tractStrs:
                tractFile.write(tt + '\n')
            tractFile.close()
        if toInt:
            tractList = map(lambda x: int(x), tractStrs)
        else:
            tractList = tractStrs

    return tractList

""" Plot a list of polygons and its outline """
def plotListPoly(polyList, outPNG='polyList.png', outer=None,
                 color=None, minX=None, minY=None, maxX=None,
                 maxY=None, xSize=20, ySize=16, dpi=120):

    """
    Right now, have to make sure that every Polygon in the list is
    simple and valid TODO
    """

    """ Set up the color """
    BLUE = '#6699cc'
    GRAY = '#999999'
    ec = GRAY
    if color is None:
        fc = BLUE
    else:
        fc = color

    fig = plt.figure(figsize=(xSize, ySize), dpi=dpi)
    ax = fig.add_subplot(111)

    if len(polyList) == 1:
        partShow = PolygonPatch(polyList[0], fc='r', ec=GRAY,
                                alpha=0.8, zorder=1)
        ax.add_patch(partShow)
    else:
        for poly in polyList:
            partShow = PolygonPatch(poly, fc=numpy.random.rand(3,1),
                                    ec=GRAY, alpha=0.8, zorder=1)
            ax.add_patch(partShow)

    """ Draw the outline of the region """
    if outer is not None:
        if outer.type is "Polygon":
            bound = outer.boundary
            if bound.type is "LineString":
                x, y = bound.xy
                ax.plot(x, y, c='k', lw=2.5)
        elif outer.type is "MultiPolygon":
            for oo in outer:
                bound = oo.boundary
                if bound.type is "LineString":
                    x, y = bound.xy
                    ax.plot(x, y, c='k', lw=3.0)

    if (minX is None) or (minY is None) or (maxX is None) or (maxY is None):
        ax.margins(0.02, 0.02, tight=True)
    else:
        raRange  = [(minX-0.1), (maxX+0.1)]
        decRange = [(minY-0.1), (maxY+0.1)]
        ax.set_xlim(*raRange)
        ax.set_ylim(*decRange)

    ax.set_xlabel(r'RA (deg)',  fontsize=25)
    ax.set_ylabel(r'DEC (deg)', fontsize=25)

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

    fig.subplots_adjust(bottom=0.1, left=0.1,
                        top=0.98, right=0.98)

    fig.savefig(outPNG, dpi=100)
    plt.close(fig)


""" Get the (RA, DEC) pairs of the four corners for certain Tract """
def coaddTractGetCorners(skyMap, tractId):

    """ Try to get the right tract information """
    numTract = len(skyMap)
    if tractId > numTract:
        raise Exception("The tractId is not correct: %d" % tractId)

    tractInfo = skyMap[tractId]
    if tractInfo.getId() != tractId:
        raise Exception("Something wrong with the SkyMap: %d - %d" % (tractId,
                                                                      tractInfo.getId()))
    else:
        corners = tractInfo.getVertexList()
        cornerRaDec = []
        for corner in corners:
            cornerRaDec.append((corner[0].asDegrees(),
                                corner[1].asDegrees()))

    """ Form a Polygon using the cornerRaDec """
    cornerPoly = Polygon(cornerRaDec)
    if not cornerPoly.boundary.is_simple:
        raise Exception("Something is wrong with the Polygon: %d" % tractId)

    """ Return this Polygon """
    return cornerPoly

""" Main Function """
def coaddTractShape(rootDir, filter, verbose=True, prefix='hsc_tract',
                    savePNG=True, xSize=18, ySize=16):

    """ Prefix of the output file """
    preFile1 = prefix + '_' + filter + '_tract'

    """ Get the skymap """
    butler = dafPersist.Butler(rootDir)
    skyMap = butler.get("deepCoadd_skyMap", immediate=True)

    """ Get the list of tract IDs """
    tractList = getTractList(rootDir, filter, imgType='deepCoadd',
                             toInt=True, prefix=prefix, toFile=True)
    nTracts = len(tractList)
    if verbose:
        print "### Will deal with %d tracts in total" % nTracts

    polyList = []
    for tractId in tractList:
        if verbose:
            print "### Deal with tract: %d" % tractId

        """ Get the Polygon for the tract """
        tractPoly = coaddTractGetCorners(skyMap, tractId)

        preFile2 = preFile1  + '_' + str(tractId)
        """ Save a .WKB file """
        tractWkb = preFile2 + '.wkb'
        coaddPS.polySaveWkb(tractPoly, tractWkb)
        """ Save a .REG file """
        tractReg = preFile2 + '.reg'
        coaddND.polySaveReg(tractPoly, tractReg, color='green')
        """ Append to the list """
        polyList.append(tractPoly)

    if nTracts > 1:
        """ Get the cascaded_union of all the Tracts """
        combPoly = cascaded_union(polyList)
        """ Save a combined .WKB file """
        coaddPS.polySaveWkb(combPoly, preFile1 + '_all.wkb')
        """ Save a combined .REG file """
        coaddND.polySaveReg(combPoly, preFile1 + '_all.reg', color='blue')
        """ Get the bounds of the combined region """
        minX, minY, maxX, maxY = combPoly.bounds

        """ It's possible that the combined Poly is Multi-part """
        if combPoly.type is "MultiPolygon":
            combParts = combPoly.geoms[:]
            nParts = len(combParts)
            for ii in range(nParts):
                combPart = combParts[ii]
                min1, min2, max1, max2 = combPart.bounds
                """ Save .wkb and .reg file for each part """
                coaddPS.polySaveWkb(combPart, preFile1 + '_part_' + str(ii+1) + '.wkb')
                coaddND.polySaveReg(combPart, preFile1 + '_part_' + str(ii+1) + '.reg',
                                    color='blue')
                """ Make a plot """
                plotListPoly(polyList, outer=combPart,
                             outPNG=preFile1 + '_part_' + str(ii+1) + '.png',
                             minX=min1, minY=min2, maxX=max1, maxY=max2,
                             xSize=16, ySize=14)
    else:
        combPoly = polyList[0]
        minX, minY, maxX, maxY = combPoly.bounds
        coaddPS.polySaveWkb(combPoly, preFile1 + '_all.wkb')
        coaddND.polySaveReg(combPoly, preFile1 + '_all.reg', color='blue')
        polyList = [combPoly]

    """ Save a PNG file """
    if savePNG:
        pngFile = preFile1 + '_all.png'
        plotListPoly(polyList, outPNG=pngFile, outer=combPoly,
                     minX=minX, minY=minY, maxX=maxX, maxY=maxY,
                     xSize=xSize, ySize=ySize)

def batchTractShapes(rootDir, xSize=20, ySize=20, prefix='hsc_tract',
                     verbose=False, savePNG=True):

    hscFilters = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y']

    for filter in hscFilters:
        if verbose:
            print "### Deal with FILTER: %s" % filter
        coaddTractShape(rootDir, filter, verbose=verbose, prefix=prefix,
                        savePNG=savePNG, xSize=xSize, ySize=ySize)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("root",   help="Root directory of data repository")
    parser.add_argument('-p', '--prefix', dest='prefix',
                        help='Prefix of the output file',
                        default='hsc_tract')
    parser.add_argument('-x', '--xsize', dest='xsize',
                        help='X-size of the output figure for all tracts',
                        default=20)
    parser.add_argument('-y', '--ysize', dest='ysize',
                        help='Y-size of the output figure for all tracts',
                        default=18)
    args = parser.parse_args()

    batchTractShapes(args.root, prefix=args.prefix, xSize=args.xsize,
                     ySize=args.ySize)
