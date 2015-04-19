#!/usr/bin/env python
# encoding: utf-8

from __future__ import division

import copy
import os
import argparse
import numpy as np

# Astropy related
from astropy.io import fits
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

# Shapely related imports
from shapely.geometry import Polygon, LineString, Point
from shapely          import wkb

def polyReadWkb(wkbName, load=True):

    wkbFile = open(wkbName, 'r')
    polyWkb = wkbFile.read().decode('hex')
    wkbFile.close()

    if load is True:
        return wkb.loads(polyWkb)
    else:
        return polyWkb

""" TODO: Need to be organized"""
def showMatchShape(wkbFile, large=None, corner=None, title='No Data Mask Plane',
                  pngName='tract_mask.png', xsize=20, ysize=18, dpi=150,
                  saveFile=True):

    fig = plt.figure(figsize=(xsize, ysize), dpi=dpi)

    ax = fig.add_subplot(111)
    fontsize = 20
    ax.minorticks_on()

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    # Set title
    ax.set_title(title, fontsize=25, fontweight='bold')
    ax.title.set_position((0.5,1.01))

    maskShow = polyReadWkb(wkbFile, load=True)
    # Outline all the mask regions
    if maskShow.type is "Polygon":
        bounds = maskShow.boundary
        if bounds.type is "LineString":
            x, y = bounds.xy
            ax.plot(x, y, c='r', lw=2.5)
        elif bounds.type is "MultiLineString":
            for bb in bounds:
                x, y = bb.xy
                ax.plot(x, y, lw=2.5, color='r')
    elif maskShow.type is "MultiPolygon":
        for ii, mask in enumerate(maskShow):
            bounds = mask.boundary
            if bounds.type is "LineString":
                x, y = bounds.xy
                ax.plot(x, y, c='r', lw=1.5)
            elif bounds.type is "MultiLineString":
                for bb in bounds:
                    x, y = bb.xy
                    ax.plot(x, y, lw=1.5, color='r')
            else:
                print " !!! Can not plot shape %d - %s !" % (ii, bounds.type)

    # highlight all the large ones
    if large is not None:
        bigShow = polyReadWkb(large, load=True)
        if bigShow.type is "Polygon":
            bounds = bigShow.boundary
            if bounds.type is "LineString":
                x, y = bounds.xy
                ax.plot(x, y, c='b', lw=2.5)
            elif bounds.type is "MultiLineString":
                for bb in bounds:
                    x, y = bb.xy
                    ax.plot(x, y, lw=2.5, color='b')
            elif bigShow.type is "MultiPolygon":
                for ii, mask in enumerate(bigShow):
                    bounds = mask.boundary
                    if bounds.type is "LineString":
                        x, y = bounds.xy
                        ax.plot(x, y, c='b', lw=2.0)
                    elif bounds.type is "MultiLineString":
                        for bb in bounds:
                            x, y = bb.xy
                            ax.plot(x, y, lw=2.0, color='b')
                    else:
                        print " !!! Can not plot shape %d - %s !" % (ii,
                                                                     bounds.type)

    # highlight all the tract corner
    if corner is not None:
        cornerShow = polyReadWkb(corner, load=True)
        if cornerShow.type is "Polygon":
            bounds = cornerShow.boundary
            if bounds.type is "LineString":
                x, y = bounds.xy
                ax.plot(x, y, c='g', lw=2.5)
            elif bounds.type is "MultiLineString":
                for bb in bounds:
                    x, y = bb.xy
                    ax.plot(x, y, lw=2.5, color='g')
            else:
                print " !!! Can not plot shape %d - %s !" % (ii, bounds.type)
        elif cornerShow.type is "MultiPolygon":
            for ii, mask in enumerate(cornerShow.geoms[:]):
                bounds = mask.boundary
                if bounds.type is "LineString":
                    x, y = bounds.xy
                    ax.plot(x, y, c='g', lw=2.5)
                elif bounds.type is "MultiLineString":
                    for bb in bounds:
                        x, y = bb.xy
                        ax.plot(x, y, lw=2.5, color='g')
                else:
                    print " !!! Can not plot shape %d - %s !" % (ii, bounds.type)
        else:
            print " !!! Not valid tract_corner Polygon"

    ax.margins(0.02, 0.02, tight=True)

    ax.set_xlabel(r'RA (deg)',  fontsize=22)
    ax.set_ylabel(r'DEC (deg)', fontsize=22)

    fig.subplots_adjust(hspace=0.1, wspace=0.1,
                        top=0.95, right=0.95)

    if saveFile:
        fig.savefig(pngName)
        plt.close(fig)
    else:
        return fig


def showMaskMatch(objMatch, acpMask, rejMask=None,
                  raField='ra', decField='dec', pngFile=None,
                  infoField1=None, infoText1=None,
                  infoField2=None, infoText2=None):
    """
    Show the matched objects on top of the masks
    """
    if pngFile is None:
        pngFile = 'object_matched.png'
    fig = showMatchShape(acpMask, large=rejMask,
                        pngName=pngFile, xsize=20, ysize=18, dpi=150,
                        saveFile=True)

def coaddMaskMatch(inCat, acpMask, raField=None, decField=None,
                   showMatch=True, infoField1=None, infoText1=None,
                   infoField2=None, infoText2=None, outCat=None,
                   rejMask=None, verbose=True):

    """ Read in the catalog """
    if not os.path.isfile(inCat):
        raise Exception("Can not find the input catalog: %s !" % inCat)

    """ Load the accept mask """
    if not os.path.isfile(acpMask):
        raise Exception("Can not find the accept mask: %s !" % acpMask)
    else:
        acpRegs = polyReadWkb(acpMask)
    """ Load the rejection mask """
    if rejMask is not None:
        if not os.path.isfile(rejMask):
            raise Exception("Can not find the rejection mask: %s !" % rejMask)
            rejRegs = None
        else:
            rejRegs = polyReadWkb(rejMask)

    """ Try to find the Ra and Dec columns """
    """ TODO: Allow different format in the future """
    catData = Table.read(inCat, format='fits')
    catCols = catData.colnames
    if verbose:
        print "### There are %d objects to be matched" % len(catData)

    """ Form an array of (RA, DEC) """
    """ TODO: Here, column name is case sensitive; should try a few combination"""
    """ Match the upper case RA, DEC with catCols"""
    if raField is None:
        raField = 'ra'
    if decField is None:
        decField = 'dec'
    try:
        raArr = catData.field(raField)
    except KeyError as key:
        raise Exception(" !!! Can not find the column for RA: %s" % raField)
    try:
        decArr = catData.field(decField)
    except KeyError as key:
        raise Exception(" !!! Can not find the column for dec: %s" % decField)

    """ Do the match """
    objInside = map(lambda x, y: acpRegs.contains(Point(x, y)), raArr, decArr)
    if rejRegs is not None:
        objMasked = map(lambda x, y: rejRegs.contains(Point(x, y)), raArr, decArr)
        objUseful = map(lambda x, y: x and (not y), objInside, objMasked)
    else:
        objUseful = objInside

    """ Extract a sub-table of matched objects """
    objMatch = catData[np.asarray(objUseful)]
    if verbose:
        print "### Returns %d matched objects" % len(objMatch)

    """ Save the matched results to a new FITS catalog """
    if outCat is None:
        prefix = os.path.splitext(inCat)[0]
        outCat = prefix + '_matched.fits'
    if os.path.isfile(outCat):
        os.remove(outCat)
    objMatch.write(outCat, format='fits')

    """ Visualize the results """
    pngFile = os.path.splitext(outCat)[0] + '.png'
    if showMatch:
        showMaskMatch(objMatch, acpMask, rejMask=rejMask,
                      raField=raField, decField=decField,
                      infoField1=infoField1, infoText1=infoText1,
                      infoField2=infoField2, infoText2=infoText2,
                      pngFile=pngFile)
    return objMatch

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("input",  help="The input FITS catalog for match")
    parser.add_argument("accept", help="'Accept' mask shows the footprint")
    parser.add_argument('-r', '--ra', dest='raField',
                        help='The column name for RA',
                        default=None)
    parser.add_argument('-d', '--dec', dest='decField',
                        help='The column name for DEC',
                        default=None)
    parser.add_argument('-o', '--output', dest='outCat',
                        help='Name of the output catalog',
                        default=None)
    args = parser.parse_args()

    coaddMaskMatch(args.input, args.accept, raField=args.raField,
                   decField=args.decField, outCat=args.outCat)
