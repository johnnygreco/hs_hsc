#!/usr/bin/env python

from __future__ import division

import copy
import os
import argparse
import numpy as np

# Astropy related
from astropy.io import fits

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

#
from coaddPatchNoData import polyReadWkb

def showMaskMatch(objMatch, acpMask, rejMask=rejMask,
                      infoField1=None, infoText1=None,
                  inforField2=None, infoText2=None):
    """
    Show the matched objects on top of the masks
    """


def coaddMaskMatch(inCat, acpMask, raField=None, decField=None,
                   showMatch=True, infoField1=None, infoText1=None,
                   infoField2=None, infoText2=None, outCat=None,
                   rejMask=None):

    """ Read in the catalog """
    if not os.path.isfile(inCat):
        raise Exception("Can not find the input catalog: %s !" % inCat)

    """ Load the accept mask """
    if not os.path.isfile(acpMask):
        raise Exception("Can not find the accept mask: %s !" % acpMask)
    """ Load the rejection mask """
    if rejMask is not None:
        if not os.path.isfile(rejMask):
            raise Exception("Can not find the rejection mask: %s !" % rejMask)

    """ Try to find the Ra and Dec columns """
    catHdu = fits.open(inCat)
    catData = catHdu[1].data

    """ Form an array of (RA, DEC) """

    """ Do the match """

    """ Extract a sub-table of matched objects """

    """ Save the matched results to a new FITS catalog """

    """ Visualize the results """
    if showMatch:
        showMaskMatch(objMatch, acpMask, rejMask=rejMask,
                      infoField1=infoField1, infoText1=infoText1,
                      infoField2=infoField2, infoText2=infoText2)

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
