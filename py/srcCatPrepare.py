#!/usr/bin/env python
"""
srcCatPrepare.py
Generate a copy of the src-HSC catalog that is easier to read
"""

import lsst.daf.persistence as dafPersist
from lsst.afw.table import SourceCatalog, SchemaMapper
import lsst.afw.geom
import lsst.pex.exceptions
import numpy as np
import argparse
import re
import collections
import astropy.table
import lsst.afw.geom.ellipses

def getMag(flux, fluxerr, zeropoint):
    """
    Return the magnitude and error
    """
    mag, magerr = -2.5 * np.log10(flux), 2.5/np.log(10.0)*fluxerr/flux
    return (mag.T + zeropoint).T, magerr

def getMagBatch(calExp, fluxArr, ferrArr=None):
    """
    Return the magnitude and error for a list of flux
    """
    calCalib = calExp.getCalib()
    calCalib.setThrowOnNegativeFlux(False)
    if ferrArr is None:
        magArr = calCalib.getMagnitude(fluxArr)
        return magArr
    else:
        magArr, merrArr = calCalib.getMagnitude(fluxArr, ferrArr)
        return magArr, merrArr

def getS2N(fluxArr, ferrArr):
    """
    Get the S/N of the detections
    TODO: Is this the right way
    """
    return (fluxArr / ferrArr)

def getEllipse(quad):
    """
    Returns the semi-major axis, axes ratio and PA for a given quadrupole moment
    """
    e = lsst.afw.geom.ellipses.Axes(quad)
    return e.getA(), e.getB()/e.getA(), e.getTheta() * 180.0/np.pi


def getSrcData(butler, dataId):
    """
    Get the src catalog, exposure, metadata from the butler
    """
    srcCat = butler.get('deepCoadd_src', dataId, immediate=True,
                       flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS)
    calExp = butler.get('deepCoadd',    dataId)
    calMd  = butler.get('deepCoadd_md', dataId)
    return srcCat, calExp, calMd


def srcCatPrepare(rootDir, tract, patch, filt):

    # Make a butler and specify the dataId
    butler = dafPersist.Butler(rootDir)
    dataId = {'tract': tract, 'patch': patch, 'filter': filt}
    # Return the src catalog, exposure, and the metadata
    srcCat, calExp, calMd = getSrcData(butler, dataId)

    # Get the WCS and pixel scale




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('root',  help="Root directory of data repository")
    parser.add_argument('tract', type=int, help="Visit to show")
    parser.add_argument('patch', help="patch to show")
    parser.add_argument('-f', '--filter', dest='filt', help="Filter",
                       default='HSC-I')
    args = parser.parse_args()

    srcCatPrepare(args.root, args.tract, args.patch, args.filt)
