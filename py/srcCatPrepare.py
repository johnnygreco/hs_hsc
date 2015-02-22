#!/usr/bin/env python
"""
srcCatPrepare.py
Generate a copy of the src-HSC catalog that is easier to read
"""

import lsst.daf.persistence as dafPersist
from lsst.afw.table import SourceCatalog, SchemaMapper
import lsst.pex.exceptions
import lsst.afw.geom as afwGeom
import numpy as np
import argparse
import re
import collections
import astropy.table

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

def getS2n(fluxArr, ferrArr):
    """
    Get the S/N of the detections
    TODO: Is this the right way
    """
    return (fluxArr / ferrArr)

def getEllipse(quad):
    """
    Returns the semi-major axis, axes ratio and PA for a given quadrupole moment
    """
    e = afwGeom.ellipses.Axes(quad)
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

def getSrcParams(srcCat, calExp, calMd):

    # Get the WCS and pixel scale
    wcs   = calExp.getWcs()
    # Get the calibration information
    calib = calExp.getCalib()
    # don't raise an exception when we encounter a negative or NaN flux
    calib.setThrowOnNegativeFlux(False)

    # Get the magnitudes for different methods of photometry
    # Get the PSF magnitude and its error
    psfMag, psfMerr = calib.getMagnitude(srcCat.getPsfFlux(),
                                         srcCat.getPsfFluxErr())
    psfS2n = getS2n(srcCat.getPsfFlux(), srcCat.getPsfFluxErr())
    # Get the Kron magnitude and its error
    kroMag, kroMerr = calib.getMagnitude(srcCat.get('flux.kron'),
                                         srcCat.get('flux.kron.err'))
    kroS2n = getS2n(srcCat.get('flux.kron'), srcCat.get('flux.kron.err'))
    # Get the Gaussian magnitude and its error
    gauMag, gauMerr = calib.getMagnitude(srcCat.get('flux.gaussian'),
                                         srcCat.get('flux.gaussian.err'))
    gauS2n = getS2n(srcCat.get('flux.gaussian'),
                    srcCat.get('flux.gaussian.err'))
    # Get the ExpModel magnitude and its error
    expMag, expMerr = calib.getMagnitude(srcCat.get('cmodel.exp.flux'),
                                         srcCat.get('cmodel.exp.flux.err'))
    expS2n = getS2n(srcCat.get('cmodel.exp.flux'),
                    srcCat.get('cmodel.exp.flux.err'))
    # Get the DevModel magnitude and its error
    devMag, devMerr = calib.getMagnitude(srcCat.get('cmodel.dev.flux'),
                                         srcCat.get('cmodel.dev.flux.err'))
    devS2n = getS2n(srcCat.get('cmodel.dev.flux'),
                    srcCat.get('cmodel.dev.flux.err'))
    # Get the CmodelModel magnitude and its error
    modMag, modMerr = calib.getMagnitude(srcCat.get('cmodel.flux'),
                                         srcCat.get('cmodel.flux.err'))
    modS2n = getS2n(srcCat.get('cmodel.flux'),
                    srcCat.get('cmodel.flux.err'))

    # Get the (Ra, Dec) and shape information for each records
    srcRa    = []
    srcDec   = []
    sdssR    = []
    sdssBa   = []
    sdssPa   = []
    expR     = []
    expBa    = []
    expPa    = []
    devR     = []
    devBa    = []
    devPa    = []

    for record in srcCat:

        # Get the (Ra, Dec) coordinates of the objects
        coord = record.get('coord')
        srcRa.append(coord.getRa().asDegrees())
        srcDec.append(coord.getDec().asDegrees())
        # Get the SDSS Shape
        rr1, ba1, pa1 = getEllipse(record.get('shape.sdss'))
        sdssR.append(rr1)
        sdssBa.append(ba1)
        sdssPa.append(pa1)
        # Get the exp Shape
        rr2, ba2, pa2 = getEllipse(record.get('cmodel.exp.ellipse'))
        expR.append(rr2)
        expBa.append(ba2)
        expPa.append(pa2)
        # Get the dev Shape
        rr3, ba3, pa3 = getEllipse(record.get('cmodel.dev.ellipse'))
        devR.append(rr3)
        devBa.append(ba3)
        devPa.append(pa3)

    srcRa   = np.array(srcRa)
    srcDec  = np.array(srcDec)
    sdssR   = np.array(sdssR)
    sdssBa  = np.array(sdssBa)
    sdssPa  = np.array(sdssPa)
    expR    = np.array(expR)
    expBa   = np.array(expBa)
    expPa   = np.array(expPa)
    devR    = np.array(devR)
    devBa   = np.array(devBa)
    devPa   = np.array(devPa)

    srcParams = np.array([(srcCat.get('id')), (srcCat.get('parent')),
                          (srcCat.get('deblend.nchild')),
                          (srcCat.get('classification.extendedness')),
                          (srcRa), (srcDec),
                          (psfMag), (psfMerr), (psfS2n),
                          (kroMag), (kroMerr), (kroS2n),
                          (gauMag), (gauMerr), (gauS2n),
                          (expMag), (expMerr), (expS2n),
                          (devMag), (devMerr), (devS2n),
                          (modMag), (modMerr), (modS2n),
                          (sdssR), (sdssBa), (sdssPa),
                          (expR), (expBa), (expPa),
                          (devR), (devBa), (devPa),
                          (srcCat.get('cmodel.fracDev')),
                          (srcCat.get('detect.is-patch-inner')),
                          (srcCat.get('detect.is-tract-inner')),
                          (srcCat.get('detect.is-primary'))
                          ],
                         dtype=[('id', int), ('parent', int),
                                ('nchild', int),
                                ('extended', int),
                                ('ra', float), ('dec', float),
                                ('psfMag', float), ('psfMerr', float),
                                ('psfS2n', float),
                                ('kroMag', float), ('kroMerr', float),
                                ('kroS2n', float),
                                ('gauMag', float), ('gauMerr', float),
                                ('gauS2n', float),
                                ('expMag', float), ('expMerr', float),
                                ('expS2n', float),
                                ('devMag', float), ('devMerr', float),
                                ('devS2n', float),
                                ('modMag', float), ('modMerr', float),
                                ('modS2n', float),
                                ('sdssR',  float), ('sdssBa', float),
                                ('sdssPa', float),
                                ('expR',  float), ('expBa', float),
                                ('expPa', float),
                                ('devR',  float), ('devBa', float),
                                ('devPa', float),
                                ('fracDev', float),
                                ('patch_inner', bool), ('tract_inner', bool),
                                ('is_primary',  bool)
                               ])

    return srcParams


def srcCatPrepare(rootDir, tract, patch, filt, prefix):

    # Make a butler and specify the dataId
    butler = dafPersist.Butler(rootDir)
    dataId = {'tract': tract, 'patch': patch, 'filter': filt}

    # Get the prefix of the output files
    prefix = prefix + '-' + str(tract) + '-' + patch + '-' + filt

    # Return the src catalog, exposure, and the metadata
    srcCat, calExp, calMd = getSrcData(butler, dataId)

    # Return a numpy array of useful information
    srcParams = getSrcParams(srcCat, calExp, calMd)

    print srcParams[0]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('root',  help="Root directory of data repository")
    parser.add_argument('tract', type=int, help="Visit to show")
    parser.add_argument('patch', help="patch to show")
    parser.add_argument('-f', '--filter', dest='filt', help="Filter",
                       default='HSC-I')
    parser.add_argument('-p', '--prefix', dest='outfile',
                        help='Prefix of the output file', default='hsc_coadd_src')
    args = parser.parse_args()

    srcCatPrepare(args.root, args.tract, args.patch, args.filt, args.outfile)
