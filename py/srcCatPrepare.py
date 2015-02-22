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
import astropy.io.fits as fits

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

def getSrcParams(srcCat, calExp, calMd, outFits):

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

    outTab = astropy.table.Table()

    outTab.add_column(astropy.table.Column(
        name='id', data=srcCat.get('id')))
    outTab.add_column(astropy.table.Column(
        name='nchild', data=srcCat.get('deblend.nchild')))
    outTab.add_column(astropy.table.Column(
        name='extend', data=srcCat.get('classification.extendedness')))

    outTab.add_column(astropy.table.Column(
        name='Ra', data=srcRa))
    outTab.add_column(astropy.table.Column(
        name='Dec', data=srcDec))

    outTab.add_column(astropy.table.Column(
        name='psfMag', data=psfMag))
    outTab.add_column(astropy.table.Column(
        name='psfMer', data=psfMerr))
    outTab.add_column(astropy.table.Column(
        name='psfS2n', data=psfS2n))
    outTab.add_column(astropy.table.Column(
        name='kroMag', data=kroMag))
    outTab.add_column(astropy.table.Column(
        name='kroMer', data=kroMerr))
    outTab.add_column(astropy.table.Column(
        name='kroS2n', data=kroS2n))
    outTab.add_column(astropy.table.Column(
        name='gauMag', data=gauMag))
    outTab.add_column(astropy.table.Column(
        name='gauMer', data=gauMerr))
    outTab.add_column(astropy.table.Column(
        name='gauS2n', data=gauS2n))
    outTab.add_column(astropy.table.Column(
        name='expMag', data=expMag))
    outTab.add_column(astropy.table.Column(
        name='expMer', data=expMerr))
    outTab.add_column(astropy.table.Column(
        name='expS2n', data=expS2n))
    outTab.add_column(astropy.table.Column(
        name='devMag', data=devMag))
    outTab.add_column(astropy.table.Column(
        name='devMer', data=devMerr))
    outTab.add_column(astropy.table.Column(
        name='devS2n', data=devS2n))
    outTab.add_column(astropy.table.Column(
        name='modMag', data=modMag))
    outTab.add_column(astropy.table.Column(
        name='modMer', data=modMerr))
    outTab.add_column(astropy.table.Column(
        name='modS2n', data=modS2n))

    outTab.add_column(astropy.table.Column(
        name='sdssR', data=sdssR))
    outTab.add_column(astropy.table.Column(
        name='sdssBa', data=sdssBa))
    outTab.add_column(astropy.table.Column(
        name='sdssPa', data=sdssPa))
    outTab.add_column(astropy.table.Column(
        name='expR', data=expR))
    outTab.add_column(astropy.table.Column(
        name='expBa', data=expBa))
    outTab.add_column(astropy.table.Column(
        name='expPa', data=expPa))
    outTab.add_column(astropy.table.Column(
        name='devR', data=devR))
    outTab.add_column(astropy.table.Column(
        name='devBa', data=devBa))
    outTab.add_column(astropy.table.Column(
        name='devPa', data=devPa))
    outTab.add_column(astropy.table.Column(
        name='fracDev', data=srcCat.get('cmodel.fracDev')))

    outTab.add_column(astropy.table.Column(
        name='is_primary', data=srcCat.get('detect.is-primary')))
    outTab.add_column(astropy.table.Column(
        name='patch_inner', data=srcCat.get('detect.is-patch-inner')))
    outTab.add_column(astropy.table.Column(
        name='tract_inner', data=srcCat.get('detect.is-tract-inner')))

    outTab.add_column(astropy.table.Column(
        name='deblend_many_peaks', data=srcCat.get('deblend.too-many-peaks')))
    outTab.add_column(astropy.table.Column(
        name='deblend_parent_too_big', data=srcCat.get('deblend.parent-too-big')))
    outTab.add_column(astropy.table.Column(
        name='flag_badcentroid', data=srcCat.get('flags.badcentroid')))
    outTab.add_column(astropy.table.Column(
        name='flag_pix_edge', data=srcCat.get('flags.pixel.edge')))
    outTab.add_column(astropy.table.Column(
        name='flag_pix_intr', data=srcCat.get('flags.pixel.interpolated.center')))
    outTab.add_column(astropy.table.Column(
        name='flag_pix_satu', data=srcCat.get('flags.pixel.saturated.center')))
    outTab.add_column(astropy.table.Column(
        name='flag_pix_cr', data=srcCat.get('flags.pixel.cr.center')))
    outTab.add_column(astropy.table.Column(
        name='flag_pix_susp', data=srcCat.get('flags.pixel.suspect.center')))
    outTab.add_column(astropy.table.Column(
        name='flag_pix_bad', data=srcCat.get('flags.pixel.bad')))

    outTab.add_column(astropy.table.Column(
        name='shape_sdss_flags', data=srcCat.get('shape.sdss.flags')))
    outTab.add_column(astropy.table.Column(
        name='shape_sdss_unweighted', data=srcCat.get('shape.sdss.flags.unweighted')))
    outTab.add_column(astropy.table.Column(
        name='shape_sdss_unweightedbad', data=srcCat.get('shape.sdss.flags.unweightedbad')))

    outTab.add_column(astropy.table.Column(
        name='flux_gau_flags', data=srcCat.get('flux.gaussian.flags')))
    outTab.add_column(astropy.table.Column(
        name='flux_psf_flags', data=srcCat.get('flux.psf.flags')))

    outTab.add_column(astropy.table.Column(
        name='flux_kro_flags', data=srcCat.get('flux.kron.flags')))
    outTab.add_column(astropy.table.Column(
        name='flux_kro_rad', data=srcCat.get('flux.kron.radius')))
    outTab.add_column(astropy.table.Column(
        name='flux_kro_flags_radius', data=srcCat.get('flux.kron.flags.radius')))
    outTab.add_column(astropy.table.Column(
        name='flux_kro_flags_edge', data=srcCat.get('flux.kron.flags.edge')))
    outTab.add_column(astropy.table.Column(
        name='flux_kro_flags_smallRadius', data=srcCat.get('flux.kron.flags.smallRadius')))
    outTab.add_column(astropy.table.Column(
        name='flux_kro_flags_usedMinimumRadius',
        data=srcCat.get('flux.kron.flags.usedMinimumRadius')))
    outTab.add_column(astropy.table.Column(
        name='flux_kro_flags_usedPsfRadius',
        data=srcCat.get('flux.kron.flags.usedPsfRadius')))

    outTab.add_column(astropy.table.Column(
        name='cmodel_iniflux_flag', data=srcCat.get('cmodel.initial.flux.flags')))
    outTab.add_column(astropy.table.Column(
        name='cmodel_exp_fluxflag', data=srcCat.get('cmodel.exp.flux.flags')))
    outTab.add_column(astropy.table.Column(
        name='cmodel_exp_niter', data=srcCat.get('cmodel.exp.nIters')))
    outTab.add_column(astropy.table.Column(
        name='cmodel_dev_fluxflag', data=srcCat.get('cmodel.dev.flux.flags')))
    outTab.add_column(astropy.table.Column(
        name='cmodel_dev_niter', data=srcCat.get('cmodel.dev.nIters')))
    outTab.add_column(astropy.table.Column(
        name='cmodel_flux_flag', data=srcCat.get('cmodel.flux.flags')))
    outTab.add_column(astropy.table.Column(
        name='cmodel_noshape_flag', data=srcCat.get('cmodel.flags.noShape')))
    outTab.add_column(astropy.table.Column(
        name='cmodel_maxbadpixfrac_flag',
        data=srcCat.get('cmodel.flags.maxBadPixelFraction')))

    outTab.write(outFits, format='fits', overwrite=True)

def srcCatPrepare(rootDir, tract, patch, filt, prefix):

    # Make a butler and specify the dataId
    butler = dafPersist.Butler(rootDir)
    dataId = {'tract': tract, 'patch': patch, 'filter': filt}

    # Get the prefix of the output files
    prefix = prefix + '-' + str(tract) + '-' + patch + '-' + filt
    outFits = prefix + '.fits'

    # Return the src catalog, exposure, and the metadata
    srcCat, calExp, calMd = getSrcData(butler, dataId)

    # Return a numpy array of useful information
    getSrcParams(srcCat, calExp, calMd, outFits)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('root',  help="Root directory of data repository")
    parser.add_argument('tract', type=int, help="Visit to show")
    parser.add_argument('patch', help="patch to show")
    parser.add_argument('-f', '--filter', dest='filt', help="Filter",
                       default='HSC-I')
    parser.add_argument('-p', '--prefix', dest='prefix',
                        help='Prefix of the output file', default='hsc_coadd_src')
    args = parser.parse_args()

    srcCatPrepare(args.root, args.tract, args.patch, args.filt, args.prefix)
