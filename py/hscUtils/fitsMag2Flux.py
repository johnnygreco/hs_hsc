#!/usr/bin/env python
# encoding: utf-8
"""Convert the HSC magnitudes in a FITS table to fluxes."""

from __future__ import division, print_function

import os
import copy
import argparse
import numpy as np

from astropy.table import Table, Column

# Personal
# import galSBP
import hscUtils as hUtil

HSC_FILTERS_SHORT = ['g', 'r', 'i', 'z', 'y']
HSC_FILTERS = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-I', 'HSC-Y']
"""
Extinction correction factor for HSC
A\_lambda = Coeff * E(B-V)
"""
A_G = 3.233
A_R = 2.291
A_I = 1.635
A_Z = 1.261
A_Y = 1.076
"""
For output
"""
COM = '#' * 100
SEP = '-' * 100
WAR = '!' * 100


def hscFlux2AB(flux, zero=27.0):
    """
    Convert HSC flux in unit of ADU to AB magnitude.

    So far, constant zeropoint is applied to the calibration
    """
    return (-2.5 * np.log10(flux) + zero)


def hscMag2Flux(mag, unit='maggy'):
    """
    Convert HSC AB magnitude into physical flux.

    Three units can be used here:
    unit='maggy/nanomaggy/jy'
    """
    flux = 10.0 ** (-0.4 * mag)

    if unit.lower().strip() == 'jy':
        return (flux * 3631.0)
    elif unit.lower().strip() == 'maggy':
        return flux
    elif unit.lower().strip() == 'nanomaggy':
        return (flux * 1.0E-9)
    else:
        raise Exception("## Wrong unit, should be jy/maggy/nanomaggy")


def hscMaggy2AB(flux):
    """Convert flux in unit of Maggies into AB magnitude."""
    return (np.log10(flux) / -0.4)


def hscMaggyErr2ABErr(flux, fluxErr, ivar=False):
    """Convert (flux, fluxErr) into AB magnitude error."""
    if ivar:
        fluxErr = np.sqrt(1.0 / fluxErr)

    return (2.5 * np.log10((flux + fluxErr) / flux))


def hscMagerr2Ivar(flux, magErr):
    """Get the inverse variance of flux estimates from Flux and magErr."""
    fluxErr = flux * ((10.0 ** (magErr/2.5)) - 1.0)

    return (1.0 / (fluxErr ** 2.0))


def hscMagerr2Fluxerr(flux, magErr):
    """Get the inverse variance of flux estimates from Flux and magErr."""
    fluxErr = flux * ((10.0 ** (magErr/2.5)) - 1.0)

    return fluxErr


def hscFluxSNR2Ivar(flux, snr):
    """Estimate inverse variance of flux error using HSC flux and SNR."""
    fluxErr = flux * snr

    return (1.0 / (fluxErr ** 2.0))


def getLuminosity(mag, redshift, extinction=None,
                  amag_sun=None):
    """Get the absolute magnitude or luminosity."""
    distmod = hUtil.cosmoDistMod(redshift)
    absMag = (mag - distmod)
    if extinction is not None:
        absMag -= extinction
    if amag_sun is not None:
        absMag = ((amag_sun - absMag) / 2.5)

    return absMag


def getExtinction(ra, dec, a_lambda=None):
    """
    Estimate the Galactic extinction for HSC filters.

    Parameters:
        ra, dec : The input coordinates can be arrays
    """
    # Check the input, if it's scalar, convert to array
    if not hasattr(ra, "__len__") or not hasattr(dec, "__len__"):
        ra = np.asrray([ra])
        dec = np.asarray([dec])
    if len(ra) != len(dec):
        raise Exception("## The RA and DEC should contain " +
                        "same number of elements!")
    # First try mwdust from Jo Bovy
    try:
        import mwdust
        sfd = mwdust.SFD(sf10=True)

        from astropy.coordinates import SkyCoord
        coords = SkyCoord(ra, dec, frame='icrs', unit='deg')
        galactic = coords.galactic
        l, b = galactic.l, galactic.b
        ebv = sfd(l, b, 0)
    except ImportError:
        try:
            # Then try sncosmo
            from sncosmo import SFD98Map
            dustDir = os.environ.get('DUST_DIR')
            if (not os.path.isfile(os.path.join(dustDir,
                                   'SFD_dust_4096_ngp.fits'))) or (
                not os.path.isfile(os.path.join(dustDir,
                                   'SFD_dust_4096_sgp.fits'))):
                print('# DUST_DIR : %s' % dustDir)
                raise Exception("# Can not find the SFD dust map!")
            else:
                sfd = SFD98Map(dustDir)
                ebv = sfd.get_ebv((ra, dec))
        except ImportError:
            raise Exception("# Both mwdust and sncosmo are not available")
    if a_lambda is not None:
        return (ebv * a_lambda)
    else:
        return ebv


def run(incat, idCol='ID', magType='cmodel', snType='kron', redCol='Z',
        raCol='RA', decCol='DEC', agCol='a_g', arCol='a_r', aiCol='a_i',
        azCol='a_z', ayCol='a_y', snError=True, verbose=True, extCor=True):
    """
    Read in an input catalog, convert the HSC magnitudes to useful fluxes.

    Parameters:
    """
    if not os.path.isfile(incat):
        raise Exception("## Can not find the input catalog : %s" % incat)
    else:
        inTab = Table.read(incat, format='fits')
        """Name of the output table"""
        outCat = incat.replace('.fits', ('_flux_' + magType + '.fits'))
        outTab = copy.deepcopy(inTab)
        colNames = inTab.colnames
        empty = np.zeros(len(incat))
        if verbose:
            print(SEP)
            print("## There are %d galaxies in the input catalog" % len(incat))
    """Check Extinction Correction"""
    if extCor:
        if ((agCol not in colNames) or (arCol not in colNames) or
           (aiCol not in colNames) or (azCol not in colNames) or
           (ayCol not in colNames)):
            if verbose:
                print(WAR)
                print("## Can not find information for Galactic Extinction!")
            try:
                ebv = getExtinction(inTab[raCol], inTab[decCol])
                aG = ebv * A_G
                aR = ebv * A_R
                aI = ebv * A_I
                aZ = ebv * A_Z
                aY = ebv * A_Y
                outTab.add_column(Column(name='a_g', data=aG))
                outTab.add_column(Column(name='a_r', data=aR))
                outTab.add_column(Column(name='a_i', data=aI))
                outTab.add_column(Column(name='a_z', data=aZ))
                outTab.add_column(Column(name='a_y', data=aY))
            except Exception:
                print(WAR)
                print("## Galactic Extinctions are not corrected!!")
                aG, aR, aI, aZ, aY = empty, empty, empty, empty, empty
        else:
            aG, aR, aI = inTab[agCol], inTab[arCol], inTab[aiCol]
            aZ, aY = inTab[azCol], inTab[ayCol]
    else:
        aG, aR, aI, aZ, aY = empty, empty, empty, empty, empty
    """Magnitude for photometry"""
    gCol = 'gmag_' + magType.strip()
    rCol = 'rmag_' + magType.strip()
    iCol = 'imag_' + magType.strip()
    zCol = 'zmag_' + magType.strip()
    yCol = 'ymag_' + magType.strip()
    """Convert magnitudes to fluxes"""
    fluxG = hscMag2Flux(inTab[gCol] - aG, unit='maggy')
    fluxR = hscMag2Flux(inTab[rCol] - aR, unit='maggy')
    fluxI = hscMag2Flux(inTab[iCol] - aI, unit='maggy')
    fluxZ = hscMag2Flux(inTab[zCol] - aZ, unit='maggy')
    fluxY = hscMag2Flux(inTab[yCol] - aY, unit='maggy')
    """Add columns for fluxes"""
    Maggies = np.dstack((fluxG, fluxR, fluxI, fluxZ, fluxY))[0]
    outTab.add_column(Column(name='Maggies', data=Maggies))

    """Error of fluxes"""
    if not snError:
        gErrCol = gCol + '_err'
        rErrCol = rCol + '_err'
        iErrCol = iCol + '_err'
        zErrCol = zCol + '_err'
        yErrCol = yCol + '_err'
        ivarG = hscMagerr2Ivar(fluxG, inTab[gErrCol])
        ivarR = hscMagerr2Ivar(fluxR, inTab[rErrCol])
        ivarI = hscMagerr2Ivar(fluxI, inTab[iErrCol])
        ivarZ = hscMagerr2Ivar(fluxZ, inTab[zErrCol])
        ivarY = hscMagerr2Ivar(fluxY, inTab[yErrCol])
    else:
        sigmaG = hscMag2Flux(inTab['gmag_' + snType] - aG, unit='maggy')
        sigmaR = hscMag2Flux(inTab['rmag_' + snType] - aR, unit='maggy')
        sigmaI = hscMag2Flux(inTab['imag_' + snType] - aI, unit='maggy')
        sigmaZ = hscMag2Flux(inTab['zmag_' + snType] - aZ, unit='maggy')
        sigmaY = hscMag2Flux(inTab['ymag_' + snType] - aY, unit='maggy')
        snrG = sigmaG / hscMagerr2Fluxerr(sigmaG,
                                          inTab['gmag_' + snType + '_err'])
        snrR = sigmaR / hscMagerr2Fluxerr(sigmaR,
                                          inTab['rmag_' + snType + '_err'])
        snrI = sigmaI / hscMagerr2Fluxerr(sigmaI,
                                          inTab['imag_' + snType + '_err'])
        snrZ = sigmaZ / hscMagerr2Fluxerr(sigmaZ,
                                          inTab['zmag_' + snType + '_err'])
        snrY = sigmaY / hscMagerr2Fluxerr(sigmaY,
                                          inTab['ymag_' + snType + '_err'])
        ivarG = hscFluxSNR2Ivar(fluxG, snrG)
        ivarR = hscFluxSNR2Ivar(fluxR, snrR)
        ivarI = hscFluxSNR2Ivar(fluxI, snrI)
        ivarZ = hscFluxSNR2Ivar(fluxZ, snrZ)
        ivarY = hscFluxSNR2Ivar(fluxY, snrY)
    """Add columns for invariance of fluxes"""
    IvarMaggies = np.dstack((ivarG, ivarR, ivarI, ivarZ, ivarY))[0]
    outTab.add_column(Column(name='Ivars', data=IvarMaggies))

    """Also get the absolute magnitudes"""
    if redCol in colNames:
        absmagG = getLuminosity(inTab[gCol], inTab[redCol],
                                extinction=aG)
        absmagR = getLuminosity(inTab[rCol], inTab[redCol],
                                extinction=aR)
        absmagI = getLuminosity(inTab[iCol], inTab[redCol],
                                extinction=aI)
        absmagZ = getLuminosity(inTab[zCol], inTab[redCol],
                                extinction=aZ)
        absmagY = getLuminosity(inTab[yCol], inTab[redCol],
                                extinction=aY)
        outTab.add_column(Column(name=('absG_' + magType),
                                 data=absmagG))
        outTab.add_column(Column(name=('absR_' + magType),
                                 data=absmagR))
        outTab.add_column(Column(name=('absI_' + magType),
                                 data=absmagI))
        outTab.add_column(Column(name=('absZ_' + magType),
                                 data=absmagZ))
        outTab.add_column(Column(name=('absY_' + magType),
                                 data=absmagY))
    else:
        print(WAR)
        print("## Can not find columns of redshift!")

    """Save the output catalog"""
    outTab.write(outCat, format='fits', overwrite=True)

    return outTab


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("incat", help="Input catalog")
    parser.add_argument('--id', '--idCol', dest='idCol',
                        help='Column for galaxy ID',
                        default='ID')
    parser.add_argument('--m', '--magType', dest='magType',
                        help='Type of magnitude for photometry',
                        default='cmodel')
    parser.add_argument('--sn', '--snType', dest='snType',
                        help='Type of magnitude for S/N',
                        default='kron')
    parser.add_argument('--z', '--redCol', dest='redCol',
                        help='Column for redshift',
                        default='Z')
    parser.add_argument('--ra', '--raCol', dest='raCol',
                        help='Column for RA',
                        default='RA')
    parser.add_argument('--dec', '--decCol', dest='decCol',
                        help='Column for DEC',
                        default='DEC')
    parser.add_argument('--ag', '--agCol', dest='agCol',
                        help='Column for extinction of HSC-G',
                        default='a_g')
    parser.add_argument('--ar', '--arCol', dest='arCol',
                        help='Column for extinction of HSC-R',
                        default='a_r')
    parser.add_argument('--ai', '--aiCol', dest='aiCol',
                        help='Column for extinction of HSC-I',
                        default='a_i')
    parser.add_argument('--az', '--azCol', dest='azCol',
                        help='Column for extinction of HSC-Z',
                        default='a_z')
    parser.add_argument('--ay', '--ayCol', dest='ayCol',
                        help='Column for extinction of HSC-Y',
                        default='a_y')
    parser.add_argument('--se', '--snError', dest='snError',
                        action="store_true",
                        default=True)
    parser.add_argument('--ext', '--extCor', dest='extCor',
                        action="store_true",
                        default=False)
    parser.add_argument('-v', '--verbose', dest='verbose',
                        action="store_true",
                        default=True)

    args = parser.parse_args()

    run(args.incat, idCol=args.idCol, magType=args.magType, snType=args.snType,
        zCol=args.zCol, raCol=args.raCol, decCol=args.decCol,
        agCol=args.agCol, arCol=args.arCol, aiCol=args.aiCol, azCol=args.azCol,
        ayCol=args.ayCol, snError=args.snError, verbose=args.verbose,
        extCor=args.extCor)
