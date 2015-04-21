#!/usr/bin/env python

import os
import numpy
import argparse
import warnings

from astropy.io import fits

import coaddImageCutout as cdCutout
import coaddColourImage as cdColor


def decideCutoutSize(z):

    """
    Decide the typical cutout size for certain redshift
    """

    if (z <= 0.15):
        return 800
    if (z > 0.15) and (z < 0.25):
        return 500
    if (z > 0.25) and (z < 0.35):
        return 400
    if (z > 0.35) and (z < 0.45):
        return 350
    if (z > 0.45):
        return 300


def parseInputCatalog(list, sizeDefault=300, idField='id',
                     raField='ra', decField='dec', sizeField='cutout_size',
                     zField=None, zCutoutSize=False):

    # Read in the catalog
    hduList = fits.open(list)
    cat = hduList[1].data

    # Try to get the ID, Ra, Dec
    try:
        id = cat.field(idField)
    except KeyError:
        raise Exception('Can not find the ID field')

    try:
        ra = cat.field(raField)
    except KeyError:
        raise Exception('Can not find the RA field')

    try:
        dec = cat.field(decField)
    except KeyError:
        raise Exception('Can not find the DEC field')

    if zField is not None:
        try:
            redshift = cat.field(zField)
        except KeyError:
            raise Exception('Can not find the REDSHIFT field')
    else:
        redshift = None

    if zCutoutSize and (redshift is not None):
        size = map(lambda x: decideCutoutSize(x), redshift)
        size = numpy.asarray(size)
    else:
        try:
            size = cat.field(sizeField)
        except KeyError:
            warnings.warn("### No field name for cutout size is provided !")
            nObjs = len(id)
            size = numpy.empty(nObjs)
            size.fill(sizeDefault)

    return id, ra, dec, size, redshift


def coaddBatchCutout(root, inCat, size=100, filter='HSC-I', prefix='coadd_cutout',
                     idField='id', raField='ra', decField='dec', colorFilters='gri',
                     sizeField='cutout_size', zCutoutSize=False, zField=None,
                     verbose=True, noColor=False, onlyColor=False):
    """
    Givin an input catalog with RA, DEC information, generate HSC
    coadd cutout images.

    Also have the option to generate (or just generate) a 3-band
    color image
    """

    if not os.path.isdir(root):
        raise Exception("Wrong root directory for data! %s" % root)

    if not os.path.exists(inCat):
        raise Exception("Can not find the input catalog! %s" % list)
    else:
        id, ra, dec, size, z = parseInputCatalog(inCat, sizeDefault=size,
                                                 idField=idField, raField=raField,
                                                 decField=decField,
                                                 sizeField=sizeField, zField=zField,
                                                 zCutoutSize=zCutoutSize)

    if not onlyColor:
        logFile = prefix + '_match_status.lis'
        logMatch = open(logFile, 'w')

    nObjs = len(id)
    if verbose:
        print "### Will try to get cutout image for %d objects" % nObjs

    for i in range(nObjs):

        if verbose:
            print "### %d -- ID: %d ; RA: %10.5f DEC %10.5f ; Size: %d" % (i+1, id[i], ra[i], dec[i], size[i])

        # New prefix
        newPrefix = prefix + '_' + str(id[i]).strip()

        # Cutout Image
        if not onlyColor:
            coaddFound, noData, partialCut = cdCutout.coaddImageCutout(root, ra[i], dec[i],
                                                                       size[i], saveMsk=True,
                                                                       filt=filter,
                                                                       prefix=newPrefix)
            if coaddFound:
                if not noData:
                    if not partialCut:
                        matchStatus = 'Full'
                    else:
                        matchStatus = 'Part'
                else:
                    matchStatus = 'NoData'
            else:
                matchStatus = 'Outside'
            logMatch.write(str(id[i]) + '   ' + matchStatus + '\n')

        # Color Image
        if not noColor:
            if (zField is not None) and (z is not None):
                info = "z= %5.3f" % z[i]
            else:
                info = None
            if (matchStatus is 'Full') or (matchStatus is 'Part'):
                name = str(id[i])
                cdColor.coaddColourImage(root, ra[i], dec[i], size[i], filt=colorFilters,
                                         prefix=newPrefix, name=name,
                                         info=info )

    if not onlyColor:
        logMatch.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("list", help="The input catalog for cutout")
    parser.add_argument("-s", '--size', dest='size', type=int,
                        help="Half size of the cutout box", default=200)
    parser.add_argument('-f', '--filter', dest='filt', help="Filter",
                       default='HSC-I')
    parser.add_argument('-p', '--prefix', dest='prefix',
                        help='Prefix of the output file',
                        default='hsc_coadd_cutout')
    args = parser.parse_args()

    coaddBatchCutout(args.root, args.list, size=args.size,
                     filter=args.filt, prefix=args.prefix)
