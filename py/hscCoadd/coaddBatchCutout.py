#!/usr/bin/env python

import os
import numpy
import argparse
import warnings

from astropy.io import fits

import coaddImageCutout as cdCutout
import coaddColourImage as cdColor

# HSC Pipeline
import lsst.daf.persistence  as dafPersist


def decideCutoutSize(z, safe=False):

    """
    Decide the typical cutout size for certain redshift
    """

    if (z <= 0.15):
        if safe:
            return 800
        else:
            return 1000
    if (z > 0.15) and (z < 0.25):
        if safe:
            return 500
        else:
            return 650
    if (z > 0.25) and (z < 0.35):
        if safe:
            return 400
        else:
            return 500
    if (z > 0.35) and (z < 0.45):
        if safe:
            return 350
        else:
            return 400
    if (z > 0.45):
        if safe:
            return 300
        else:
            return 350


def parseInputCatalog(list, sizeDefault=300, idField='id',
                     raField='ra', decField='dec', sizeField='cutout_size',
                     zField=None, zCutoutSize=False, infoField1=None,
                     infoField2=None):

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

    # Extra information 1
    if infoField1 is not None:
        try:
            info2 = cat.field(infoField1)
        except KeyError:
            print "### Can not find field: %s in the catalog !" % infoField1
            info2 = None
        else:
            if isinstance(info2[0], (float, int, numpy.number)):
                info2 = map(lambda x: "{:10.3f}".format(x).strip(), info2)
    else:
        info2 = None

    # Extra information 2
    if infoField2 is not None:
        try:
            info3 = cat.field(infoField2)
        except KeyError:
            print "### Can not find field: %s in the catalog !" % infoField2
        else:
            if isinstance(info3[0], (float, int, numpy.number)):
                info3 = map(lambda x: "{:10.3f}".format(x).strip(), info3)
    else:
        info3 = None

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

    return id, ra, dec, size, redshift, info2, info3


def coaddBatchCutout(root, inCat, size=100, filter='HSC-I', prefix='coadd_cutout',
                     idField='id', raField='ra', decField='dec', colorFilters='gri',
                     sizeField='cutout_size', zCutoutSize=False, zField=None,
                     verbose=True, noColor=False, onlyColor=False,
                     infoField1=None, infoField2=None,
                     min=-0.0, max=0.72, Q=15, stitch=False):
    """
    Givin an input catalog with RA, DEC information, generate HSC
    coadd cutout images.

    Also have the option to generate (or just generate) a 3-band
    color image
    """

    if not os.path.isdir(root):
        raise Exception("Wrong root directory for data! %s" % root)

    if os.path.exists(inCat):
        id, ra, dec, size, z, extr1, extr2  = parseInputCatalog(inCat, sizeDefault=size,
                                                 idField=idField, raField=raField,
                                                 decField=decField, zField=zField,
                                                 zCutoutSize=zCutoutSize,
                                                 infoField1=infoField1,
                                                 infoField2=infoField2)
    else:
        raise Exception("### Can not find the input catalog: %s" % inCat)


    if not onlyColor:
        logFile = prefix + '_match_status.lis'
        logMatch = open(logFile, 'w')

    nObjs = len(id)
    if verbose:
        print "### Will try to get cutout image for %d objects" % nObjs

    for i in range(nObjs):

        if verbose:
            print "### %d -- ID: %s ; " % (i+1, str(id[1])) + \
                    "RA: %10.5f DEC %10.5f ; Size: %d" % (ra[i], dec[i], size[i])

        # New prefix
        newPrefix = prefix + '_' + str(id[i]).strip()

        # Cutout Image
        if not onlyColor:
            coaddFound, noData, partialCut = cdCutout.coaddImageCutout(root, ra[i], dec[i],
                                                                       size[i],
                                                                       saveMsk=True,
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
        # Whether put redshift on the image
        if (zField is not None) and (z is not None):
            info1 = "z=%5.3f" % z[i]
        else:
            info1 = None
        # Extra information
        if (infoField1 is not None) and (extr1 is not None):
            info2 = str(extr1[i]).strip()
        else:
            info2 = None
        if (infoField2 is not None) and (extr2 is not None):
            info3 = str(extr2[i]).strip()
        else:
            info3 = None

        if onlyColor:
            name = str(id[i])
            if stitch:
                cdColor.coaddColourImageFull(root, ra[i], dec[i], size[i],
                                            filt=colorFilters,
                                            prefix=newPrefix, name=name,
                                            info1=info1, info2=info2, info3=info3,
                                            min=min, max=max, Q=Q)
            else:
                cdColor.coaddColourImage(root, ra[i], dec[i], size[i],
                                         filt=colorFilters,
                                         prefix=newPrefix, name=name,
                                         info1=info1, info2=info2, info3=info3,
                                         min=min, max=max, Q=Q)
        elif (matchStatus is 'Full') or (matchStatus is 'Part'):
            name = str(id[i])
            if stitch:
                cdColor.coaddColourImageFull(root, ra[i], dec[i], size[i],
                                            filt=colorFilters,
                                            prefix=newPrefix, name=name,
                                            info1=info1, info2=info2, info3=info3,
                                            min=min, max=max, Q=Q)
            else:
                cdColor.coaddColourImage(root, ra[i], dec[i], size[i],
                                         filt=colorFilters,
                                         prefix=newPrefix, name=name,
                                         info1=info1, info2=info2, info3=info3,
                                         min=min, max=max, Q=Q)
    if not onlyColor:
        logMatch.close()


def coaddBatchCutFull(root, inCat, size=100, filter='HSC-I', prefix='coadd_cutout',
                      idField='id', raField='ra', decField='dec', colorFilters='gri',
                      sizeField='cutout_size', zCutoutSize=False, zField=None,
                      verbose=True, noColor=False, onlyColor=False,
                      infoField1=None, infoField2=None,
                      min=-0.0, max=0.72, Q=15):
    """
    Givin an input catalog with RA, DEC information, generate HSC
    coadd cutout images.

    Also have the option to generate (or just generate) a 3-band
    color image
    """

    if not os.path.isdir(root):
        raise Exception("Wrong root directory for data! %s" % root)
    else:
        try:
            butler = dafPersist.Butler(root)
        except Exception:
            warnings.warn("### Get not load the Butler!")

    if os.path.exists(inCat):
        id, ra, dec, size, z, extr1, extr2  = parseInputCatalog(inCat, sizeDefault=size,
                                                 idField=idField, raField=raField,
                                                 decField=decField, zField=zField,
                                                 zCutoutSize=zCutoutSize,
                                                 infoField1=infoField1,
                                                 infoField2=infoField2)
    else:
        raise Exception("### Can not find the input catalog: %s" % inCat)


    if not onlyColor:
        logFile = prefix + '_match_status.lis'
        logMatch = open(logFile, 'w')

    nObjs = len(id)
    if verbose:
        print "### Will try to get cutout image for %d objects" % nObjs

    for i in range(nObjs):

        if verbose:
            print "### %d -- ID: %s ; " % (i+1, str(id[1])) + \
                    "RA: %10.5f DEC %10.5f ; Size: %d" % (ra[i], dec[i], size[i])

        # New prefix
        newPrefix = prefix + '_' + str(id[i]).strip()

        # Cutout Image
        if not onlyColor:
            coaddFound = cdCutout.coaddImageCutFull(root, ra[i], dec[i], size[i],
                                                     savePsf=True, saveSrc=True,
                                                     visual=True, filt=filter,
                                                     prefix=newPrefix, butler=butler)
            if coaddFound:
                matchStatus = 'Found'
            else:
                matchStatus = 'NoData'

            logMatch.write(str(id[i]) + '   ' + matchStatus + '\n')

        # Color Image
        # Whether put redshift on the image
        if (zField is not None) and (z is not None):
            info1 = "z=%5.3f" % z[i]
        else:
            info1 = None
        # Extra information
        if (infoField1 is not None) and (extr1 is not None):
            info2 = str(extr1[i]).strip()
        else:
            info2 = None
        if (infoField2 is not None) and (extr2 is not None):
            info3 = str(extr2[i]).strip()
        else:
            info3 = None

        if onlyColor:
            name = str(id[i])
            cdColor.coaddColourImageFull(root, ra[i], dec[i], size[i],
                                         filt=colorFilters,
                                         prefix=newPrefix, name=name,
                                         info1=info1, info2=info2, info3=info3,
                                         min=min, max=max, Q=Q)

        elif (matchStatus is 'Found'):
            name = str(id[i])
            cdColor.coaddColourImageFull(root, ra[i], dec[i], size[i],
                                         filt=colorFilters,
                                         prefix=newPrefix, name=name,
                                         info1=info1, info2=info2, info3=info3,
                                         min=min, max=max, Q=Q)

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
