#!/usr/bin/env python
"""Generate HSC full cutout in batch mode."""

import os
import numpy
import argparse
import warnings

from astropy.io import fits
# HSC Pipeline
import lsst.daf.persistence as dafPersist

import coaddImageCutout as cdCutout
import coaddColourImage as cdColor

COM = '#' * 100
SEP = '-' * 100
WAR = '!' * 100
HSC_FILTERS = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y']


def decideCutoutSize(z, safe=False):
    """
    Decide the typical cutout size for certain redshift.

    Parameters:
        safe  : True will make the cutout smaller
    """
    if (z <= 0.15):
        if safe:
            return 1000
        else:
            return 1200
    if (z > 0.15) and (z < 0.25):
        if safe:
            return 650
        else:
            return 750
    if (z > 0.25) and (z < 0.35):
        if safe:
            return 550
        else:
            return 600
    if (z > 0.35) and (z < 0.45):
        if safe:
            return 350
        else:
            return 400
    if (z > 0.45) and (z < 0.65):
        if safe:
            return 300
        else:
            return 350
    if (z > 0.65):
        if safe:
            return 200
        else:
            return 250


def parseInputCatalog(list, sizeDefault=300, idField='id',
                      raField='ra', decField='dec', sizeField='cutout_size',
                      zField=None, zCutoutSize=False, infoField1=None,
                      infoField2=None, safe=False):
    """
    Parse the input catalog.

    Parameters:
    """
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
            print WAR
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
            print WAR
            print "### Can not find field: %s in the catalog !" % infoField2
            info3 = None
        else:
            if isinstance(info3[0], (float, int, numpy.number)):
                info3 = map(lambda x: "{:10.3f}".format(x).strip(), info3)
    else:
        info3 = None

    if zCutoutSize and (redshift is not None):
        size = map(lambda x: decideCutoutSize(x, safe=safe), redshift)
        size = numpy.asarray(size)
    else:
        try:
            size = cat.field(sizeField)
        except KeyError:
            # warnings.warn("### No field name for cutout size is provided !")
            nObjs = len(id)
            size = numpy.empty(nObjs)
            size.fill(sizeDefault)

    return id, ra, dec, size, redshift, info2, info3


def coaddBatchCutout(root, inCat, size=100, filter='HSC-I',
                     prefix='coadd_cutout', idField='id',
                     raField='ra', decField='dec', colorFilters='gri',
                     sizeField='cutout_size', zCutoutSize=False,
                     zField=None, verbose=True, noColor=False,
                     onlyColor=False, infoField1=None, infoField2=None,
                     clean=False, min=-0.0, max=0.72, Q=15, stitch=False,
                     noName=False):
    """
    Generate HSC coadd cutout images in batch mode.

    Also have the option to generate (or just generate) a 3-band
    color image
    """
    if not os.path.isdir(root):
        raise Exception("Wrong root directory for data! %s" % root)

    if os.path.exists(inCat):
        result = parseInputCatalog(inCat, sizeDefault=size, idField=idField,
                                   raField=raField, decField=decField,
                                   zField=zField, zCutoutSize=zCutoutSize,
                                   infoField1=infoField1,
                                   infoField2=infoField2)
        id, ra, dec, size, z, extr1, extr2 = result
    else:
        raise Exception("### Can not find the input catalog: %s" % inCat)

    if not onlyColor:
        logFile = prefix + '_match_status_' + filter.strip() + '.lis'
        logMatch = open(logFile, 'w')

    nObjs = len(id)
    if verbose:
        print SEP
        print "### Will try to get cutout image for %d objects" % nObjs
        print SEP

    for i in range(nObjs):

        if verbose:
            print "### %d -- ID: %s ; " % (i+1, str(id[i])) + \
                  "RA: %10.5f DEC %10.5f ; Size: %d" % (ra[i], dec[i], size[i])

        # New prefix
        newPrefix = prefix + '_' + str(id[i]).strip()

        # Cutout Image
        if not onlyColor:
            tempOut = cdCutout.coaddImageCutout(root, ra[i], dec[i], size[i],
                                                saveMsk=True, filt=filter,
                                                prefix=newPrefix)
            coaddFound, noData, partialCut = tempOut
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
            if noName:
                name = None
            else:
                name = str(id[i])
            if stitch:
                if not clean:
                    cdColor.coaddColourImageFull(root, ra[i], dec[i], size[i],
                                                 filt=colorFilters,
                                                 scaleBar=10,
                                                 prefix=newPrefix, name=name,
                                                 info1=info1, info2=info2,
                                                 info3=info3, min=min, max=max,
                                                 Q=Q)
                else:
                    cdColor.coaddColourImageFull(root, ra[i], dec[i], size[i],
                                                 filt=colorFilters,
                                                 scaleBar=10,
                                                 prefix=newPrefix, name=None,
                                                 info1=None, info2=None,
                                                 info3=None, min=min, max=max,
                                                 Q=Q)
            else:
                cdColor.coaddColourImage(root, ra[i], dec[i], size[i],
                                         filt=colorFilters,
                                         prefix=newPrefix, name=name,
                                         info1=info1, info2=info2, info3=info3,
                                         min=min, max=max, Q=Q)
        elif (matchStatus is 'Full') or (matchStatus is 'Part'):
            if noName:
                name = None
            else:
                name = str(id[i])
            if stitch:
                cdColor.coaddColourImageFull(root, ra[i], dec[i], size[i],
                                             filt=colorFilters,
                                             prefix=newPrefix, name=name,
                                             info1=info1, info2=info2,
                                             info3=info3, min=min, max=max,
                                             Q=Q)
            else:
                cdColor.coaddColourImage(root, ra[i], dec[i], size[i],
                                         filt=colorFilters,
                                         prefix=newPrefix, name=name,
                                         info1=info1, info2=info2, info3=info3,
                                         min=min, max=max, Q=Q)
    if not onlyColor:
        logMatch.close()


def coaddBatchCutFull(root, inCat, size=100, filter='HSC-I',
                      prefix='coadd_cutout', idField='id',
                      raField='ra', decField='dec', colorFilters='gri',
                      sizeField='cutout_size', zCutoutSize=False, zField=None,
                      verbose=True, noColor=False, onlyColor=False,
                      infoField1=None, infoField2=None, clean=False,
                      min=-0.0, max=0.72, Q=15, safe=False, saveSrc=False,
                      makeDir=False, noName=False,
                      imgOnly=False, allFilters=False):
    """
    Generate HSC coadd cutout images.

    Also have the option to generate (or just generate) a 3-band
    color image
    """
    if not os.path.isdir(root):
        raise Exception("Wrong root directory for data! %s" % root)
    else:
        try:
            butler = dafPersist.Butler(root)
            if verbose:
                "### Load in the Butler "
        except Exception:
            warnings.warn("### Can not load the Butler!")

    if os.path.exists(inCat):
        if verbose:
            print COM
            print "              PARSE THE INPUT CATALOG                 "
            print COM
        tempOut = parseInputCatalog(inCat, sizeDefault=size, idField=idField,
                                    raField=raField, decField=decField,
                                    zField=zField, zCutoutSize=zCutoutSize,
                                    infoField1=infoField1,
                                    infoField2=infoField2, safe=safe)
        id, ra, dec, size, z, extr1, extr2 = tempOut
    else:
        raise Exception("### Can not find the input catalog: %s" % inCat)

    if not onlyColor:
        logFile = prefix + '_match_status.lis'
        logMatch = open(logFile, 'w')

    nObjs = len(id)
    if verbose:
        print SEP
        print "### Will try to get cutout image for %d objects" % nObjs
        print SEP

    for i in range(nObjs):
        if verbose:
            print "### %d -- ID: %s ; " % (i+1, str(id[i])) + \
                  "RA: %10.5f DEC %10.5f ; Size: %d" % (ra[i], dec[i], size[i])
        # New prefix
        newPrefix = prefix + '_' + str(id[i]).strip()
        # Cutout Image
        if not onlyColor:
            if verbose:
                print "### Make the Cutout Fits Files!  "
            if not allFilters:
                filterUse = filter.strip()
                if makeDir:
                    dirLoc = (str(id[i]).strip() + '/' +
                              str(filterUse).strip() + '/')
                    if not os.path.exists(dirLoc):
                        os.makedirs(dirLoc)
                    filterPrefix = dirLoc + newPrefix
                else:
                    filterPrefix = newPrefix

                if saveSrc:
                    tempOut = cdCutout.coaddImageCutFull(root, ra[i], dec[i],
                                                         size[i], savePsf=True,
                                                         saveSrc=True,
                                                         visual=True,
                                                         filt=filterUse,
                                                         prefix=filterPrefix,
                                                         butler=butler,
                                                         imgOnly=imgOnly)
                    found, full, npatch = tempOut
                else:
                    tempOut = cdCutout.coaddImageCutFull(root, ra[i], dec[i],
                                                         size[i], savePsf=True,
                                                         saveSrc=False,
                                                         visual=True,
                                                         filt=filterUse,
                                                         prefix=filterPrefix,
                                                         butler=butler,
                                                         imgOnly=imgOnly)
                    found, full, npatch = tempOut
                if found:
                    matchStatus = 'Found'
                    if full:
                        full = 'Full'
                    else:
                        full = 'Part'
                else:
                    matchStatus = 'NoData'
                    full = 'None'

                logMatch.write(str(id[i]) + '    ' + filterUse +
                               '   ' + matchStatus +
                               '   ' + full + '   ' +
                               str(npatch) + '\n')
            else:
                for filterUse in HSC_FILTERS:
                    print "## Working on %s now" % filterUse

                    if makeDir:
                        dirLoc = (str(id[i]).strip() + '/' +
                                  str(filterUse).strip() + '/')
                        if not os.path.exists(dirLoc):
                            os.makedirs(dirLoc)
                        filterPrefix = dirLoc + newPrefix
                    else:
                        filterPrefix = newPrefix

                    if saveSrc:
                        tempOut = cdCutout.coaddImageCutFull(root,
                                                             ra[i], dec[i],
                                                             size[i],
                                                             savePsf=True,
                                                             saveSrc=True,
                                                             visual=True,
                                                             filt=filterUse,
                                                             prefix=filterPrefix,
                                                             butler=butler,
                                                             imgOnly=imgOnly)
                        found, full, npatch = tempOut
                    else:
                        tempOut = cdCutout.coaddImageCutFull(root,
                                                             ra[i], dec[i],
                                                             size[i],
                                                             savePsf=True,
                                                             saveSrc=False,
                                                             visual=True,
                                                             filt=filterUse,
                                                             prefix=filterPrefix,
                                                             butler=butler,
                                                             imgOnly=imgOnly)
                        found, full, npatch = tempOut
                    if found:
                        matchStatus = 'Found'
                        if full:
                            full = 'Full'
                        else:
                            full = 'Part'
                    else:
                        matchStatus = 'NoData'
                        full = 'None'
                    logMatch.write(str(id[i]) + '    ' + filterUse +
                                   '   ' + matchStatus +
                                   '   ' + full + '   ' +
                                   str(npatch) + '\n')

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
            if noName:
                name = None
            else:
                name = str(id[i])
            if verbose:
                print "### Generate Color Image !"
            if clean:
                cdColor.coaddColourImageFull(root, ra[i], dec[i], size[i],
                                             filt=colorFilters,
                                             prefix=newPrefix, name=None,
                                             info1=None, info2=None,
                                             info3=None,
                                             min=min, max=max, Q=Q,
                                             butler=butler)
            else:
                cdColor.coaddColourImageFull(root, ra[i], dec[i], size[i],
                                             filt=colorFilters,
                                             prefix=newPrefix, name=name,
                                             info1=info1, info2=info2,
                                             info3=info3,
                                             min=min, max=max, Q=Q,
                                             butler=butler)
        elif (matchStatus is 'Found' and not noColor):
            if noName:
                name = None
            else:
                name = str(id[i])
            if verbose:
                print "### Generate Color Image !"
            cdColor.coaddColourImageFull(root, ra[i], dec[i], size[i],
                                         filt=colorFilters,
                                         prefix=newPrefix, name=name,
                                         info1=info1, info2=info2, info3=info3,
                                         min=min, max=max, Q=Q, butler=butler)
    if not onlyColor:
        logMatch.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("incat", help="The input catalog for cutout")
    parser.add_argument("-s", '--size', dest='size', type=int,
                        help="Half size of the cutout box", default=200)
    parser.add_argument('-f', '--filter', dest='filt', help="Filter",
                        default='HSC-I')
    parser.add_argument('-p', '--prefix', dest='prefix',
                        help='Prefix of the output file',
                        default='hsc_coadd_cutout')
    parser.add_argument('-id', '--id', dest='idField',
                        help="Column name for ID", default='index')
    parser.add_argument('-ra', '--ra', dest='raField',
                        help="Column name for RA",
                        default='ra_hsc')
    parser.add_argument('-dec', '--dec', dest='decField',
                        help="Column name for DEC",
                        default='dec_hsc')
    parser.add_argument('-z', '--redshift', dest='zField',
                        help="Column name for z",
                        default='z_use')
    parser.add_argument('-cf', '--color-filters', dest='colorFilters',
                        help="Choice of filters for color images",
                        default='riz')
    parser.add_argument('-sf', '--size-field', dest='sizeField',
                        help="Column name for cutout size",
                        default='cutout_size')
    parser.add_argument('-info1', '--infoField1', dest='infoField1',
                        help="Column name for first extra information",
                        default=None)
    parser.add_argument('-info2', '--infoField2', dest='infoField2',
                        help="Column name for second extra information",
                        default=None)
    parser.add_argument('-af', '--allFilters', action="store_true",
                        dest='allFilters', default=False)
    parser.add_argument('-img', '--imgOnly', action="store_true",
                        dest='imgOnly', default=False)
    parser.add_argument('-zc', '--zCutoutSize', action="store_true",
                        dest='zCutoutSize', default=True)
    parser.add_argument('-nc', '--noColor', action="store_true",
                        dest='noColor', default=True)
    parser.add_argument('-oc', '--onlyColor', action="store_true",
                        dest='onlyColor', default=False)
    parser.add_argument('-safe', '--safe', action="store_true", dest='safe',
                        default=False)
    parser.add_argument('-clean', '--clean', action="store_true", dest='clean',
                        default=False)
    parser.add_argument('-nn', '--noName', action="store_true",
                        dest='noName', default=False)
    parser.add_argument('-v', '--verbose', action="store_true", dest='verbose',
                        default=False)
    parser.add_argument('-src', '--src', action="store_true", dest='saveSrc',
                        default=False)
    parser.add_argument('-makeDir', '--makeDir', action="store_true",
                        dest='makeDir', default=False)
    args = parser.parse_args()

    coaddBatchCutFull(args.root, args.incat, size=args.size,
                      filter=args.filt, prefix=args.prefix,
                      idField=args.idField, raField=args.raField,
                      decField=args.decField, sizeField=args.sizeField,
                      colorFilters=args.colorFilters, zField=args.zField,
                      zCutoutSize=args.zCutoutSize, noColor=args.noColor,
                      onlyColor=args.onlyColor, infoField1=args.infoField1,
                      infoField2=args.infoField2, safe=args.safe,
                      verbose=args.verbose, clean=args.clean,
                      saveSrc=args.saveSrc,
                      makeDir=args.makeDir, noName=args.noName,
                      imgOnly=args.imgOnly, allFilters=args.allFilters)
