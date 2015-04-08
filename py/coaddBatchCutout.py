#!/usr/bin/env python

import os
import numpy
import argparse
from astropy.io import fits
from coaddImageCutout import coaddImageCutout
from coaddColourImage import coaddColourImage

def parseInputCatalog(list, sizeDefault=100, idField='id',
                     raField='ra', decField='dec', sizeField='cutout_size'):

    # Read in the catalog
    hduList = fits.open(list)
    cat = hduList[1].data

    # Try to get the ID, Ra, Dec
    # TODO: provide the field name for ID, Ra, Dec, Size
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

    try:
        size = cat.field(sizeField)
    except KeyError:
        nObjs = len(id)
        size = numpy.empty(nObjs)
        size.fill(sizeDefault)

    return id, ra, dec, size


def coaddBatchCutout(root, list, size=100, filter='HSC-I', prefix='coadd_cutout',
                     idField='id', raField='ra', decField='dec',
                     sizeField='cutout_size'):

    if not os.path.isdir(root):
        raise Exception("Wrong root directory for data! %s" % root)

    if not os.path.exists(list):
        raise Exception("Can not find the input catalog! %s" % list)
    else:
        id, ra, dec, size = parseInputCatalog(list, sizeDefault=size,
                                              idField=idField, raField=raField,
                                              decField=decField,
                                              sizeField=sizeField)

    nObjs = len(id)

    for i in range(nObjs):
        # New prefix
        newPrefix = prefix + '_' + str(id[i]).strip()
        print newPrefix, ra[i], dec[i], size[i]
        # Cutout Image
        coaddImageCutout(root, ra[i], dec[i], size[i], saveMsk=True,
                         filt=filter, prefix=newPrefix)
        # Color Image
        coaddColourImage(root, ra[i], dec[i], size[i], filt='gri',
                        prefix=newPrefix)


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
