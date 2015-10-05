#!/usr/bin/env python
# encoding: utf-8

import os
import argparse
import coaddCutoutPrepare as ccp

from astropy.io import fits

def run(args):

    if os.path.isfile(args.incat):

        data = fits.open(args.incat)[1].data

        id     = (args.id)
        rerun  = (args.rerun).strip()
        prefix = (args.prefix).strip()
        filter = (args.filter).strip().upper()

        print "## Will deal with %d galaxies ! " % len(data)

        for galaxy in data:

            galID = str(galaxy[id]).strip()

            print "################################################################\n"
            galPrefix = prefix + '_' + galID + '_' + filter + '_full'
            print "################################################################\n"

            galRoot   = os.path.join(galID, filter)
            galImg    = galPrefix + '_img.fits'

            if not os.path.isdir(galRoot):
                raise Exception('### Can not find the root folder for the galaxy data !')
            if not os.path.isfile(os.path.join(galRoot, galImg)):
                raise Exception('### Can not find the cutout image of the galaxy !')

            ccp.coaddCutoutPrepare(galPrefix, root=galRoot, rerun=rerun)

            #coaddCutoutPrepare(args.prefix, root=args.root,
                               #bSizeH=args.bSizeH, bSizeC=args.bSizeC,
                               #thrH=args.thrH, thrC=args.thrC,
                               #growH=args.growH, growW=args.growW, growC=args.growC,
                               #kernel=args.kernel, central=args.central,
                               #mask=args.mask, useSigArr=args.useSigArr,
                               #noBkgC=args.noBkgC, noBkgH=args.noBkgH,
                               #minDetH=args.minDetH, minDetC=args.minDetC,
                               #debThrH=args.debThrH, debThrC=args.debThrC,
                               #debConH=args.debConH, debConC=args.debConC,
                               #combBad=args.combBad, combDet=args.combDet,
                               #rerun=args.rerun)

    else:
        raise Exception("### Can not find the input catalog: %s" % args.incat)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="Prefix of the galaxy image files")
    parser.add_argument("incat", help="The input catalog for cutout")
    parser.add_argument('-i', '--id', dest='id', help="Name of the column for galaxy ID",
                       default='ID')
    parser.add_argument('-f', '--filter', dest='filter', help="Filter",
                       default='HSC-I')
    parser.add_argument('-r', '--rerun', dest='rerun',
                        help="Name of the rerun", default='default')
    args = parser.parse_args()

    run(args)

