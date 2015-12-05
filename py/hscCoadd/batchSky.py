#!/usr/bin/env python
# encoding: utf-8

import os
import logging
import argparse
import warnings

from astropy.io import fits

import coaddCutoutSky as ccs

def run(args):

    if os.path.isfile(args.incat):

        data = fits.open(args.incat)[1].data

        id     = (args.id)
        rerun  = (args.rerun).strip()
        prefix = (args.prefix).strip()
        filter = (args.filter).strip().upper()

        """ Keep a log """
        logFile = (args.incat).replace('.fits', '_%s_sky.log' % rerun)
        logging.basicConfig(filename=logFile)

        print "## Will deal with %d galaxies ! " % len(data)

        for galaxy in data:

            galID = str(galaxy[id]).strip()

            print "################################################################\n"
            galPrefix = prefix + '_' + galID + '_' + filter + '_full'
            print "################################################################\n"

            galRoot   = os.path.join(galID, filter, rerun)
            galMsk    = galPrefix + '_mskfin.fits'

            if not os.path.isdir(galRoot):
                raise Exception('### Can not find the root folder for the galaxy data !')
            if not os.path.isfile(os.path.join(galRoot, galMsk)):
                raise Exception('### Can not find the final mask of the galaxy !')

            try:
                ccs.coaddCutoutSky(args.prefix, root=args.root,
                                   pix=args.pix,
                                   zp=args.zp,
                                   rebin=args.rebin,
                                   skyClip=args.skyClip,
                                   verbose=args.verbose,
                                   visual=args.visual)
            except Exception:
                warnings.warn('### The sky estimate is failed for %s') % args.prefix
                logging.warning('### The sky estimate is failed for %s') % args.prefix

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
    """ Optional """
    parser.add_argument('--skyclip', dest='skyClip', help='Sigma for pixel clipping',
                        type=float, default=3.0)
    parser.add_argument('--rebin', dest='rebin', help='Rebin the image by N x N pixels',
                        type=int, default=6)
    parser.add_argument('--pix', dest='pix', help='Pixel scale of the iamge',
                        type=float, default=0.168)
    parser.add_argument('--zp', dest='zp', help='Photometric zeropoint of the image',
                        type=float, default=27.0)
    parser.add_argument('--verbose', dest='verbose',
                        action="store_true", default=True)
    parser.add_argument('--visual', dest='visual',
                        action="store_true", default=True)

    args = parser.parse_args()

    run(args)

