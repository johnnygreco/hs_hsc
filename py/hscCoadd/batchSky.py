#!/usr/bin/env python
# encoding: utf-8

import os
import glob
import logging
import argparse
import warnings

from astropy.io import fits

import coaddCutoutSky as ccs

COM = '#' * 100
SEP = '-' * 100
WAR = '!' * 100


def run(args):
    """
    Run coaddCutoutPrepare in batch mode.

    Parameters:
    """
    if os.path.isfile(args.incat):

        data = fits.open(args.incat)[1].data

        id = (args.id)
        rerun = (args.rerun).strip()
        prefix = (args.prefix).strip()
        filter = (args.filter).strip().upper()

        """ Keep a log """
        logSuffix = '_%s_%s_sky.log' % (filter, rerun)
        logFile = (args.incat).replace('.fits', logSuffix)
        logging.basicConfig(filename=logFile)

        print COM
        print "## Will deal with %d galaxies ! " % len(data)
        print COM

        for galaxy in data:

            galID = str(galaxy[id]).strip()

            print SEP
            galPrefix = prefix + '_' + galID + '_' + filter + '_full'

            galRoot = os.path.join(galID, filter)
            if not os.path.isdir(galRoot):
                print WAR
                raise Exception('### Can not find the root folder' +
                                ' for the galaxy data !')
            fitsList = glob.glob(os.path.join(galRoot, '*.fits'))
            if len(fitsList) <= 3:
                print WAR
                raise Exception("### Missing data under %s" % galRoot)

            """
            Set up a rerun
            """
            galRoot = os.path.join(galRoot, rerun.strip())
            if not os.path.isdir(galRoot):
                os.makedirs(galRoot)
            """ Link the necessary files to the rerun folder """
            for fitsFile in fitsList:
                seg = fitsFile.split('/')
                link = os.path.join(galRoot, seg[-1])
                if (not os.path.islink(link)) and (not os.path.isfile(link)):
                    os.symlink(fitsFile, link)

            """
            External mask
            """
            if args.maskFilter is not None:
                mskFilter = (args.maskFilter).strip().upper()
                print "###  Use %s filter for mask \n" % mskFilter
                mskPrefix = prefix + '_' + galID + '_' + mskFilter + '_full'
                mskRoot = os.path.join(galID, mskFilter, rerun)
                galMsk = os.path.join(mskRoot, mskPrefix + '_mskfin.fits')
                if not os.path.isfile(galMsk):
                    raise Exception('### Can not find the final mask ' +
                                    'of the galaxy !')
            else:
                galMsk = None

            try:
                ccs.coaddCutoutSky(galPrefix, root=galRoot,
                                   pix=args.pix,
                                   zp=args.zp,
                                   rebin=args.rebin,
                                   skyClip=args.skyClip,
                                   verbose=args.verbose,
                                   visual=args.visual,
                                   exMask=galMsk)
            except Exception, errMsg:
                print WAR
                print str(errMsg)
                warnings.warn('### The sky estimate is failed ' +
                              'for %s' % args.prefix)
                logging.warning('### The sky estimate is failed ' +
                                'for %s' % args.prefix)
            print COM
    else:
        raise Exception("### Can not find the input catalog: %s" % args.incat)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="Prefix of the galaxy image files")
    parser.add_argument("incat", help="The input catalog for cutout")
    parser.add_argument('-i', '--id', dest='id',
                        help="Name of the column for galaxy ID",
                        default='ID')
    parser.add_argument('-f', '--filter', dest='filter', help="Filter",
                        default='HSC-I')
    parser.add_argument('-m', '--mask', dest='mask', help="Filter for Mask",
                        default=None)
    parser.add_argument('-r', '--rerun', dest='rerun',
                        help="Name of the rerun", default='default')
    parser.add_argument('-mf', '--mFilter', dest='maskFilter',
                        help="Filter for Mask", default=None)
    """ Optional """
    parser.add_argument('--skyclip', dest='skyClip',
                        help='Sigma for pixel clipping',
                        type=float, default=3.0)
    parser.add_argument('--rebin', dest='rebin',
                        help='Rebin the image by N x N pixels',
                        type=int, default=6)
    parser.add_argument('--pix', dest='pix',
                        help='Pixel scale of the iamge',
                        type=float, default=0.168)
    parser.add_argument('--zp', dest='zp',
                        help='Photometric zeropoint of the image',
                        type=float, default=27.0)
    parser.add_argument('--verbose', dest='verbose',
                        action="store_true", default=True)
    parser.add_argument('--visual', dest='visual',
                        action="store_true", default=True)

    args = parser.parse_args()

    run(args)
