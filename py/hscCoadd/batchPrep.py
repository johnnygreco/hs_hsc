#!/usr/bin/env python
# encoding: utf-8

import os
import logging
import argparse
import warnings

from astropy.io import fits

import coaddCutoutPrepare as ccp


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
        logFile = (args.incat).replace('.fits', '_%s_prep.log' % rerun)
        logging.basicConfig(filename=logFile)

        print "## Will deal with %d galaxies ! " % len(data)

        for galaxy in data:

            galID = str(galaxy[id]).strip()

            print "########################################################\n"
            galPrefix = prefix + '_' + galID + '_' + filter + '_full'
            print "########################################################\n"

            galRoot = os.path.join(galID, filter)
            galImg = galPrefix + '_img.fits'

            if not os.path.isdir(galRoot):
                raise Exception('### Can not find the root folder for the \
                        galaxy data !')
            if not os.path.isfile(os.path.join(galRoot, galImg)):
                raise Exception('### Can not find the cutout image of the \
                        galaxy !')

            try:
                if rerun == 'default':
                    ccp.coaddCutoutPrepare(galPrefix, root=galRoot,
                                           rerun=rerun,
                                           bSizeH=args.bSizeH,
                                           bSizeC=args.bSizeC,
                                           thrH=args.thrH,
                                           thrC=args.thrC,
                                           growH=args.growH,
                                           growW=args.growW,
                                           growC=args.growC,
                                           kernel=args.kernel,
                                           central=args.central,
                                           maskMethod=args.mask,
                                           growMethod=args.grow,
                                           useSigArr=args.useSigArr,
                                           noBkgC=args.noBkgC,
                                           noBkgH=args.noBkgH,
                                           minDetH=args.minDetH,
                                           minDetC=args.minDetC,
                                           debThrH=args.debThrH,
                                           debThrC=args.debThrC,
                                           debConH=args.debConH,
                                           debConC=args.debConC,
                                           combBad=args.combBad,
                                           combDet=args.combDet)
                elif rerun == 'new':
                    ccp.coaddCutoutPrepare(galPrefix, root=galRoot,
                                           rerun=rerun,
                                           bSizeH=args.bSizeH,
                                           bSizeC=args.bSizeC,
                                           thrH=args.thrH,
                                           thrC=args.thrC,
                                           growH=args.growH,
                                           growW=args.growW,
                                           growC=args.growC,
                                           kernel=args.kernel,
                                           central=args.central,
                                           maskMethod=args.mask,
                                           growMethod=args.grow,
                                           useSigArr=args.useSigArr,
                                           noBkgC=args.noBkgC,
                                           noBkgH=args.noBkgH,
                                           minDetH=args.minDetH,
                                           minDetC=args.minDetC,
                                           debThrH=args.debThrH,
                                           debThrC=args.debThrC,
                                           debConH=args.debConH,
                                           debConC=args.debConC,
                                           combBad=args.combBad,
                                           combDet=args.combDet)
                else:
                    ccp.coaddCutoutPrepare(galPrefix, root=galRoot,
                                           rerun=rerun,
                                           bSizeH=args.bSizeH,
                                           bSizeC=args.bSizeC,
                                           thrH=args.thrH,
                                           thrC=args.thrC,
                                           growH=args.growH,
                                           growW=args.growW,
                                           growC=args.growC,
                                           kernel=args.kernel,
                                           central=args.central,
                                           maskMethod=args.mask,
                                           growMethod=args.grow,
                                           useSigArr=args.useSigArr,
                                           noBkgC=args.noBkgC,
                                           noBkgH=args.noBkgH,
                                           minDetH=args.minDetH,
                                           minDetC=args.minDetC,
                                           debThrH=args.debThrH,
                                           debThrC=args.debThrC,
                                           debConH=args.debConH,
                                           debConC=args.debConC,
                                           combBad=args.combBad,
                                           combDet=args.combDet)
            except Exception:
                warnings.warn('### The cutout preparation is failed \
                        for %s' % galPrefix)
                logging.warning('### The cutout preparation is failed \
                        for %s' % galPrefix)
    else:
        raise Exception("### Can not find the input catalog: %s" % args.incat)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="Prefix of the galaxy image files")
    parser.add_argument("incat", help="The input catalog for cutout")
    parser.add_argument('-i', '--id', dest='id',
                        help="Name of the column for galaxy ID", default='ID')
    parser.add_argument('-f', '--filter', dest='filter', help="Filter",
                        default='HSC-I')
    parser.add_argument('-r', '--rerun', dest='rerun',
                        help="Name of the rerun", default='default')
    """ Optional """
    parser.add_argument('-k', dest='kernel',
                        help='SExtractor detection kernel',
                        type=int, default=4, choices=range(1, 7))
    parser.add_argument('-c', dest='central',
                        help='Method to clean the central region',
                        type=int, default=1, choices=range(1, 3))
    parser.add_argument('-m', dest='mask',
                        help='Method to grow the All object mask',
                        type=int, default=1, choices=range(1, 3))
    parser.add_argument('-g', dest='grow',
                        help='Method to grow the Final object mask',
                        type=int, default=1, choices=range(1, 2))
    parser.add_argument('--bkgH', dest='bSizeH',
                        help='Background size for the Hot Run',
                        type=int, default=10)
    parser.add_argument('--bkgC', dest='bSizeC',
                        help='Background size for the Cold Run',
                        type=int, default=80)
    parser.add_argument('--thrH', dest='thrH',
                        help='Detection threshold for the Hot Run',
                        type=float, default=2.5)
    parser.add_argument('--thrC', dest='thrC',
                        help='Detection threshold for the Cold Run',
                        type=float, default=1.2)
    parser.add_argument('--growC', dest='growC',
                        help='Ratio of Growth for the Cold Objects',
                        type=float, default=4.0)
    parser.add_argument('--growW', dest='growW',
                        help='Ratio of Growth for the Warm Objects',
                        type=float, default=3.0)
    parser.add_argument('--growH', dest='growH',
                        help='Ratio of Growth for the Hot Objects',
                        type=float, default=1.5)
    parser.add_argument('--minDetC', dest='minDetC',
                        help='Minimum pixels for Cold Detections',
                        type=float, default=8.0)
    parser.add_argument('--minDetH', dest='minDetH',
                        help='Minimum pixels for Hot Detections',
                        type=float, default=5.0)
    parser.add_argument('--debThrC', dest='debThrC',
                        help='Deblending threshold for the Cold Run',
                        type=float, default=32.0)
    parser.add_argument('--debThrH', dest='debThrH',
                        help='Deblending threshold for the Hot Run',
                        type=float, default=16.0)
    parser.add_argument('--debConC', dest='debConC',
                        help='Deblending continuum level for the Cold Run',
                        type=float, default=0.001)
    parser.add_argument('--debConH', dest='debConH',
                        help='Deblending continuum level for the Hot Run',
                        type=float, default=0.0001)
    parser.add_argument('--noBkgC', dest='noBkgC',
                        action="store_true", default=False)
    parser.add_argument('--noBkgH', dest='noBkgH',
                        action="store_true", default=False)
    parser.add_argument('--useSigArr', dest='useSigArr', action="store_true",
                        default=False)
    parser.add_argument('--combBad', dest='combBad', action="store_true",
                        default=True)
    parser.add_argument('--combDet', dest='combDet', action="store_true",
                        default=True)
    args = parser.parse_args()

    run(args)
