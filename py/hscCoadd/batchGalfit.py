#!/usr/bin/env python
# encoding: utf-8
"""Run GALFIT to fit 2-D models."""

import os
import gc
import glob
import logging
import warnings
import argparse

try:
    import psutil
    psutilOk = True
except Exception:
    psutilOk = False
from astropy.io import fits

# import coaddCutoutGalfitSimple as cGalfit
from coaddCutoutGalfitSimple import coaddCutoutGalfitSimple

COM = '#' * 100
SEP = '-' * 100
WAR = '!' * 100


def run(args):
    """
    Run coaddCutoutSbp in batch mode.

    Parameters
    """
    if psutilOk:
        proc = psutil.Process(os.getpid())
        gc.collect()
        mem0 = proc.memory_info().rss
        print SEP
        print "@@@ Initial: %i" % mem0
        print SEP
    else:
        gc.collect()

    if os.path.isfile(args.incat):
        """ Read the input catalog """
        data = fits.open(args.incat)[1].data

        id = (args.id)

        rerun = (args.rerun).strip()
        prefix = (args.prefix).strip()
        filter = (args.filter).strip().upper()

        """ Keep a log """
        logSuffix = '_%s_%s_galfit.log' % (filter, rerun)
        logFile = (args.incat).replace('.fits', logSuffix)
        logging.basicConfig(filename=logFile)

        print COM
        print "## Will deal with %d galaxies ! " % len(data)

        for index, galaxy in enumerate(data):

            print COM
            galID = str(galaxy[id]).strip()
            galPrefix = prefix + '_' + galID + '_' + filter + '_full'
            galRoot = os.path.join(galID, filter)
            print "## Deal with %s now : %i / %i" % (galID,
                                                     (index + 1),
                                                     len(data))
            print COM
            if not os.path.isdir(galRoot):
                logging.warning('### Can not find ' +
                                'ROOT folder for %s' % galRoot)
                continue
            fitsList = glob.glob(os.path.join(galRoot, '*.fits'))
            if len(fitsList) < 2:
                logging.warning('### MISSING Data in %s' % galRoot)
                continue

            galImg = galPrefix + '_img.fits'
            if (not os.path.isfile(os.path.join(galRoot, galImg)) and not
                    os.path.islink(os.path.join(galRoot, galImg))):
                logging.warning('### Can not find ' +
                                'CUTOUT IMAGE for %s' % galPrefix)
                continue

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
                """ Filter for mask """
                mskFilter = (args.maskFilter).strip().upper()
                """ Type of mask """
                mskType = (args.maskType).strip().upper()
                print "###  Use %s filter for mask \n" % mskFilter
                mskPrefix = prefix + '_' + galID + '_' + mskFilter + '_full'
                mskRoot = os.path.join(galID, mskFilter, rerun)
                galMsk = os.path.join(mskRoot, mskPrefix + '_' +
                                      mskType + '.fits')
                if not os.path.isfile(galMsk):
                    logging.warning('### Can not find the mask : %s' % galMsk)
                    continue
            else:
                galMsk = None

            print '\n' + SEP
            try:
                coaddCutoutGalfitSimple(galPrefix, root=galRoot,
                                        pix=args.pix,
                                        zp=args.zp,
                                        verbose=args.verbose,
                                        useBkg=args.useBkg,
                                        usePsf=args.usePsf,
                                        useSig=args.useSig,
                                        model=args.model,
                                        mag=args.mag,
                                        run1=args.run1,
                                        run2=args.run2,
                                        run3=args.run3,
                                        skyGrad=args.skyGrad,
                                        ser2Comp=args.ser2Comp,
                                        ser3Comp=args.ser3Comp,
                                        useF4=args.useF4,
                                        useF1=args.useF1,
                                        checkCenter=args.checkCenter,
                                        constrCen=args.constrCen,
                                        deleteAfter=args.deleteAfter,
                                        maskType=args.maskType,
                                        externalMask=galMsk,
                                        abspath=args.abspath,
                                        imax=args.imax)
                logging.info('### The Galfit Run is DONE for %s' % galPrefix)
                print SEP
            except Exception, errMsg:
                print str(errMsg)
                warnings.warn('### The Galfit Run is failed for %s' %
                              galPrefix)
                logging.warning('### The Galfit Run is FAILED for %s' %
                                galPrefix)
                print SEP + '\n'

            if psutilOk:
                mem1 = proc.memory_info().rss
                gc.collect()
                mem2 = proc.memory_info().rss
                print "@@@ Collect: %0.2f%%" % (100.0 * (mem2 - mem1) / mem0)
            else:
                gc.collect()

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
    parser.add_argument('-mf', '--maskFilter', dest='maskFilter',
                        help="Filter for Mask", default=None)
    parser.add_argument('-mt', '--maskType', dest='maskType',
                        help='Type of the mask use', default='mskfin')
    parser.add_argument('-r', '--rerun', dest='rerun',
                        help="Name of the rerun", default='default')
    parser.add_argument('-v', '--verbose', dest='verbose', action="store_true",
                        default=False)
    """ Optional """
    parser.add_argument('--model', dest='model',
                        help='Suffix of the model',
                        default=None)
    parser.add_argument('--pix', dest='pix', help='Pixel Scale',
                        type=float, default=0.168)
    parser.add_argument('--zp', dest='zp', help='Photometric zeropoint',
                        type=float, default=27.0)
    parser.add_argument('--imax', dest='imax',
                        help='Maximum number of iterations',
                        type=int, default=150)
    parser.add_argument('--mag', dest='mag', help='Total magnitude',
                        type=float, default=18.00)
    parser.add_argument('--noBkg', dest='useBkg', action="store_false",
                        default=True)
    parser.add_argument('--noPsf', dest='usePsf', action="store_false",
                        default=True)
    parser.add_argument('--noSig', dest='useSig', action="store_false",
                        default=True)
    parser.add_argument('--run1', dest='run1', action="store_true",
                        default=False)
    parser.add_argument('--run2', dest='run2', action="store_true",
                        default=False)
    parser.add_argument('--run3', dest='run3', action="store_true",
                        default=False)
    parser.add_argument('--ser2Comp', dest='ser2Comp', action="store_true",
                        default=False)
    parser.add_argument('--ser3Comp', dest='ser3Comp', action="store_true",
                        default=False)
    parser.add_argument('--skyGrad', dest='skyGrad', action="store_true",
                        default=False)
    parser.add_argument('--useF1', dest='useF1', action="store_true",
                        default=False)
    parser.add_argument('--useF4', dest='useF4', action="store_true",
                        default=False)
    parser.add_argument('--noConstrCen', dest='constrCen',
                        action="store_false", default=True)
    parser.add_argument('--noCheckCenter', dest='checkCenter',
                        action="store_false", default=True)
    parser.add_argument('--deleteAfter', dest='deleteAfter',
                        action="store_true", default=False)
    parser.add_argument('--abspath', dest='abspath',
                        action="store_true", default=False)

    args = parser.parse_args()

    run(args)
