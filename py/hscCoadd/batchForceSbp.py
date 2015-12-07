#!/usr/bin/env python
# encoding: utf-8

import os
import glob
import logging
import warnings
import argparse

from astropy.io import fits

import coaddCutoutSbp as cSbp


def run(args):
    """
    Run coaddCutoutSbp in batch mode.

    Parameters:
    """
    if os.path.isfile(args.incat):
        data = fits.open(args.incat)[1].data
        id = (args.id)
        rerun = (args.rerun).strip()
        prefix = (args.prefix).strip()
        suffix = (args.suffix).strip()
        filter = (args.filter).strip().upper()

        """ Keep a log """
        logFile = (args.incat).replace('.fits', '_%s_forcesbp.log' % rerun)
        logging.basicConfig(filename=logFile)

        print "#########################################################"
        print "## Will deal with %d galaxies ! " % len(data)

        for galaxy in data:

            galID = str(galaxy[id]).strip()

            print "#########################################################\n"
            galPrefix = prefix + '_' + galID + '_' + filter + '_full'

            galRoot = os.path.join(galID, filter)
            if not os.path.isdir(galRoot):
                raise Exception('### Can not find the root folder ' +
                                ' for the galaxy data !')
            fitsList = glob.glob(os.path.join(galRoot, '*.fits'))

            galImg = galPrefix + '_img.fits'
            if (not os.path.isfile(os.path.join(galRoot, galImg)) and not
                    os.path.islink(os.path.join(galRoot, galImg))):
                raise Exception('### Can not find the cutout image of ' +
                                'the galaxy !')

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
                    raise Exception('### Can not find the final mask of ' +
                                    'the galaxy !')
            else:
                galMsk = None

            """
            Input Ellip Binary File
            """
            """ The reference filter """
            refFilter = (args.refFilter).strip().upper()
            galRefRoot = os.path.join(galID, refFilter, rerun)
            galRefPrefix = (prefix + '_' + galID + '_' + refFilter +
                            '_full_img_ellip_')
            """ The reference model """
            inEllipPrefix = os.path.join(galRefRoot, galRefPrefix)
            refModel = (args.refModel).strip()
            inEllipBin = inEllipPrefix + refModel + '.bin'
            print "###   INPUT ELLIP BIN : %s" % inEllipBin
            if not os.path.isfile(inEllipBin):
                raise Exception('### Can not fin the input ellip ' +
                                'bin : %s' % inEllipBin)

            """
            Suffix of the output file
            """
            ellipSuffix = rerun + '_' + suffix
            try:
                cSbp.coaddCutoutSbp(galPrefix, root=galRoot,
                                    verbose=args.verbose,
                                    psf=args.psf,
                                    inEllip=inEllipBin,
                                    bkgCor=args.bkgCor,
                                    zp=args.zp,
                                    step=4,
                                    galX0=args.galX0,
                                    galY0=args.galY0,
                                    galQ0=args.galQ0,
                                    galPA0=args.galPA0,
                                    galRe=args.galRe,
                                    checkCenter=args.noCheckCenter,
                                    updateIntens=args.updateIntens,
                                    pix=args.pix,
                                    plot=args.plot,
                                    redshift=args.redshift,
                                    olthresh=args.olthresh,
                                    fracBad=args.fracBad,
                                    lowClip=args.lowClip,
                                    uppClip=args.uppClip,
                                    nClip=args.nClip,
                                    intMode=args.intMode,
                                    minIt=args.minIt,
                                    maxIt=args.maxIt,
                                    maxTry=args.maxTry,
                                    outRatio=args.outRatio,
                                    exMask=galMsk,
                                    suffix=ellipSuffix,
                                    plMask=args.plmask)

                logging.warning('### The 1-D SBP is DONE for %s' % galPrefix)
            except Exception, errMsg:
                print str(errMsg)
                warnings.warn('### The 1-D SBP is failed for %s' % galPrefix)
                logging.warning('### The 1-D SBP is FAILED for %s' % galPrefix)
            print "#########################################################"

    else:
        raise Exception("### Can not find the input catalog: %s" % args.incat)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="Prefix of the galaxy image files")
    parser.add_argument("incat", help="The input catalog for cutout")
    parser.add_argument("filter", help="Filter to analysis")
    parser.add_argument('-i', '--id', dest='id',
                        help="Name of the column for galaxy ID",
                        default='ID')
    parser.add_argument('-mf', '--mFilter', dest='maskFilter',
                        help="Filter for Mask", default=None)
    parser.add_argument('-rf', '--rFilter', dest='refFilter',
                        help="Reference filter for Ellipse run",
                        default='HSC-I')
    parser.add_argument('-ref', '--rModel', dest='refModel',
                        help="Reference ellipse binary output",
                        default='default_3')
    parser.add_argument('-r', '--rerun', dest='rerun',
                        help="Name of the rerun", default='default')
    parser.add_argument("--suffix",
                        help="Suffix of the output file",
                        default='')
    """ Optional """
    parser.add_argument("--intMode", dest='intMode',
                        help="Method for integration",
                        default='median')
    parser.add_argument('--pix', dest='pix', help='Pixel Scale',
                        type=float, default=0.168)
    parser.add_argument('--step', dest='step', help='Step size',
                        type=float, default=0.10)
    parser.add_argument('--zp', dest='zp', help='Photometric zeropoint',
                        type=float, default=27.0)
    parser.add_argument('--redshift', dest='redshift',
                        help='Photometric zeropoint',
                        type=float, default=None)
    parser.add_argument('--olthresh', dest='olthresh',
                        help='Central locator threshold',
                        type=float, default=0.30)
    parser.add_argument('--uppClip', dest='uppClip',
                        help='Upper limit for clipping',
                        type=float, default=2.0)
    parser.add_argument('--lowClip', dest='lowClip',
                        help='Upper limit for clipping',
                        type=float, default=3.0)
    parser.add_argument('--nClip', dest='nClip',
                        help='Upper limit for clipping',
                        type=int, default=3)
    parser.add_argument('--fracBad', dest='fracBad',
                        help='Outer threshold',
                        type=float, default=0.5)
    parser.add_argument('--minIt', dest='minIt',
                        help='Minimum number of iterations',
                        type=int, default=10)
    parser.add_argument('--maxIt', dest='maxIt',
                        help='Maximum number of iterations',
                        type=int, default=100)
    parser.add_argument('--maxTry', dest='maxTry',
                        help='Maximum number of attempts of ellipse run',
                        type=int, default=3)
    parser.add_argument('--galX0', dest='galX0',
                        help='Center X0',
                        type=float, default=None)
    parser.add_argument('--galY0', dest='galY0', help='Center Y0',
                        type=float, default=None)
    parser.add_argument('--galQ0', dest='galQ0', help='Input Axis Ratio',
                        type=float, default=None)
    parser.add_argument('--galPA0', dest='galPA0', help='Input Position Angle',
                        type=float, default=None)
    parser.add_argument('--galRe', dest='galRe',
                        help='Input Effective Radius in pixel',
                        type=float, default=None)
    parser.add_argument('--outRatio', dest='outRatio',
                        help='Increase the boundary of SBP by this ratio',
                        type=float, default=1.2)
    parser.add_argument('--verbose', dest='verbose', action="store_true",
                        default=True)
    parser.add_argument('--psf', dest='psf', action="store_true",
                        help='Ellipse run on PSF', default=True)
    parser.add_argument('--plot', dest='plot', action="store_true",
                        help='Generate summary plot', default=True)
    parser.add_argument('--bkgCor', dest='bkgCor', action="store_true",
                        help='Background correction', default=False)
    parser.add_argument('--noCheckCenter', dest='noCheckCenter',
                        action="store_false",
                        help='Check if the center is off', default=True)
    parser.add_argument('--updateIntens', dest='updateIntens',
                        action="store_true",
                        default=True)
    parser.add_argument('--plmask', dest='plmask', action="store_true",
                        default=False)

    args = parser.parse_args()

    run(args)
