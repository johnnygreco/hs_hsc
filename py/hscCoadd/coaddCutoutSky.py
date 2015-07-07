#!/usr/bin/env python
# encoding: utf-8

from __future__ import division

import os
import copy
import argparse
import numpy as np
import scipy

# Astropy
from astropy.io import fits
from astropy    import units as u
from astropy.stats import sigma_clip
# AstroML
from astroML.plotting import hist

# SEP
import sep

# Cubehelix color scheme
import cubehelix  # Cubehelix color scheme from https://github.com/jradavenport/cubehelix
# For high-contrast image
cmap1 = cubehelix.cmap(start=0.5, rot=-0.8, gamma=1.0,
                       minSat=1.2, maxSat=1.2,
                       minLight=0.0, maxLight=1.0)
cmap1.set_bad('k',1.)
# For Mask
cmap2 = cubehelix.cmap(start=2.0, rot=-1.0, gamma=2.5,
                       minSat=1.2, maxSat=1.2,
                       minLight=0.0, maxLight=1.0, reverse=True)
# For Sigma
cmap3 = cubehelix.cmap(start=0.5, rot=-0.8, gamma=1.2,
                       minSat=1.2, maxSat=1.2,
                       minLight=0.0, maxLight=1.0)

# Matplotlib related
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['figure.figsize'] = 12, 10
mpl.rcParams['xtick.major.size'] = 8.0
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.size'] = 4.0
mpl.rcParams['xtick.minor.width'] = 1.5
mpl.rcParams['ytick.major.size'] = 8.0
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.size'] = 4.0
mpl.rcParams['ytick.minor.width'] = 1.5
mpl.rc('axes', linewidth=2)
import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.patches import Ellipse

# Personal
import hscUtils as hUtil
import coaddCutoutPrepare as cdPrep


def readCutout(prefix, root=None):

    # Get the names of necessary input images
    imgFile = prefix + '_img.fits'
    mskFile = prefix + '_mskall.fits'

    if root is not None:
        imgFile = os.path.join(root, imgFile)
        mskFile = os.path.join(root, mskFile)
    if (not os.path.isfile(imgFile)) or (not os.path.isfile(mskFile)):
        print imgFile
        print mskFile
        raise Exception("### Can not find the Image or BadMask File!")
    else:
        imgHdu = fits.open(imgFile)
        imgArr = imgHdu[0].data
        # Header
        imgHead = imgHdu[0].header
        # All objects mask
        mskHdu = fits.open(mskFile)
        mskArr = mskHdu[0].data

    imgArrV = imgArr.view('float32')
    #imgArr = cdPrep.imgByteSwap(imgArr)
    #mskArr = cdPrep.imgByteSwap(mskArr)

    return imgArr, imgHead, mskArr


def showSkyHist(skypix, skypix2=None, skypix3=None,
                sbExpt=None, pngName='skyhist.png', skyAvg=None, skyStd=None,
                skyMed=None, skySkw=None):
    """
    Plot the distribution of sky pixels

    """
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    fig.subplots_adjust(hspace=0.1, wspace=0.1,
                        top=0.95, right=0.95)
    fontsize = 18
    ax.minorticks_on()

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    counts1, bins1, patches1 = hist(skypix, bins='knuth', ax=ax, alpha=0.4,
                                    color='cyan', histtype='stepfilled', normed=True)
    if skypix2 is not None:
        counts2, bins2, patches2 = hist(skypix2, bins='knuth', ax=ax, alpha=0.9,
                                        color='k', histtype='step', normed=True,
                                        linewidth=2)
    if skypix3 is not None:
        counts3, bins3, patches3 = hist(skypix3, bins='knuth', ax=ax, alpha=0.8,
                                        color='k', histtype='step', normed=True,
                                        linewidth=2, linestyle='dashed')
    # Horizontal line
    ax.axvline(0.0, linestyle='-', color='k', linewidth=1.5)

    # Basic properties of the sky pixels
    skyMin = np.nanmin(skypix)
    skyMax = np.nanmax(skypix)
    if skyAvg is None:
        skyAvg = np.nanmean(skypix)
    if skyStd is None:
        skyStd = np.nanstd(skypix)
    if skyMed is None:
        skyMed = np.nanmedian(skypix)
    if skySkw is None:
        skySkw = scipy.stats.skew(skypix)
    # Highligh the mode of sky pixel distribution
    ax.axvline(skyMed, linestyle='--', color='b', linewidth=1.5)

    ax.set_xlabel('Pixel Value', fontsize=20)
    ax.set_xlim(skyMed - 4.0 * skyStd, skyMed + 5.0 * skyStd)
    # Show a few information
    ax.text(0.7, 0.9, "Min : %8.4f" % skyMin, fontsize=21, transform=ax.transAxes)
    ax.text(0.7, 0.8, "Max : %8.4f" % skyMax, fontsize=21, transform=ax.transAxes)
    ax.text(0.7, 0.7, "Avg : %8.4f" % skyAvg, fontsize=21, transform=ax.transAxes)
    ax.text(0.7, 0.6, "Std : %8.4f" % skyStd, fontsize=21, transform=ax.transAxes)
    ax.text(0.7, 0.5, "Med : %8.4f" % skyMed, fontsize=21, transform=ax.transAxes)
    ax.text(0.7, 0.4, "Skew: %8.4f" % skySkw, fontsize=21, transform=ax.transAxes)
    if sbExpt is not None:
        ax.text(0.7, 0.3, "S.B : %8.5f" % sbExpt, fontsize=21, transform=ax.transAxes)

    fig.savefig(pngName)
    plt.close(fig)


def getIsophoteSky(imgArr, mskArr):
    """
    Estimating the background value within isophotes around the galaxy

    """
    # TODO
    a=1


def getRadBoxSky(imgArr, mskArr):
    """
    Estimating the background by thowing random boxes or apertures on the
    image

    This could also be used to estimate the detection limit of the image

    """


def getGlobalSky(imgArr, mskAll, skyClip=3, zp=27.0, pix=0.168,
                 rebin=4, prefix='coadd_sky', suffix=None,
                 verbose=True, visual=True):
    """
    Estimating the global sky background level by using the mean
    of a rebined image

    This could also be used to estimate the expect surface brightness
    limit of the image

    """
    # Estimate the global background level
    if verbose:
        print "### ESTIMATING THE GLOBAL BACKGROUND AND SURFACE BRIGHTNESS LIMIT"
    dimX, dimY = imgArr.shape
    # Pixel values of all pixels that are not masked out (before rebinned)
    pixels = imgArr[mskAll == 0].flatten()
    pixNoMsk = sigma_clip(pixels, skyClip, 3)

    # Rebin image
    dimBinX = int((dimX-1) / rebin)
    dimBinY = int((dimY-1) / rebin)
    print "###   REBIN IMAGE "
    imgBin = hUtil.congrid(imgArr, (dimBinX, dimBinY), method='nearest')
    print "###   REBIN MASK "
    mskBin = hUtil.congrid(mskAll, (dimBinX, dimBinY), method='neighbour')
    # Get all the pixels that are not masked out
    pixels = imgBin[mskBin == 0].flatten()
    pixNoMskBin = sigma_clip(pixels, skyClip, 3)
    numSkyPix = pixNoMskBin.shape[0]
    if verbose:
        print "### Global Background After Rebin the Image "
        print "###     N Pixels: %10d" % numSkyPix
    # Get the basic statistics of the global sky
    skyAvg, skyStd = np.nanmean(pixNoMskBin), np.nanstd(pixNoMskBin)
    skyMed = np.nanmedian(pixNoMskBin)
    skySkw = scipy.stats.skew(pixNoMskBin)
    sbExpt = cdPrep.getSbpValue(3.0 * skyStd, pix*rebin, pix*rebin, zp=zp)
    if verbose:
        print "###    Median Sky: %8.5f" % skyMed
        print "###      Mean Sky: %8.5f" % skyAvg
        print "###    StdDev Sky: %8.5f" % skyStd
        print "###      Skewness: %8.5f" % skySkw
        print "###    SB. Expect: %8.2f" % sbExpt
    if visual:
        skyPNG = prefix + '_' + suffix + 'skyhist.png'
        showSkyHist(pixNoMskBin, skypix2=pixNoMsk, sbExpt=sbExpt,
                pngName=skyPNG, skyAvg=skyAvg, skyMed=skyMed, skyStd=skyStd,
                skySkw=skySkw)

    # Save a txt file summary
    skyTxt = prefix + '_' + suffix + 'sky.dat'
    text_file = open(skyTxt, "w")
    text_file.write("IMAGE: %s \n" % prefix)
    text_file.write("REBIN: %3d \n" % rebin)
    text_file.write("NSKYPIX: %10d \n" % numSkyPix)
    text_file.write("SKYMED: %10.6f \n" % skyMed)
    text_file.write("SKYAVG: %10.6f \n" % skyAvg)
    text_file.write("SKYSTD: %10.6f \n" % skyStd)
    text_file.write("SKYSKW: %10.6f \n" % skySkw)
    text_file.write("SBEXPT: %10.6f \n" % sbExpt)
    text_file.close()

def coaddCutoutSky(prefix, root=None, verbose=True, skyClip=3.0,
                   pix=0.168, zp=27.0, rebin=6, visual=True):
    """
    doc
    """

    # 0. Get necessary information
    # Read the input cutout image
    imgArr, imgHead, mskArr = readCutout(prefix, root=root)
    if verbose:
        print "##########################################################################"
        print "### DEAL WITH IMAGE : %s" % (prefix + '_img.fits')
    # Necessary information
    if verbose:
        print "###    The pixel scale in X/Y directions " + \
                "are %7.4f / %7.4f arcsecs" % (pix, pix)
        print "###    The photometric zeropoint is %6.2f " % zp
        print "###    A %3d x %3d binning will be applied" % (rebin, rebin)
        print "###    A %4.1f sigma-clipping will be applied" % skyClip
    # Get rid of the NaN pixels, if there is any
    mskArr[np.isnan(imgArr)] = 1

    # 1. Global Background Estimation
    suffixGlob = 'rebin' + str(rebin).strip() + '_'
    getGlobalSky(imgArr, mskArr, skyClip=skyClip, zp=zp, pix=pix,
                 rebin=rebin, prefix=prefix, suffix=suffixGlob,
                 visual=visual, verbose=verbose)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="Prefix of the cutout image files")
    parser.add_argument('-r', '--root', dest='root', help='Path to the image files',
                        default=None)
    parser.add_argument('--skyclip', dest='skyClip', help='Sigma for pixel clipping',
                       type=float, default=3.0)
    parser.add_argument('--rebin', dest='rebin', help='Rebin the image by N x N pixels',
                       type=int, default=6)
    parser.add_argument('--pix', dest='pix', help='Pixel scale of the iamge',
                       type=float, default=0.168)
    parser.add_argument('--zp', dest='zp', help='Photometric zeropoint of the image',
                       type=float, default=27.0)
    parser.add_argument('--verbose', dest='verbose', action="store_true", default=False)
    parser.add_argument('--visual', dest='visual', action="store_true", default=False)

    args = parser.parse_args()

    coaddCutoutSky(args.prefix, root=args.root, pix=args.pix, zp=args.zp,
                   rebin=args.rebin, skyClip=args.skyClip, verbose=args.verbose,
                   visual=args.visual)
