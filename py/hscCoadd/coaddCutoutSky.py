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


def showSkyHist(skypix, skypix2=None, pngName='skyhist.png'):
    """
    Plot the distribution of sky pixels

    """
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    fig.subplots_adjust(hspace=0.1, wspace=0.1,
                        top=0.95, right=0.95)
    fontsize = 18
    ax.minorticks_on()

    ax.set_xlim(-0.7, 1.0)

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    counts1, bins2, patches3 = hist(skypix, bins='knuth', ax=ax, alpha=0.4,
                                    color='cyan', histtype='stepfilled', normed=True)
    counts1, bins2, patches3 = hist(skypix2, bins='knuth', ax=ax, alpha=0.9,
                                    color='k', histtype='step', normed=True, linewidth=2)
    ax.axvline(0.0, linestyle='-', color='k', linewidth=1.5)
    ax.axvline(np.nanmedian(skypix), linestyle='--', color='b', linewidth=1.5)

    ax.set_xlabel('Pixel Value', fontsize=20)
    # TODO: Adjust axes range ; Add sky information

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


def getGlobalSky(imgArr, mskArr, skyClip=3, zp=27.0, pix=0.168,
                 rebin=4, prefix='coadd_sky', suffix=None):
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
    # Get the basic statistics of the global sky
    meanSky1, stdSky1 = np.nanmean(pixNoMsk), np.nanstd(pixNoMsk)
    medSky1 = np.nanmedian(pixNoMsk)
    sbExp1 = getSbpValue(stdSky1, pixX, pixY, zp=zp)
    if verbose:
        print "### Global Background Before Rebin the Image "
        print "###    Median Sky: %8.5f" % medSky1
        print "###      Mean Sky: %8.5f" % meanSky1
        print "###    StdDev Sky: %8.5f" % stdSky1
        print "###    SB. Expect: %8.2f" % sbExp1
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
    #pixNoMskBin = pixels
    # Get the basic statistics of the global sky
    meanSky2, stdSky2 = np.nanmean(pixNoMskBin), np.nanstd(pixNoMskBin)
    medSky2 = np.nanmedian(pixNoMskBin)
    sbExp2 = cdPrep.getSbpValue(stdSky2, pixX*rebin, pixY*rebin, zp=photZP)
    if verbose:
        print "### Global Background After Rebin the Image "
        print "###    Median Sky: %8.5f" % medSky2
        print "###      Mean Sky: %8.5f" % meanSky2
        print "###    StdDev Sky: %8.5f" % stdSky2
        print "###    SB. Expect: %8.2f" % sbExp2
    if visual:
        skyPNG = prefix + '_' + suffix + 'skyhist.png'
        showSkyHist(pixNoMskBin, skypix2=pixNoMsk, pngName=skyPNG)

