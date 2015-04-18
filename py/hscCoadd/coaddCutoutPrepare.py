#!/usr/bin/env python

from __future__ import division

import os
import copy
import argparse
import numpy as np

from astropy.io import fits
from astropy    import units as u

import cubehelix  # Cubehelix color scheme from https://github.com/jradavenport/cubehelix

import matplotlib as mpl
import matplotlib.pyplot as plt
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


def zscale(img, contrast=0.25, samples=500):

    # Image scaling function form http://hsca.ipmu.jp/hscsphinx/scripts/psfMosaic.html
    ravel = img.ravel()
    if len(ravel) > samples:
        imsort = np.sort(np.random.choice(ravel, size=samples))
    else:
        imsort = np.sort(ravel)

    n = len(imsort)
    idx = np.arange(n)

    med = imsort[n/2]
    w = 0.25
    i_lo, i_hi = int((0.5-w)*n), int((0.5+w)*n)
    p = np.polyfit(idx[i_lo:i_hi], imsort[i_lo:i_hi], 1)
    slope, intercept = p

    z1 = med - (slope/contrast)*(n/2-n*w)
    z2 = med + (slope/contrast)*(n/2-n*w)

    return z1, z2


def readCutoutHeader(imgHead):

    # Get the pixel scale of the image
    try:
        pixScaleX = (np.abs(imgHead['CD1_1']) * 3600.0)
        pixScaleY = (np.abs(imgHead['CD2_2']) * 3600.0)
    except KeyError:
        pixScaleX = pixScaleY = 0.168
    else:
        pixScaleX = pixScaleY = 0.168
    print "The pixel scale in X/Y directions are %7.4f / %7.4f \
            arcsecs" % (pixScaleX, pixScaleY)

    # Get the image size
    imgSizeX = imgHead['NAXIS1']
    imgSizeY = imgHead['NAXIS2']
    print "The image size in X/Y directions are %d / %d \
            pixels" % (imgSizeX, imgSizeY)
    print "      %10.2f / %10.2f arcsecs" % (imgSizeX * pixScaleX, imgSizeY * pixScaleY)

    return pixScaleX, pixScaleY, imgSizeX, imgSizeY


def readCutoutImage(file):

    # Open the HDU list
    cutout = fits.open(file)
    # Image Data
    imgData = cutout[1].data
    # Mask Data
    mskData = cutout[2].data
    # Variance Data
    sigData = np.sqrt(cutout[3].data)

    # Get the header information of the image
    imgHead = cutout[1].header

    return imgData, mskData, sigData, imgHead


def showCutoutImg(imgData, mskData, sigData, prefix):

    # Output PNG file name
    outPNG = prefix + '_summary.png'

    fig = plt.figure(figsize=(14, 14))

    # TODO

    fig.savefig(outPNG)

def showPsfImg(prefix):

    # The name of the PSF image
    psfFile = prefix + '_psf.fits'

    if os.path.isfile(psfFile):
        psfImg = fits.open(psfFile)[0].data
    # TODO


def coaddCutoutPrepare(cutout, showImg=True, showPSF=True):

    # Get the prefix of the cutout image
    prefix = os.path.splitext(cutout)[0]

    # Read the input cutout image
    if os.path.isfile(cutout):
        imgData, mskData, sigData, imgHead = readCutoutImage(cutout)
        pixX, pixY, dimX, dimY = readCutoutHeader(imgHead)
    else:
        raise Exception("Can not find the input cutout image ! %s" % cutout )

    # If necessary, show a summary figure of the cutout image
    if showImg:
        showCutoutImg(imgData, mskData, sigData, prefix)

    # Show PSF image
    if showPSF:
        showPsfImg(prefix)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("cutout", help="Name of the cutout image file")
    parser.add_argument('-i', '--info', dest='info',
                        help='Information to show on the image',
                        default=None)
    args = parser.parse_args()

    coaddCutoutPrepare(args.cutout, info=args.info)
