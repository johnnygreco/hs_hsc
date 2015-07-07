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
import galSBP

def readSbpInput(prefix, root=None):
    """
    doc
    """
    # Get the names of necessary input images
    imgFile = prefix + '_img.fits'
    mskFile = prefix + '_mskfin.fits'

    if root is not None:
        imgFile = os.path.join(root, imgFile)
        mskFile = os.path.join(root, mskFile)

    if not os.path.isfile(imgFile):
        raise Exception("### Can not find the input cutout image : %s !" % imgFile)
    if not os.path.isfile(mskFile):
        raise Exception("### Can not find the input mask image : %s !" % mskFile)

    # Image
    imgHdu = fits.open(imgFile)
    imgArr = imgHdu[0].data
    imgHead = imgHdu[0].header
    # Mask
    mskHdu = fits.open(mskFile)
    mskArr = mskHdu[0].data
    mskHead = mskHdu[0].header

    return imgFile, imgArr, imgHead, mskFile, mskArr, mskHead


def readPsfInput(prefix, root=None):
    """
    doc
    """
    psfFile = prefix + '_psf.fits'
    psfFile = os.path.join(root, psfFile)

    if root is not None:
        psfFile = os.path.join(root, psfFile)

    if not os.path.isfile(psfFile):
        raise Exception("### Can not find the input psf image : %s !" % psfFile)

    # PSF
    psfHdu = fits.open(psfFile)
    psfArr = psfHdu[0].data
    psfDimX, psfDimY = psfArr.shape

    return psfFile, psfDimX, psfDimY


def readInputSky(skyFile):
    """
    doc
    """

    if not os.path.isfile(skyFile):
        raise Exception("### Can not find the input sky summary : %s !" % skyFile)

    skySum = open(skyFile, 'r').readlines()
    skyMed = float(skySum[3].split(':')[1].strip())
    skyAvg = float(skySum[4].split(':')[1].strip())
    skyStd = float(skySum[5].split(':')[1].strip())
    skySkw = float(skySum[6].split(':')[1].strip())

    return skyMed, skyAvg, skyStd


def imgSameSize(img1, img2):
    """
    doc
    """
    dimX1, dimY1 = img1.shape
    dimX2, dimY2 = img2.shape
    if (dimX1 == dimX2) and (dimY1 = dimY2):
        return True
    else:
        return False


def prepareSbpInput(imgHead, mskHead):
    """
    doc
    """


def coaddCutoutSbp(prefix, root=None, verbose=True, psf=True, inEllip=None,
        ellStage=3, bkg=True, zp=None):
    """
    doc
    """

    """ 0. Organize Input Data """
    # Read in the input image, mask, psf, and their headers
    imgFile, imgArr, imgHead, mskFile, mskArr, mskHead = readSbpInput(prefix,
            root=root)
    if not imgSameSize(imgArr, mskArr):
        raise Exception("### The Image and Mask need to have EXACTLY same dimensions!")





