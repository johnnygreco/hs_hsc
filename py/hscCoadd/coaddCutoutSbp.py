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
    if (dimX1 == dimX2) and (dimY1 == dimY2):
        return True
    else:
        return False


def coaddCutoutSbp(prefix, root=None, verbose=True, psf=True, inEllip=None,
        ellStage=3, bkg=True, zp=27.0):
    """
    doc
    """

    """ 0. Organize Input Data """
    # Read in the input image, mask, psf, and their headers
    imgFile, imgArr, imgHead, mskFile, mskArr, mskHead = readSbpInput(prefix,
            root=root)
    if not imgSameSize(imgArr, mskArr):
        raise Exception("### The Image and Mask need to have EXACTLY same dimensions!")

    """ 1. Prepare the Input for SBP """
    galX, galY = mskHead['GAL_CENX'], mskHead['GAL_CENY']
    galQ, galPA = mskHead['GAL_Q'], mskHead['GAL_PA']
    galR50 = mskHead['GAL_R50']

    maxR = mskHead['NAXIS1'] if mskHead['NAXIS1'] >= mskHead['NAXIS2'] else mskHead['NAXIS2']
    maxR *= np.sqrt(2.0)

    if mskHead['MSK_R20'] != 1:
        print "### The central region is masked out"
    else:
        if inEllip is None:
            """ # Start with Stage 1 """
            maxSma1 = maxR * 0.8
            ellOut1 = galSBP.galSBP(imgFile, mskFile, galX=galX, galY=galY,
                                    maxSma=maxSma1, iniSma=int(galR50),
                                    galQ=galQ, galPA=galPA, stage=1, zpPhoto=zp)
            galX0 = ellOut1['avg_x0'][0]
            galY0 = ellOut1['avg_y0'][0]

            """ # Start with Stage 2 """
            maxSma2 = maxR
            ellOut2 = galSBP.galSBP(imgFile, mskFile, galX=galX0, galY=galY0,
                                    maxSma=maxSma2, iniSma=int(galR50),
                                    galQ=galQ, galPA=galPA, stage=2, zpPhoto=zp)
            galQ0  = ellOut1['avg_q'][0]
            galPA0 = ellOut1['avg_pa'][0]

            """ # Start with Stage 3 """
            maxSma3 = maxR
            ellOut3 = galSBP.galSBP(imgFile, mskFile, galX=galX0, galY=galY0,
                                    maxSma=maxSma3, iniSma=int(galR50),
                                    galQ=galQ0, galPA=galPA0, stage=3, zpPhoto=zp)
        else:
            """ # Run Ellipse in Forced Photometry Mode """
            ellOut4 = galSBP.galSBP(imgFile, mskFile, galX=galX, galY=galY,
                                    inEllip=inEllip, maxSma=maxSma1, iniSma=int(galR50),
                                    galQ=galQ, galPA=galPA, stage=1, zpPhoto=zp)





