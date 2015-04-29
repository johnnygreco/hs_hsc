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
# For Variance
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


def mergerObjCat(cat1, cat2, tol=1.5):
    """
    Merge two objects detections from SEP
    """

def getConvKernel(kernel):
    """
    Convolution kernel for the SEP detections
    """
    if kernel is 1:
        # Tophat_3.0_3x3
        convKer = np.asarray([[0.560000, 0.980000, 0.560000],
                              [0.980000, 1.000000, 0.980000],
                              [0.560000, 0.980000, 0.560000]])
    else:
        raise Exception("### More options will be available in the future")

    return convKer


def getEll2Plot(objects, radius=None):
    """
    Generate the ellipse shape for each object to plot
    """
    x  = objects['x'].copy()
    y  = objects['y'].copy()
    pa = objects['theta'].copy() # in unit of radian

    if radius is not None:
        a = radius.copy()
        b = radius.copy()*(objects['b'].copy()/objects['a'].copy())
    else:
        a  = objects['a'].copy()
        b  = objects['b'].copy()

    ells = [Ellipse(xy=np.array([x[i], y[i]]),
                    width=np.array(2.0 * b[i]),
                    height=np.array(2.0 * a[i]),
                    angle=np.array(pa[i]*180.0/np.pi + 90.0))
            for i in range(x.shape[0])]

    return ells


def getSbpValue(flux, pixX, pixY, zp=None):
    """
    Convert flux into surface brightness value

    TODO: Right now only support log-magnitude,
    In the future, should also support asinh-magnitude
    See:
    http://www.astro.washington.edu/users/ajc/ssg_page_working/elsst/opsim.shtml?lightcurve_mags
    """
    sbp = -2.5 * np.log10(flux/(pixX * pixY))
    if zp is not None:
        sbp += zp
    return sbp


def readCutoutHeader(imgHead):

    # Get the pixel scale of the image
    try:
        pixScaleX = pixScaleY = imgHead['PIXEL']
    except:
        print "### Pixel scale keyword is not available in the header"
        print "### Default 0.168 arcsec/pixel value is adopted"
        pixScaleX = pixScaleY = 0.168
    # Get the image size
    imgSizeX = imgHead['NAXIS1']
    imgSizeY = imgHead['NAXIS2']
    # Get the photometric zeropoint
    try:
        photZP = imgHead['PHOTZP']
    except:
        print "### PHOZP keyword is not available in the header"
        print "### Default value of 27.0 is adopted!"
        photZP = 27.0
    # Total exptime
    try:
        expTot = imgHead['TOTEXPT']
    except:
        print "### TOTEXPT keyword is not available in the header"
        print "### Use 1.0 sec instead"
        expTot = 1.0

    return pixScaleX, pixScaleY, imgSizeX, imgSizeY, photoZP, expTot


def imgByteSwap(data):
    """
    Byte Swap before sending image to SEP
    """
    dataCopy = copy.deepcopy(data)
    return dataCopy.byteswap(True).newbyteorder()


def sepGetBkg(img, mask=None, bkgSize=None, bkgFilter=None):
    """
    Wrapper of SEP.Background function
    """
    if bkgSize is None:
        dimX, dimY = img.shape
        bkgX = imt(dimX / 15)
        bkgY = imt(dimY / 15)
    else:
        bkgX = bkgY = int(bkgSize)
    if bkgFilter is None:
        bkgFilter = 4

    bkg = sep.Background(img, mask, bkgX, bkgY, bkgFilter, bkgFilter)
    # Subtract the Background off
    bkg.subfrom(img)

    return bkg, img


def readCutoutImage(prefix, root=None):

    # Get the names of necessary input images
    imgFile = prefix + '_img.fits'
    mskFile = prefix + '_bad.fits'
    detFile = prefix + '_det.fits'
    if root is not None:
        imgFile = os.path.join(root, imgFile)
        mskFile = os.path.join(root, mskFile)
        detFile = os.path.join(root, detFile)
    if (not os.path.isfile(imgFile)) or (not os.path.isfile(mskFile)):
        raise Exception("### Can not find the Image or BadMask File!")
    else:
        imgHdu = fits.open(imgFile)
        imgArr = imgHdu[0].data
        # Header
        imgHead = imgHdu[0].header
        # Bad mask
        mskHdu = fits.open(mskFile)
        mskArr = mskHdu[0].data
    # Optional detection plane
    if not os.path.isfile(detFile):
        print "### Can not find the coadd DetectionPlane file!"
        detArr = None
    else:
        detHdu = fits.open(detFile)
        detArr = detHdu[0].data

    return imgArr, imgHead, mskArr, detArr


def coaddCutoutPrepare(prefix, root=None, srcCat=None, verbose=True,
                       bSizeH=8, bSizeC=80, thrH=3.5, thrH=1.5,
                       growC=7.0, growW=3.5, growH=4.0,
                       galX=None, galY=None, galR1=None, galR2=None, galR3=None,
                       galQ=None, galPA=None):


    """
    The structure of the cutout has been changed.  Now the cutout procedure
    will generate separated files for Image, Bad Mask, Detection Plane, and
    Variance (also Sigma) images.  Souce catalogs can also be made available.

    Right now, this new format is only available for the coaddImageCutFull()
    function; coaddImageCutout() will be modified later to also adopt this
    format
    """

    # Read the input cutout image
    imgArr, imgHead, mskArr, detArr = readCutoutImage(prefix, root=root)
    if verbose:
        print "### Deal with image: %s" % (prefix + '_img.fits')
    if detArr is None:
        detFound = False
    else:
        detFound = True
    pixX, pixY, dimX, dimY, photZP, expTot = readCutoutHeader(imgHead)
    if verbose:
        print "### The pixel scale in X/Y directions " + \
                "are %7.4f / %7.4f arcsecs" % (pixX, pixY)
        print "### The image size in X/Y directions " + \
                "are %d / %d pixels" % (dimX, dimY)
        print "      %10.2f / %10.2f arcsecs" % (dimX * pixX, dimY * pixY)
        print "### The photometric zeropoint is %6.2f " % photoZP
        print "### The total exposure time for the center is %6.1f secs" % expTot

    """
    Construct "background" images with different size and filters using SEP,
    and subtract these background off before extract objects

    The SEP detections will be run in two-modes:
        Cold: relative global background; low-detection threshold
        Hot:  very local background; median-detection threshold
    """
    # Cold Run
    fSizeC = int(bSizeC / 4)
    imgC = imgByteSwap(imgArr)
    bkgC, imgSubC = sepGetBkg(imgC, bkgSize=bSizeC, bkgFilter=fSizeC)
    if verbose:
        print "### Cold Background -- Avg: %9.5f " % bkgC.globalback + \
                "  RMS: %9.5f " % bkgC.globalrms
    # Hot Run
    fSizeH = int(bSizeH / 2)
    imgH = imgByteSwap(imgArr)
    bkgH, imgSubH = sepGetBkg(imgH, bkgSize=bSizeH, bkgFilter=fSizeH)
    if verbose:
        print "### Hot Background  -- Avg: %9.5f " % bkgH.globalback + \
                "  RMS: %9.5f " % bkgH.globalrms

    """
    Use SEP to extract information of detected objects
    """
    """ Cold Run """
    # Parameters for cold run
    thrC *= bkgC.globalrms  # Cold run detection threshold
    minDetC = 5             # Minimun adjacent pixels for detection for cold run
    debThrC = 20            # Deblending threshold for cold run (32)
    debConC = 0.002         # Minimum deblending contrast ratio for cold run (0.005)
    convKerC = getConvKernel(1)
    # Cold Run
    objC = sep.extract(imgSubC, thrC, minarea=minDetC, conv=convKerC,
                       deblend_nthresh=debThrC, deblend_cont=debConC)
    if verbose:
        print "### %d objects have been detected in the Cold Run" % objC['x'].shape[0]

    """ Hot Run """
    # Parameters for hot run
    thrH *= bkgH.globalrms  # Hot run detection threshold
    minDetH = 5             # Minimun adjacent pixels for detection for hot run
    debThrH = 20            # Deblending threshold for hot run (32)
    debConH = 0.002         # Minimum deblending contrast ratio for hot run (0.005)
    convKerH = getConvKernel(1)
    # Hot Run
    objH = sep.extract(imgSubH, thrH, minarea=minDetH, conv=convKerH,
                       deblend_nthresh=debThrH, deblend_cont=debConH)
    if verbose:
        print "### %d objects have been detected in the Hot Run" % objH['x'].shape[0]





if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="Prefix of the cutout image files")
    parser.add_argument('-r', '--root', dest='root', help='Path to the image files',
                        default=None)
    args = parser.parse_args()

    coaddCutoutPrepare(args.prefix, root=args.root)
