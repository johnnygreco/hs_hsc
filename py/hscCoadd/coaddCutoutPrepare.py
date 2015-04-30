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


def saveSEPObjects(objs, output='sep_objects.fits'):
    """
    Save the properties of objects that are extracted by SEP into
    a file

    TODO: To be finished
    Two options:
        1. Numpy binary file: http://docs.scipy.org/doc/numpy/reference/generated/numpy.save.html
        2. Hickle/HDF5: https://github.com/telegraphic/hickle
    """
    nObjs = len(objs)


def adaptiveMask(objC, a=None, b=None, c=None):
    """
    Scale the size of mask for different objects according to their
    flux and major axis radii

    XXX This is still very empirical!
    We adopt the form of: Ratio = a*((log(A)-thrA) <= 0 ? 0) + b*logRho + c
    logRho = log10(Flux/(A*B));  and thrA=log10(seeing_fwhm/(2.0*pixscale))
    """
    # Corresponding to the pixel size of a seeing disk
    # Use as a threshold to select "point-ish" objects
    thrA = np.log10(seeing / (2.0 * pix))
    # "Flux density" like properties
    logRho = np.log10(objC['cflux']/(objC['a'] * objC['b']))
    # Logrithmic difference between the major axis radii and
    # the threshold size
    logA = (np.log10(np.sqrt(objC['a'])) - thrA )
    # Ignore all the "point-ish" objects
    logA[(logA < 0) | (logRho < 1.5)] = 0
    # Empirical function to get the mask size ratio
    rAdRatio = logA * a + logRho * b + c
    if verbose:
        print "### rAdRatio "
        print np.nanmin(rAdRatio)
        print np.nanmax(rAdRatio)
        print np.median(rAdRatio)

    return rAdRatio


def mergerObjCat(objCold, objHot, tol=1.0):
    """
    Merge the object lists from cold and hot run

    tol = 1.0  : Tolerance of central difference in unit of pixel
    """
    # Make a copy of the objects array
    objC = copy.deepcopy(objCold)
    objH = copy.deepcopy(objHot)
    # The central coordinates of each objects
    objCX = objC['xpeak']
    objCY = objC['ypeak']
    objHX = objH['xpeak']
    objHY = objH['ypeak']
    # Get the minimum separation between each Hot run object and all the
    # Cold run objects
    minDistH = np.asarray(map(lambda x, y: np.min(np.sqrt((x-objCX)**2.0
                                                         + (y-objCY)**2.0)),
                             objHX, objHY))
    # Locate the matched objects
    indMatchH = np.where(minDistH < tol)
    # Delete the matched objects from the Hot run list
    objHnew = copy.deepcopy(objH)
    objHnew = np.delete(objHnew, indMatchH)

    # Try to account for the difference in size and area of the same
    # object from Cold and Hot run
    objMatchH = objH[(minDistH < tol) & (np.log10(objH['npix']) < 2.6)]
    # Do something similar to the Cold run list
    minDistC = np.asarray(map(lambda x, y: np.min(np.sqrt((x-objHX)**2.0
                                                         + (y-objHY)**2.0)),
                             objCX, objCY))
    indMatchC = np.where(minDistC < tol)
    objMatchC = objC[(minDistC < tol) & (np.log10(objC['npix']) < 2.6)]

    # This is designed to be only a rough match
    # Only works for not very big galaxies, will fail for stars
    aRatio = (np.nanmedian(objMatchC['a']) /
              np.nanmedian(objMatchH['a']))
    objHnew['a'] *= aRatio
    objHnew['b'] *= aRatio

    nRatio = (np.nanmedian(objMatchC['npix']) /
              np.nanmedian(objMatchH['npix']))
    objHnew['npix'] *=  nRatio

    objComb = np.concatenate((objC, objHnew))

    return objComb, objHnew

def getConvKernel(kernel):
    """
    Convolution kernel for the SEP detections
    """
    if kernel is 1:
        # Tophat_3.0_3x3
        convKer = np.asarray([[0.560000, 0.980000, 0.560000],
                              [0.980000, 1.000000, 0.980000],
                              [0.560000, 0.980000, 0.560000]])
    elif kernel is 2:
        # Topcat_4.0_5x5
        convKer = np.asarray([[0.000000, 0.220000, 0.480000, 0.220000, 0.000000],
                              [0.220000, 0.990000, 1.000000, 0.990000, 0.220000],
                              [0.480000, 1.000000, 1.000000, 1.000000, 0.480000],
                              [0.220000, 0.990000, 1.000000, 0.990000, 0.220000],
                              [0.000000, 0.220000, 0.480000, 0.220000, 0.000000]])
    elif kernel is 3:
        # Topcat_5.0_5x5
        convKer = np.asarray([[0.150000, 0.770000, 1.000000, 0.770000, 0.150000],
                              [0.770000, 1.000000, 1.000000, 1.000000, 0.770000],
                              [1.000000, 1.000000, 1.000000, 1.000000, 1.000000],
                              [0.770000, 1.000000, 1.000000, 1.000000, 0.770000],
                              [0.150000, 0.770000, 1.000000, 0.770000, 0.150000]])
    elif kernel is 4:
        # Gaussian_3.0_5x5
        convKer = np.asarray([[0.092163, 0.221178, 0.296069, 0.221178, 0.092163],
                              [0.221178, 0.530797, 0.710525, 0.530797, 0.221178],
                              [0.296069, 0.710525, 0.951108, 0.710525, 0.296069],
                              [0.221178, 0.530797, 0.710525, 0.530797, 0.221178],
                              [0.092163, 0.221178, 0.296069, 0.221178, 0.092163]])
    elif kernel is 5:
        # Gaussian_4.0_7x7
        convKer = np.asarray([[0.047454, 0.109799, 0.181612, 0.214776, 0.181612, 0.109799, 0.047454],
                              [0.109799, 0.254053, 0.420215, 0.496950, 0.420215, 0.254053, 0.109799],
                              [0.181612, 0.420215, 0.695055, 0.821978, 0.695055, 0.420215, 0.181612],
                              [0.214776, 0.496950, 0.821978, 0.972079, 0.821978, 0.496950, 0.214776],
                              [0.181612, 0.420215, 0.695055, 0.821978, 0.695055, 0.420215, 0.181612],
                              [0.109799, 0.254053, 0.420215, 0.496950, 0.420215, 0.254053, 0.109799],
                              [0.047454, 0.109799, 0.181612, 0.214776, 0.181612, 0.109799, 0.047454]])
    elif kernel is 6:
        # Gaussian_5.0_9x9
        convKer = np.asarray([[0.030531, 0.065238, 0.112208, 0.155356, 0.173152, 0.155356, 0.112208, 0.065238, 0.030531],
                              [0.065238, 0.139399, 0.239763, 0.331961, 0.369987, 0.331961, 0.239763, 0.139399, 0.065238],
                              [0.112208, 0.239763, 0.412386, 0.570963, 0.636368, 0.570963, 0.412386, 0.239763, 0.112208],
                              [0.155356, 0.331961, 0.570963, 0.790520, 0.881075, 0.790520, 0.570963, 0.331961, 0.155356],
                              [0.173152, 0.369987, 0.636368, 0.881075, 0.982004, 0.881075, 0.636368, 0.369987, 0.173152],
                              [0.155356, 0.331961, 0.570963, 0.790520, 0.881075, 0.790520, 0.570963, 0.331961, 0.155356],
                              [0.112208, 0.239763, 0.412386, 0.570963, 0.636368, 0.570963, 0.412386, 0.239763, 0.112208],
                              [0.065238, 0.139399, 0.239763, 0.331961, 0.369987, 0.331961, 0.239763, 0.139399, 0.065238],
                              [0.030531, 0.065238, 0.112208, 0.155356, 0.173152, 0.155356, 0.112208, 0.065238, 0.030531]])
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
                       growC=7.0, growW=3.5, growH=4.0, kernel=4,
                       galX=None, galY=None, galR1=None, galR2=None, galR3=None,
                       galQ=None, galPA=None, visual=True):
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
    # Whether the center of the galaxy is provided; If not, assume that
    # galaxy center is located at the image center
    if galX is None:
        galX = (dimX / 2.0)
    if galY is None:
        galY = (dimY / 2.0)

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
    convKerC = getConvKernel(kernel)
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
    convKerH = getConvKernel(kernel)
    # Hot Run
    objH = sep.extract(imgSubH, thrH, minarea=minDetH, conv=convKerH,
                       deblend_nthresh=debThrH, deblend_cont=debConH)
    if verbose:
        print "### %d objects have been detected in the Hot Run" % objH['x'].shape[0]

    #
    # Merge the object lists from Cold and


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="Prefix of the cutout image files")
    parser.add_argument('-r', '--root', dest='root', help='Path to the image files',
                        default=None)
    args = parser.parse_args()

    coaddCutoutPrepare(args.prefix, root=args.root)
