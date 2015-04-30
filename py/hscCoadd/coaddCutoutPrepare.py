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


def showSEPImage(image, contrast=0.2, size=10, cmap=cmap1,
                 title='Image', pngName='sep.png', titleInside=True,
                 ellList1=None, ellList2=None, ellList3=None,
                 ellColor1='b', ellColor2='r', ellColor3='g',
                 ax=None, mask=None):
    """
    Visualization of the results

    """
    fig = plt.figure(figsize=(size, size))
    fig.subplots_adjust(hspace=0.0, wspace=0.0,
                        bottom=0.08, left=0.08,
                        top=0.92, right=0.98)
    ax = fig.add_axes([0.000, 0.002, 0.996, 0.996])
    fontsize = 16
    ax.minorticks_on()

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    ax.set_title(title, fontsize=28, fontweight='bold', color='r')
    if not titleInside:
        ax.title.set_position((0.5, 1.01))
    else:
        ax.title.set_position((0.5, 0.90))

    imcopy = copy.deepcopy(image)
    imin, imax = hUtil.zscale(imcopy, contrast=contrast, samples=500)
    if mask is not None:
        imcopy[mask > 0] = np.nan

    ax.imshow(np.arcsinh(imcopy), interpolation="none",
               vmin=imin, vmax=imax, cmap=cmap)

    if ellList1 is not None:
        for e in ellList1:
            ax.add_artist(e)
            e.set_clip_box(ax.bbox)
            e.set_alpha(0.8)
            e.set_edgecolor(ellColor1)
            e.set_facecolor('none')
            e.set_linewidth(1.5)

    if ellList2 is not None:
        for e in ellList2:
            ax.add_artist(e)
            e.set_clip_box(ax.bbox)
            e.set_alpha(0.8)
            e.set_edgecolor(ellColor2)
            e.set_facecolor('none')
            e.set_linewidth(1.5)

    if ellList3 is not None:
        for e in ellList3:
            ax.add_artist(e)
            e.set_clip_box(ax.bbox)
            e.set_alpha(0.8)
            e.set_edgecolor(ellColor3)
            e.set_facecolor('none')
            e.set_linewidth(1.5)

    fig.savefig(pngName)
    plt.close(fig)


def addFlag(dictName, flagName, flagValue):
    """
    Add a flag to a dictionary

    """
    # If the dictionary for flags has not been defined, make a new one
    try:
        dictName
    except NameError:
        dictName = np.array([], dtype=[('name', 'a20'), ('value', 'i1')])
    # Assign a new flag
    newFlag = (str(flagName), flagValue)
    newDict = np.insert(dictName, 0, newFlag)

    return newDict


def objList2Reg(objs, regName='ds9.reg', color='Blue'):
    """
    Save the Object List to DS9 Region File
    """
    # DS9 region file header
    head1 = '# Region file format: DS9 version 4.1\n'
    head2 = 'global color=%s width=2\n' % color
    head3 = 'image\n'
    # Open the output file and write the header
    regFile = open(regName, 'w')
    regFile.write(head1)
    regFile.write(head2)
    regFile.write(head3)
    # Save the data
    # The format for ellipse is: "ellipse x y radius radius angle"
    for obj in objs:
        ellLine = 'ellipse ' + str(obj['x']) + ' ' + str(obj['y']) + ' ' + \
                  str(obj['a'] / 2.0) + ' ' + str(obj['b'] / 2.0) + ' ' + \
                  str(obj['theta'] * 180.0 / np.pi) + '\n'
        regFile.write(ellLine)
    # Close the file
    regFile.close()


def saveFits(img, fitsName, head=None, clobber=True):
    """
    Save an image to FITS file

    """
    imgHdu = fits.PrimaryHDU(img)
    if head is not None:
        imgHdu.header = head
    imgHdu.writeto(fitsName, clobber=clobber)


def saveSEPObjects(objs, prefix='sep_objects', csv=True,
                   reg=True, verbose=True, color='red'):
    """
    Save the properties of objects that are extracted by SEP into
    a cPickle file or .csv file, .deg file

    """
    # 1. Save a .pkl file
    pklFile = prefix + '.pkl'
    hUtil.saveToPickle(objs, pklFile)
    if os.path.isfile(pklFile):
        if verbose:
            print "### Save object list to .pkl file: %s" % pklFile
    else:
        raise Exception("### Something is wrong with the .pkl file")
    # 2. Save a .csv file
    if csv:
        csvFile = prefix + '.csv'
        hUtil.saveToCSV(objs, csvFile)
        if os.path.isfile(csvFile):
            if verbose:
                print "### Save object list to .csv file: %s" % csvFile
        else:
            raise Exception("### Something is wrong with the .csv file")
    # 3. Save a .deg file
    if reg:
        regFile = prefix + '.reg'
        objList2Reg(objs, regName=regFile, color=color)
        if os.path.isfile(regFile):
            if verbose:
                print "### Save object list to .reg file: %s" % regFile
        else:
            raise Exception("### Something is wrong with the .reg file")


def adaptiveMask(objC, a=2.0, b=1.5, c=4.0, seeing=1.0,
                 pix=0.168, verbose=False):
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


def combObjCat(objCold, objHot, tol=1.0):
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


def getFluxRadius(img, objs, maxSize=16.0, subpix=5, byteswap=True):
    """
    Given the original image, the detected objects, using SEP
    to measure different flux radius: R20, R50, R90

    """
    # Make a copy of the image, and byteswap it if necessary
    imgOri = copy.deepcopy(img)
    if byteswap:
        imgOri = imgOri.byteswap(True).newbyteorder()
    # Get the flux radius
    rflux, flag = sep.flux_radius(imgOri, objs['x'], objs['y'],
                                  maxSize*objs['a'], [0.2, 0.5, 0.9],
                                  normflux=objs['cflux'], subpix=subpix)
    r20 = np.array([rr[0] for rr in rflux])
    r50 = np.array([rr[1] for rr in rflux])
    r90 = np.array([rr[2] for rr in rflux])

    return r20, r50, r90


def objDistTo(objs, cenX, cenY, usePeak=False, convol=False):
    """
    Get the distance of objects from SEP to a reference point on the image

    """
    if usePeak:
        if convol:
            xc, yc = objs['xcpeak'], objs['ycpeak']
        else:
            xc, yc = objs['xpeak'], objs['ypeak']
    else:
        xc, yc = objs['x'], objs['y']

    return np.sqrt((xc - cenX)**2 + (yc - cenY)**2)


def readCutoutHeader(imgHead, pixDefault=0.168,
                     zpDefault=27.0):
    """
    Read the pixel scale, image size, and photometric zeropoint form
    the image header

    TODO: Make it more generic, right now it is only for HSC
     * pixel scale can be read from the WCS information
    XXX: Right now, the TOTEXPT is not working for HSC
    """
    # Get the pixel scale of the image
    try:
        pixScaleX = pixScaleY = imgHead['PIXEL']
    except:
        print "### Pixel scale keyword is not available in the header"
        print "### Default %6.3f arcsec/pixel value is adopted" % pixDefault
        pixScaleX = pixScaleY = pixDefault
    # Get the image size
    imgSizeX = imgHead['NAXIS1']
    imgSizeY = imgHead['NAXIS2']
    # Get the photometric zeropoint
    try:
        photZP = imgHead['PHOTZP']
    except:
        print "### PHOZP keyword is not available in the header"
        print "### Default value of %5.2f is adopted!" % zpDefault
        photZP = zpDefault
    # Total exptime
    try:
        expTot = imgHead['TOTEXPT']
    except:
        print "### TOTEXPT keyword is not available in the header"
        print "### Use 1.0 sec instead"
        expTot = 1.0

    return pixScaleX, pixScaleY, imgSizeX, imgSizeY, photZP, expTot


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
                       bSizeH=8, bSizeC=80, thrH=3.5, thrC=1.5, mask=1,
                       growC=4.0, growW=3.0, growH=1.5, kernel=4, central=1,
                       galX=None, galY=None, galR1=None, galR2=None, galR3=None,
                       galQ=None, galPA=None, visual=True, suffix='',
                       badMsk=None, skyClip=3.0, rebin=4):
    """
    The structure of the cutout has been changed.  Now the cutout procedure
    will generate separated files for Image, Bad Mask, Detection Plane, and
    Variance (also Sigma) images.  Souce catalogs can also be made available.

    Right now, this new format is only available for the coaddImageCutFull()
    function; coaddImageCutout() will be modified later to also adopt this
    format
    """

    # 0. Get necessary information
    # Read the input cutout image
    imgArr, imgHead, mskArr, detArr = readCutoutImage(prefix, root=root)
    if verbose:
        print "#############################################################"
        print "### DEAL WITH IMAGE : %s" % (prefix + '_img.fits')
    if detArr is None:
        detFound = False
    else:
        detFound = True
    pixX, pixY, dimX, dimY, photZP, expTot = readCutoutHeader(imgHead)
    if verbose:
        print "###    The pixel scale in X/Y directions " + \
                "are %7.4f / %7.4f arcsecs" % (pixX, pixY)
        print "###    The image size in X/Y directions " + \
                "are %d / %d pixels" % (dimX, dimY)
        print "###                      %10.2f / %10.2f arcsecs" % (dimX * pixX, dimY * pixY)
        print "###    The photometric zeropoint is %6.2f " % photZP
        #print "### The total exposure time for the center is %6.1f secs" % expTot
    # Whether the center of the galaxy is provided; If not, assume that
    # galaxy center is located at the image center
    if galX is None:
        galX = (dimX / 2.0)
    if galY is None:
        galY = (dimY / 2.0)
    if verbose:
        print "###    The Galaxy Center is assumed at %6.1f, %6.1f" % (galX, galY)
    # Suffix
    if (suffix is not '') and (suffix[-1] is not '_'):
        suffix = suffix + '_'

    # 1. Get the backgrounds
    """
    Construct "background" images with different size and filters using SEP,
    and subtract these background off before extract objects

    The SEP detections will be run in two-modes:
        Cold: relative global background; low-detection threshold
        Hot:  very local background; median-detection threshold
    """
    if verbose:
        print "### 1. BACKGROUND SUBTRACTION USING SEP"
    # Cold Background Run
    fSizeC = int(bSizeC / 4)
    imgC = imgByteSwap(imgArr)
    bkgC, imgSubC = sepGetBkg(imgC, bkgSize=bSizeC, bkgFilter=fSizeC)
    if verbose:
        print "###    Cold Background -- Avg: %9.5f " % bkgC.globalback + \
              "       RMS: %9.5f " % bkgC.globalrms
    # Hot Background Run
    fSizeH = int(bSizeH / 2)
    imgH = imgByteSwap(imgArr)
    bkgH, imgSubH = sepGetBkg(imgH, bkgSize=bSizeH, bkgFilter=fSizeH)
    if verbose:
        print "###    Hot Background  -- Avg: %9.5f " % bkgH.globalback + \
              "       RMS: %9.5f " % bkgH.globalrms
    if visual:
        # Fig.a
        bkgPNG1 = prefix + '_' + suffix + 'bkgC.png'
        showSEPImage(bkgC.back(), contrast=0.3, title='Background - Cold Run',
                     pngName=bkgPNG1)
        # Fig.b
        bkgPNG2 = prefix + '_' + suffix + 'bkgH.png'
        showSEPImage(bkgH.back(), contrast=0.3, title='Background - Hot Run',
                     pngName=bkgPNG2)

    # 2. Object detections
    """
    Use SEP to extract information of detected objects
    """
    """ Cold Run """
    if verbose:
        print "### 2. DETECT OBJECTS USING SEP"
    # Parameters for cold run
    if verbose:
        print "###     Cold run detection threshold: %4.1f" % thrC
    thrC *= bkgC.globalrms  # Cold run detection threshold
    minDetC = 5             # Minimun adjacent pixels for detection for cold run
    debThrC = 20            # Deblending threshold for cold run (32)
    debConC = 0.002         # Minimum deblending contrast ratio for cold run (0.005)
    convKerC = getConvKernel(kernel)
    # Cold Run
    objC = sep.extract(imgSubC, thrC, minarea=minDetC, conv=convKerC,
                       deblend_nthresh=debThrC, deblend_cont=debConC)
    if verbose:
        print "###     %d objects have been detected in the Cold Run" % objC['x'].shape[0]
    # Save objects list to different format of files
    prefixC = prefix + '_' + suffix + 'objC'
    saveSEPObjects(objC, prefix=prefixC, color='Blue')
    # Calculate the object-galaxy center distance
    cenDistC = objDistTo(objC, galX, galY)

    """ Hot Run """
    # Parameters for hot run
    if verbose:
        print "###     Hot run detection threshold: %4.1f" % thrH
    thrH *= bkgH.globalrms  # Hot run detection threshold
    minDetH = 5             # Minimun adjacent pixels for detection for hot run
    debThrH = 20            # Deblending threshold for hot run (32)
    debConH = 0.002         # Minimum deblending contrast ratio for hot run (0.005)
    convKerH = getConvKernel(kernel)
    # Hot Run
    objH = sep.extract(imgSubH, thrH, minarea=minDetH, conv=convKerH,
                       deblend_nthresh=debThrH, deblend_cont=debConH)
    if verbose:
        print "###     %d objects have been detected in the Hot Run" % objH['x'].shape[0]
    # Save objects list to different format of files
    prefixH = prefix + '_' + suffix + 'objH'
    saveSEPObjects(objH, prefix=prefixH, color='Red')
    # Calculate the object-galaxy center distance
    cenDistH = objDistTo(objH, galX, galY)
    if visual:
        # Fig.c
        objPNG1 = prefix + '_' + suffix + 'objC.png'
        objEllC = getEll2Plot(objC, radius=(objC['a']*growC))
        showSEPImage(imgSubC, contrast=0.06, title='Detections - Cold Run',
                     pngName=objPNG1, ellList1=objEllC, ellColor1='b')
        # Fig.d
        objPNG2 = prefix + '_' + suffix + 'objH.png'
        objEllH = getEll2Plot(objH, radius=(objH['a']*growW))
        showSEPImage(imgSubH, contrast=0.10, title='Detections - Hot Run',
                     pngName=objPNG2, ellList1=objEllH, ellColor1='r')

    # 3. Merge the objects from Cold and Hot runs together
    if verbose:
        print "### 3. COMBINE OBJECTS FROM COLD AND HOT RUN "
    # Merge the object lists from Cold and
    objComb, objHnew = combObjCat(objC, objH, tol=1.0)
    # Also save the combined object lists
    prefixComb = prefix + '_' + suffix + 'objComb'
    saveSEPObjects(objComb, prefix=prefixComb, color='Green')
    # Calculate the object-galaxy center distance
    cenDistComb = objDistTo(objComb, galX, galY)
    if verbose:
        print "###    %d objects are left in the combined list" % len(objComb)
    if visual:
        # Fig.e
        objPNG3 = prefix + '_' + suffix + 'objComb.png'
        objEllComb = getEll2Plot(objComb, radius=(objComb['a']*growC))
        showSEPImage(imgSubC, contrast=0.06, title='Detections - Combined',
                     pngName=objPNG3, ellList1=objEllComb, ellColor1='orange')

    # 4. Extract Different Flux Radius: R20, R50, R90 for every objects
    if verbose:
        print "### 4. EXTRACTING R20, R50, R90 OF EACH OBJECTS "
    r20, r50, r90 = getFluxRadius(imgArr, objComb, maxSize=16.0, subpix=5)
    if visual:
        # Fig.f
        objPNG4 = prefix + '_' + suffix + 'objRad.png'
        objEllR20 = getEll2Plot(objComb, radius=r20)
        objEllR50 = getEll2Plot(objComb, radius=r50)
        objEllR90 = getEll2Plot(objComb, radius=r90)
        showSEPImage(imgSubC, contrast=0.20, title='Flux Radius: R20/R50/R90',
                     pngName=objPNG4, ellList1=objEllR20, ellColor1='r',
                     ellList2=objEllR50, ellColor2='orange',
                     ellList3=objEllR90, ellColor3='b')

    # 5. Mask all objects on the image
    if verbose:
        print "### 5. MASKING OUT ALL OBJECTS ON THE IMAGE "
    mskAll = np.zeros(imgSubC.shape, dtype='uint8')
    if mask == 1:
        # TODO: This is still not idea, even using flux radius, should take
        #       the compactness (R90/R50) and (R50/R20) into account
        # Make a ALL_OBJECT mask using the r90
        sep.mask_ellipse(mskAll, objComb['x'], objComb['y'],
                         objComb['a'], objComb['b'],
                         objComb['theta'], r=growC)
        # Make a new objList, and change the size
        objMskAll = copy.deepcopy(objComb)
        objMskAll['a'] *= growC
        objMskAll['b'] *= growC
    elif mask == 2:
        # Grow the cold run detections using the adaptive method
        # Use the new size to mask out all objects
        adGrowC = adaptiveMask(objComb, a=2.2)
        # Build a ALL_OBJECT Mask
        sep.mask_ellipse(mskAll, objComb['x'], objComb['y'],
                objComb['a'], objComb['b'],
                objComb['theta'], r=adGrowC)
        # Make a new objList, and change the size
        objMskAll = copy.deepcopy(objComb)
        objMskAll['a'] *= adGrowC
        objMskAll['b'] *= adGrowC
    else:
        raise Exception("mask == 1 or mask == 2")
    # Save the mask to FITS
    mskAllFile = prefix + '_' + suffix + 'mskall.fits'
    saveFits(mskAll, mskAllFile, head=imgHead)
    # Save the Objlist using the growed size
    prefixM = prefix + '_' + suffix + 'mskall'
    saveSEPObjects(objMskAll, prefix=prefixM, color='Blue')
    if visual:
        # Fig.f
        mskPNG1 = prefix + '_' + suffix + 'mskall.png'
        showSEPImage(imgSubC, contrast=0.75, title='Mask - All Objects',
                     pngName=mskPNG1, mask=mskAll)
    # TODO: Use this mask to estimate the average background level

    # 6. Find the central galaxy, get its shape and radius, remove it
    if verbose:
        print "### 6. CLEAR THE CENTRAL REGION AROUND THE GALAXY"
        print "###  6.1. LOCATE THE CENTRAL GALAXY IN THE OBJECT LIST"
    # Find the main galaxy from the cold run
    # The index of the "central" galaxy
    # TODO: Check whether or not it is truely the central galaxy
    cenObjIndex = np.argmin(cenDistComb)
    # Get its shape and size
    if verbose:
        print "###  6.2. ESTIMATE THE B/A AND PA OF THE GALAXY"
    if galQ is None:
        galQ  = (objComb[cenObjIndex]['b'] / objComb[cenObjIndex]['a'])
    if galPA is None:
        galPA = (objComb[cenObjIndex]['theta'] * 180.0 / np.pi)
    if verbose:
        print "###    (b/a) of the galaxy: %6.2f" % galQ
        print "###      PA  of the galaxy: %6.1f" % galPA
    # TODO: Make it possible to scale r20 to galR1, r50 to galR2,
    #       r90 to galR3
    if verbose:
        print "###  6.3. ESTIMATING THE GAL_R1/R2/R3"
    if galR1 is None:
        galR1 = (r50[cenObjIndex] * 1.5)
    if galR2 is None:
        galR2 = (r90[cenObjIndex] * 1.0)
    if galR3 is None:
        galR3 = (r90[cenObjIndex] * 4.0)
    if verbose:
        print "###    galR1: %7.2f" % galR1
        print "###    galR2: %7.2f" % galR2
        print "###    galR3: %7.2f" % galR3
    # TODO: Check if these R1, R2, R3 are reasonbale, compared to image size

    # 7. Remove the central object (or clear the central region)
    #    Separate the objects into different group and mask them out using
    #    different growth ratio
    if verbose:
        print "###  6.4. CLEAR A REGION AROUND CENTRAL GALAXY"
    objNoCen  = copy.deepcopy(objComb)
    r90NoCen  = copy.deepcopy(r90)
    distNoCen = copy.deepcopy(cenDistComb)
    if central == 1:
        # Remove the central objects from the list and r90 array
        objNoCen  = np.delete(objNoCen,  cenObjIndex)
        r90NoCen  = np.delete(r90NoCen,  cenObjIndex)
        distNoCen = np.delete(distNoCen, cenObjIndex)
    elif central == 2:
        # Remove all objects within certain radii to the center of galaxy
        # TODO: Now, circular distance is used; Should take the
        #       geometry of the galaxy into account
        indCen = np.where(cenDistComb < (growW * r50))
        if verbose:
            print "###    %d objects are found in the central region" % len(indCen)
        objNoCen  = np.delete(objNoCen,  indCen)
        r90NoCen  = np.delete(r90NoCen,  indCen)
        distNoCen = np.delete(distNoCen, indCen)

    # 8. Separate the rest objects into different groups according to
    #    their distance to the central galaxy
    if verbose:
        print "### 7. GENERATING THE FINAL MASK"
    # Index of objects in different groups
    indG1 = (distNoCen <= galR2)
    indG2 = (distNoCen > galR2) & (distNoCen < galR3)
    indG3 = (distNoCen > galR3)
    # Isolate them into different group
    objG1 = objNoCen[indG1]
    objG2 = objNoCen[indG2]
    objG3 = objNoCen[indG3]
    # Generating final mask
    mskG1 = mskG2 = mskG3 = np.zeros(imgArr.shape, dtype='uint8')
    sep.mask_ellipse(mskG1, objG1['x'], objG1['y'], r90NoCen[indG1],
                    (r90NoCen[indG1] * objG1['b'] / objG1['a']),
                     objG1['theta'], r=growH)
    sep.mask_ellipse(mskG2, objG2['x'], objG2['y'], r90NoCen[indG2],
                    (r90NoCen[indG2] * objG2['b'] / objG2['a']),
                     objG2['theta'], r=growW)
    sep.mask_ellipse(mskG3, objG3['x'], objG3['y'], r90NoCen[indG3],
                    (r90NoCen[indG3] * objG3['b'] / objG3['a']),
                     objG3['theta'], r=growC)
    mskFinal = (mskG1 | mskG2 | mskG3)
    # Save the mask to FITS file
    # TODO: Add an option to combine with HSC BAD MASK
    mskFinFile = prefix + '_' + suffix + 'mskfin.fits'
    saveFits(mskFinal, mskFinFile, head=imgHead)
    # Save the Objlist using the growed size
    prefixF = prefix + '_' + suffix + 'mskfin'
    objFin = copy.deepcopy(objNoCen)
    objFin[indG1]['a'] *= growH
    objFin[indG1]['b'] *= growH
    objFin[indG2]['a'] *= growW
    objFin[indG2]['b'] *= growW
    objFin[indG3]['a'] *= growC
    objFin[indG3]['b'] *= growC
    saveSEPObjects(objFin, prefix=prefixM, color='Blue')
    if visual:
        # Fig.g
        mskPNG2 = prefix + '_' + suffix + 'mskfin.png'
        showSEPImage(imgArr, contrast=0.75, title='Mask - Final',
                     pngName=mskPNG2, mask=mskFinal)

    # 9. Estimate the global background level
    if verbose:
        print "### 8: ESTIMATING BACKGROUND AND SURFACE BRIGHTNESS LIMIT"
    # Pixel values of all pixels that are not masked out (before rebinned)
    pixels = imgArr[mskAll == 0].flatten()
    pixNoMsk = sigma_clip(pixels, skyClip, 3)
    # Get the basic statistics of the global sky
    meanSky1, stdSky1 = np.nanmean(pixNoMsk), np.nanstd(pixNoMsk)
    medSky1 = np.nanmedian(pixNoMsk)
    sbExp1 = getSbpValue(stdSky1, pixX, pixY, zp=photZP)
    if verbose:
        print "###  8.1: Before Rebin the Image "
        print "###    Median Sky: %8.5f" % medSky1
        print "###      Mean Sky: %8.5f" % meanSky1
        print "###    StdDev Sky: %8.5f" % stdSky1
        print "###    SB. Expect: %8.2f" % sbExp1
    # Rebin image
    # XXX: This is still very specific for HSC cutout
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
    sbExp2 = getSbpValue(stdSky2, pixX*rebin, pixY*rebin, zp=photZP)
    if verbose:
        print "###  8.2: After Rebin the Image "
        print "###    Median Sky: %8.5f" % medSky2
        print "###      Mean Sky: %8.5f" % meanSky2
        print "###    StdDev Sky: %8.5f" % stdSky2
        print "###    SB. Expect: %8.2f" % sbExp2
    if visual:
        skyPNG = prefix + '_' + suffix + 'skyhist.png'
        showSkyHist(pixNoMskBin, skypix2=pixNoMsk, pngName=skyPNG)

    # 10. Visualize the detected objects, and find the ones need to be fit
    if verbose:
        print "### 10. SELECTING OBJECTS NEED TO BE FIT"
    # TODO: Plot the basic properties of the detected objects
    # log(flux) v.s. log(R50) v.s. log(R90/R50) v.s. log(cenDist)
    if visual:
        a = 1


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="Prefix of the cutout image files")
    parser.add_argument('-r', '--root', dest='root', help='Path to the image files',
                        default=None)
    args = parser.parse_args()

    coaddCutoutPrepare(args.prefix, root=args.root)
