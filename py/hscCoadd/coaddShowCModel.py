#!/usr/bin/env python
# encoding: utf-8

from __future__ import division

import os
import copy
import argparse

import numpy as np

# Matplotlib related
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
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
plt.ioff()

# Astropy
from astropy.io import fits
from astropy import units as u
from astropy.stats import sigma_clip
from astropy.wcs import WCS

# Colormap
from palettable.colorbrewer.sequential import Greys_3 as pcmap1
cmap5 = pcmap1.mpl_colormap
from palettable.colorbrewer.diverging  import RdYlGn_11 as pcmap2
cmap6 = pcmap2.mpl_colormap

# Personal
import hscUtils as hUtil


def showCmodel(imgData, xUse, yUse, ellipse, colorCmod, figSize=14, fontSize=14,
              filter='HSC-I', ellipName='Exponential', showSource=True,
              mag0=24.5, mag1=18.0, figName='showCmodel.png'):

    rEllip, eEllip, paEllip = srcMoments2Ellip(ellipse)
    ellipPlot = getEll2Plot(xUse, yUse, rEllip, eEllip, paEllip)

    fig = plt.figure(figsize=(figSize, figSize))
    fig.subplots_adjust(hspace=0.0, wspace=0.0, left=0.03, bottom=0.03,
                        top=0.95, right=0.995)
    ax = fig.add_subplot(1,1,1)
    fontsize = fontSize
    ax.minorticks_on()

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    ax.set_title('%s-band Image - %s' % (filter, ellipName),
                fontsize=(fontSize+13), fontweight='bold')
    ax.title.set_position((0.5, 1.01))

    imin, imax = hUtil.zscale(imgData, contrast=0.10, samples=500)
    ax.imshow(np.arcsinh(imgData), interpolation="none",
               vmin=imin, vmax=imax, cmap=cmap5)

    if showSource:
        ax.scatter(xUse, yUse, marker='+', s=25, c='r')

    for (e, c) in zip(ellipPlot, colorCmod):
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
        e.set_alpha(0.8)
        e.set_edgecolor(cmap6(int(c)))
        e.set_facecolor('none')
        e.set_linewidth(1.5)

    cax = fig.add_axes([0.14, 0.18, 0.21, 0.02])
    norm = mpl.colors.Normalize(vmin=mag1, vmax=mag0)
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap6, norm=norm, orientation='horizontal')
    cbar.set_label('cModel Magnitude (mag)', fontsize=(fontSize+3))

    ax.set_xlim(0, imgData.shape[1]-1)
    ax.set_ylim(0, imgData.shape[0]-1)

    fig.savefig(figName)
    plt.close(fig)


def srcToRaDec(src):
    """
    Return a list of (RA, DEC)
    """
    raArr  = np.asarray([coord[0] * 180.0 / np.pi for coord in src['coord']])
    decArr = np.asarray([coord[1] * 180.0 / np.pi for coord in src['coord']])

    return raArr, decArr


def toColorArr(data, bottom=None, top=None):
    """
    Convert a data array to "color array" (between 0 and 1)
    """
    if top is not None:
        data[data >= top]    = top
    if bottom is not None:
        data[data <= bottom] = bottom

    return ((data - np.nanmin(data)) / (np.nanmax(data) - np.nanmin(data))) * 255.0


def getEll2Plot(x, y, re, ell, theta):

    a = (re * 2.0)
    b = (re * (1.0 - ell) * 2.0)
    pa = (theta + 90.0)

    ells = [Ellipse(xy=np.array([x[i], y[i]]),
                    width=np.array(b[i]),
                    height=np.array(a[i]),
                    angle=np.array(pa[i]))
            for i in range(x.shape[0])]

    return ells


def srcMoments2Ellip(ellip):
    """
    Translate The 2nd Moments into Elliptical Shape
    """

    Ixx, Iyy, Ixy = ellip[:,0], ellip[:,1], ellip[:,2]

    e1 = (Ixx - Iyy) / (Ixx + Iyy)
    e2 = (2.0 * Ixy / (Ixx + Iyy))
    ell = np.sqrt(e1 ** 2.0 + e2 ** 2.0)

    theta = (0.5 * np.arctan2(e2, e1)) * 180.0 / np.pi
    r2 = np.sqrt(Ixx + Iyy)

    return r2, ell, theta


def getCoaddData(prefix, galId, root=None, filter='HSC-I', ref=False,
                 noParent=False):
    """
    Return the coadd image and source catalog
    """

    if root is None:
        root = ''
    loc = os.path.join(root, galId, filter)
    imgFile = os.path.join(loc, prefix + '_' + galId + '_' + filter + '_full_img.fits')
    if not ref:
        catFile = os.path.join(loc, prefix + '_' + galId + '_' + filter + '_full_meas.fits')
    else:
        catFile = os.path.join(loc, prefix + '_' + galId + '_' + filter + '_full_ref.fits')
    if not os.path.isfile(imgFile) or not os.path.isfile(catFile):
        print "## Image   : %s" % imgFile
        print "## Catalog : %s" % catFile
        raise Exception("### Can not find image or source catalog !")
        return
    else:
        imgData = fits.open(imgFile)[0].data
        imgHead = fits.open(imgFile)[0].header
        catData = fits.open(catFile)[1].data

        if noParent:
            print "Number of detections (with parents) : %6d" % len(catData)
            catUse = catData[(catData['deblend_nchild'] == 0)]
            print "Number of detections : %6d" % len(catUse)
            return imgData, imgHead, catUse
        else:
            return imgData, imgHead, catData


def catFlux2Mag(cat, flux, zeropoint=27.0):
    """
    Conver flux into magnitude
    """
    return (-2.5 * np.log10(cat[flux]) + zeropoint)


def coaddShowCModel(prefix, galId, filter='HSC-I', ref=False, root=None, noParent=False,
                    zeropoint=27.0, mag0=24.5, mag1=18.0, fontSize=14, figSize=14,
                    showSource=True, verbose=True):
    """
    Visualize the cModel fitting results for given coadd image
    """
    galId = str(galId).strip()
    # Get the coadd data, both image and catalog
    imgData, imgHead, catData = getCoaddData(prefix, galId, filter=filter, root=root,
                                             ref=ref)
    # Objects with useful cModel information
    catCModel = catData[(np.isfinite(catData['cmodel_flux'])) &
                    (np.isfinite(catData['cmodel_exp_flux'])) &
                    (np.isfinite(catData['cmodel_dev_flux'])) &
                    (catData['cmodel_flux'] > 0.0) &
                    (catData['deblend_nchild'] == 0) &
                    (catData['classification_extendedness'] >= 0.5)]
    if verbose:
        print "### %d objects with useful cModel information" % len(catCModel)
    # Convert (RA, DEC) into (X, Y)
    wcs = WCS(imgHead)
    raCmod, decCmod = srcToRaDec(catCModel)
    xyCmod = wcs.wcs_world2pix((catCModel['coord'] * 180.0 / np.pi), 1)
    xCmod, yCmod = xyCmod[:,0], xyCmod[:,1]
    # Get cModelMagnitude
    magCmod = catFlux2Mag(catCModel, 'cmodel_flux', zeropoint=zeropoint)
    # Convert into color array
    colorCmod = toColorArr(magCmod, top=mag0, bottom=mag1)
    #
    ellipList = ['cmodel_exp_ellipse', 'cmodel_dev_ellipse', 'shape_sdss']
    ellipType = ['Exponential', 'de Vacouleur', 'SDSS Shape']
    for (ii, ellipName) in enumerate(ellipList):
        if verbose:
            print "### Start to work on %s model !" % ellipName
        ellipse = catCModel[ellipName]
        loc = os.path.join(root, galId, filter)
        pngFile = os.path.join(loc, prefix + '_' + galId + '_' + filter + \
                              '_' + ellipName + '.png')
        showCmodel(imgData, xCmod, yCmod, ellipse, colorCmod,
                  figSize=figSize, fontSize=fontSize, filter=filter,
                  ellipName=ellipType[ii], showSource=showSource,
                  mag0=mag0, mag1=mag1, figName=pngFile)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="Prefix of the image files")
    parser.add_argument("id", help="Galaxy Id")
    parser.add_argument("-f", "--filter", dest="filter",
            default='HSC-I', help="HSC Filter")
    parser.add_argument("-r", "--root", dest="root",
            default='', help="Root to the cutout images")
    parser.add_argument("-z", "--zeropoint", type=float, default=27.0,
            dest='zeropoint', help="Photometric zeropoint")
    parser.add_argument("--mag0", type=float, default=24.5,
            dest='mag0', help="Faint limit of cModel magnitude")
    parser.add_argument("--mag1", type=float, default=18.0,
            dest='mag1', help="Bright limit of cModel magnitude")
    parser.add_argument("--fontsize", type=float, default=18.0,
            dest='fontsize', help="Font size")
    parser.add_argument("--figsize", type=float, default=14.0,
            dest='figsize', help="Figure size")
    parser.add_argument("--ref", help="Whether to use the deepCoadd_ref catalog",
            dest='ref', action='store_true', default=False)
    parser.add_argument("--nosource", help="Don't show the locations of detections",
            dest='nosource', action='store_false', default=True)
    args = parser.parse_args()

    coaddShowCModel(args.prefix, args.id, filter=args.filter, ref=args.ref,
                    root=args.root, zeropoint=args.zeropoint,
                    mag0=args.mag0, mag1=args.mag1, fontSize=args.fontsize,
                    figSize=args.figsize, showSource=args.nosource)
