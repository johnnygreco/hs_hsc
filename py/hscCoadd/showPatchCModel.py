from __future__ import division

import re
import copy
import argparse
import collections
import numpy as np
from distutils.version import StrictVersion

# Matplotlib default settings
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

import astropy.table
from astropy.stats import sigma_clip

from palettable.colorbrewer.sequential import Greys_3 as pcmap1
cmap1 = pcmap1.mpl_colormap
from palettable.colorbrewer.diverging  import RdYlGn_11 as pcmap2
cmap2 = pcmap2.mpl_colormap

import lsst.daf.persistence as dafPersist
import lsst.afw.geom.ellipses
import lsst.afw.geom
import lsst.pex.exceptions
from lsst.afw.table import SourceCatalog, SchemaMapper


def getMag(flux, fluxerr, zeropoint):

    """
    return the magnitude and error
    """
    mag, magerr = -2.5 * np.log10(flux), 2.5/np.log(10.0)*fluxerr/flux
    return (mag.T + zeropoint).T, magerr


def getEllipse(quad):

    """
    returns the semi-major axis, axes ratio and PA for a given quadrupole moment
    """
    e = lsst.afw.geom.ellipses.Axes(quad)
    return e.getA(), e.getB()/e.getA(), e.getTheta() * 180.0/np.pi


def getAstroTable(src, mags=True, zeropoint=27.0):

    """
    returns an astropy table with all the src entries
    if the entries are complex objects, it breaks them down:
      ellipse entries are broken into
           ellipse_a = semi-major axis
           ellipse_q = axis ratio (always < 1)
           ellipse_theta = rotation of semi-major axis from chip x-axis in degrees
    if mags is True, returns the magnitudes for all the flux columns
    """

    tab = astropy.table.Table()
    for name in src.schema.getNames():
        #for reasons I don't understand a lookup by name is much slower than a lookup by key
        nameKey = src.schema.find(name).getKey()
        try:
            tab.add_column(astropy.table.Column(name=name,
                                                data=src.get(nameKey)))
        except lsst.pex.exceptions.LsstException:
            if type(src[0].get(nameKey)) is lsst.afw.geom.ellipses.ellipsesLib.Quadrupole:
                reff, q, theta = zip(*[getEllipse(s.get(nameKey)) for s in src])
                tab.add_column(astropy.table.Column(name=name+'_a', data=reff))
                tab.add_column(astropy.table.Column(name=name+'_q', data=q))
                tab.add_column(astropy.table.Column(name=name+'_theta', data=theta))
            elif type(src[0].get(nameKey)) is lsst.afw.coord.coordLib.IcrsCoord:
                x, y= zip(*[(s.get(nameKey).getRa().asDegrees(),
                             s.get(nameKey).getDec().asDegrees()) for s in src])
                tab.add_column(astropy.table.Column(name=name+'_ra', data=x))
                tab.add_column(astropy.table.Column(name=name+'_dec', data=y))
            else:
                tab.add_column(astropy.table.Column(name=name,
                                                    data=np.array([s.get(nameKey) for s in src])))
            #report angles in degrees
        if isinstance(src[0].get(nameKey), lsst.afw.geom.Angle):
            tab.remove_column(name)
            tab.add_column(astropy.table.Column(data=[s.get(nameKey).asDegrees()
                                                      for s in src],
                                                dtype=float, name=name))

    if mags:
        #this is a horrible hack, but I don't think we can use the slots, since
        #not all the fluxes end up in the slots
        for col in tab.colnames:
            if (re.match('^flux\.[a-z]+$', col) or
                re.match('^flux\.[a-z]+.apcorr$', col) or
                re.match('^cmodel.+flux$', col) or
                re.match('^cmodel.+flux.apcorr$', col)):
                mag, magerr = getMag(tab[col], tab[col+'.err'],
                                     zeropoint if not re.search('apcorr', col) else 0.0)

                tab.add_column(astropy.table.Column(name=re.sub('flux', 'mag', col),
                                                    data=mag))
                tab.add_column(astropy.table.Column(name=re.sub('flux', 'mag', col+'.err'),
                                                    data=magerr))

    return tab


def getCalexpCat(root, tract, patch, filter, calexp=False):

    # make a butler and specify your dataId
    butler = dafPersist.Butler(root)
    dataId = {'tract': tract, 'patch':patch, 'filter':filter}

    if calexp:
        dataType = "deepCoadd_calexp"
    else:
        dataType = "deepCoadd"

    # get the exposure from the butler
    # Ugly work around in case the before and after Reruns are from different hscPipe
    try:
        exposure = butler.get(dataType, dataId, immediate=True)
    except:
        try:
            exposure = butler.get('deepCoadd', dataId, immediate=True)
        except:
            raise

    # get the measurement catalog
    srcCat = butler.get('deepCoadd_meas', dataId, immediate=True)

    # convert to a numpy ndarray
    return exposure, srcCat


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


def srcToRaDec(src):

    """
    Return a list of (RA, DEC)
    """

    raArr  = np.asarray([coord[0] * 180.0 / np.pi for coord in src.get('coord')])
    decArr = np.asarray([coord[1] * 180.0 / np.pi for coord in src.get('coord')])

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


def main(root, tract, patch, filter, zeropoint=27.0, calexp=False, dev=False,
        prefix="patchCModel"):

    """
    Show the image of a Patch, and overplot the shape of cModels results
    """

    expPatch, srcPatch = getCalexpCat(root, tract, patch, filter, calexp=calexp)

    """
    Image Array
    """
    imgPatch = expPatch.getMaskedImage().getImage().getArray()
    imin, imax = zscale(imgPatch, contrast=0.10, samples=500)

    """
    Catalog
    """
    print "# There are %d sources measured!" % (len(srcPatch))

    print "# Get xUse, yUse"
    x0, y0 = expPatch.getXY0()
    xUse, yUse = (srcPatch.getX() - x0), (srcPatch.getY() - y0)

    print "# Convert the source to AstroPy table"
    tabPatch = getAstroTable(srcPatch, mags=True, zeropoint=zeropoint)
    #print tabPatch.colnames

    print "# Get the cModel magnitude"
    cmodelMag = tabPatch['cmodel.mag']
    cmodelColor = toColorArr(cmodelMag, top=23.5, bottom=18.5)
    print "# Get the color array for cModel magnitde: %d -- %d" % (np.nanmin(cmodelColor),
            np.nanmax(cmodelColor))

    print "# Convert the Re, q, PA into Ellipse"
    if not dev:
        rExp = tabPatch['cmodel.exp.ellipse_a']
        eExp = (1.0 - tabPatch['cmodel.exp.ellipse_q'])
        pExp = tabPatch['cmodel.exp.ellipse_theta']
        print "# Min/Max Theta : %7.3f - %7.3f" % (np.nanmin(pExp), np.nanmax(pExp))
        modEllip = getEll2Plot(xUse, yUse, rExp, eExp, pExp)
        modStr = 'Exp'
    else:
        rDev = tabPatch['cmodel.dev.ellipse_a']
        eDev = (1.0 - tabPatch['cmodel.dev.ellipse_q'])
        pDev = tabPatch['cmodel.dev.ellipse_theta']
        print "# Min/Max Theta : %7.3f - %7.3f" % (np.nanmin(pDev), np.nanmax(pDev))
        modEllip = getEll2Plot(xUse, yUse, rDev, eDev, pDev)
        modStr = 'Dev'

    """ Fig 1 """
    fig = plt.figure(figsize=(30, 30))
    fig.subplots_adjust(hspace=0.0, wspace=0.0,
                        left=0.03, bottom=0.03,
                        top=0.95, right=0.99)
    ax = fig.add_subplot(1,1,1)
    fontsize = 16
    ax.minorticks_on()

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    ax.set_title('%s - %s - %s Image - cModel/%s' % (tract, patch, filter, modStr),
            fontsize=38, fontweight='bold')
    ax.title.set_position((0.5, 1.01))

    # Grey scale image
    print "# Plot the Grey scale image"
    ax.imshow(np.arcsinh(imgPatch), interpolation="none",
           vmin=imin, vmax=imax,
           cmap=cmap1)

    # Scatter points for detections
    print "# Scatter plot the detections"
    ax.scatter(xUse, yUse, marker='+', s=22, c='r', alpha=0.3)

    # Ellipses for cModel
    print "# Overplot the ellipses for cModel"
    for (e, c) in zip(modEllip, cmodelColor):
        if np.isfinite(c):
            ax.add_artist(e)
            e.set_clip_box(ax.bbox)
            e.set_alpha(0.8)
            e.set_edgecolor(cmap2(int(c)))
            e.set_facecolor('none')
            e.set_linewidth(1.3)

    # Add the color bar
    print "# Overplot the color bar"
    cax = fig.add_axes([0.07, 0.09, 0.18, 0.025])
    norm = mpl.colors.Normalize(vmin=18.5, vmax=23.5)
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap2, norm=norm,
                                     orientation='horizontal')
    cbar.set_label('cModel Magnitude (mag)', fontsize=24)

    ax.set_xlim(0, imgPatch.shape[1]-1)
    ax.set_ylim(0, imgPatch.shape[0]-1)

    fig.savefig("%s_%s_%s-%s-%s.png" % (prefix, modStr, tract, patch, filter))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("root",   help="Root directory of data repository")
    parser.add_argument("tract", type=int, help="Tract to show")
    parser.add_argument("patch", help="Patch to show")
    parser.add_argument("filter", help="Filter to show")
    parser.add_argument("-p", "--prefix", dest="prefix",
            default='patchCModel', help="Prefix of the file name")
    parser.add_argument("-z", "--zeropoint", type=float, default=27.0,
            dest='zeropoint', help="Photometric zeropoint")
    parser.add_argument("-c", "--calexp", help="Whether to use the deepCoadd_calexp",
            dest='calexp', action='store_true', default=False)
    parser.add_argument("-d", "--dev", help="Whether to use the cModel.dev",
            dest='dev', action='store_true', default=False)
    args = parser.parse_args()

    main(args.root, args.tract, args.patch, args.filter,
            zeropoint=args.zeropoint, calexp=args.calexp, dev=args.dev,
            prefix=args.prefix)
