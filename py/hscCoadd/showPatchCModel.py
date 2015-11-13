from __future__ import division

import copy
import numpy as np

# Matplotlib default settings
rcdef = plt.rcParams.copy()
pylab.rcParams['figure.figsize'] = 12, 10
pylab.rcParams['xtick.major.size'] = 8.0
pylab.rcParams['xtick.major.width'] = 1.5
pylab.rcParams['xtick.minor.size'] = 4.0
pylab.rcParams['xtick.minor.width'] = 1.5
pylab.rcParams['ytick.major.size'] = 8.0
pylab.rcParams['ytick.major.width'] = 1.5
pylab.rcParams['ytick.minor.size'] = 4.0
pylab.rcParams['ytick.minor.width'] = 1.5
rc('axes', linewidth=2)

from matplotlib.patches import Ellipse

from astropy.io import fits
from astropy import units as u
from astropy.stats import sigma_clip
from astropy.wcs import WCS

import cubehelix  # Cubehelix color scheme from https://github.com/jradavenport/cubehelix

cmap1 = cubehelix.cmap(start=0.5, rot=-0.8, gamma=1.0,
                       minSat=1.2, maxSat=1.2,
                       minLight=0.0, maxLight=1.0)
cmap2 = cubehelix.cmap(start=2.0, rot=-1.0, gamma=2.5,
                       minSat=1.2, maxSat=1.2,
                       minLight=0.0, maxLight=1.0, reverse=True)
cmap3 = cubehelix.cmap(start=0.5, rot=-0.8, gamma=1.2,
                       minSat=1.2, maxSat=1.2,
                       minLight=0.0, maxLight=1.0)
cmap4 = cubehelix.cmap(start=0.5, rot=-0.8, gamma=0.7,
                       minSat=1.2, maxSat=1.2,
                       minLight=0.0, maxLight=1.0)

from palettable.colorbrewer.sequential import Greys_3 as pcmap1
cmap5 = pcmap1.mpl_colormap

from palettable.colorbrewer.diverging  import RdYlGn_11 as pcmap2
cmap6 = pcmap2.mpl_colormap

import lsst.daf.persistence as dafPersist
import lsst.afw.geom.ellipses as geomEllip
from lsst.afw.table import SourceCatalog, SchemaMapper

from distutils.version import StrictVersion


def getCalexpCat(root, tract, patch, filter):

    # make a butler and specify your dataId
    butler = dafPersist.Butler(root)
    dataId = {'tract': tract, 'patch':patch, 'filter':filter}

    pipeVersion = dafPersist.eupsVersions.EupsVersions().versions['hscPipe']
    if StrictVersion(pipeVersion) >= StrictVersion('3.9.0'):
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

    # get the maskedImage from the exposure, and the image from the mimg
    mimg = exposure.getMaskedImage()
    img = mimg.getImage()

    # get the measurement catalog
    srcCat = butler.get('deepCoadd_meas', dataId, immediate=True)

    # convert to a numpy ndarray
    return img.getArray(), srcCat


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


def main(root, tract, patch, filter):

    """
    Show the image of a Patch, and overplot the shape of cModels results
    """

    imgPatch, srcPatch = getExpArray(root, tract, patch, filter)

    imin, imax = zscale(imgPatch, contrast=0.10, samples=500)

    """ Fig 1 """
    fig = plt.figure(figsize=(20, 20))
    fig.subplots_adjust(hspace=0.1, wspace=0.1,
                        top=0.95, right=0.95)
    ax = gca()
    fontsize = 14
    ax.minorticks_on()

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    ax.set_title('i-band Image - cModel/Exp', fontsize=25, fontweight='bold')
    ax.title.set_position((0.5,1.01))

    ax.imshow(np.arcsinh(imgData), interpolation="none",
           vmin=imin, vmax=imax,
           cmap=cmap5)

    ax.set_xlim(0, imgData.shape[1]-1)
    ax.set_ylim(0, imgData.shape[0]-1)


    fig.savefig("patchCmodel_%s-%s-%s.png" % (tract, patch, filter))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("root",   help="Root directory of data repository")
    parser.add_argument("tract", type=int, help="Tract to show")
    parser.add_argument("patch", help="Patch to show")
    parser.add_argument("filter", help="Filter to show")
    args = parser.parse_args()

    main(args.root, args.tract, args.patch, args.filter)
