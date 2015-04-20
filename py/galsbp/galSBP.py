#!/usr/bin/env python
# encoding: utf-8

from __future__ import division

import os
import copy
import subprocess
import numpy as np

# Matplotlib default settings
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

# Astropy related
from astropy.io import fits
from astropy import units as u
from astropy.table import Table, Column

# Color table
try:
    import cubehelix  # Cubehelix color scheme
    cmap = cubehelix.cmap(start=0.5, rot=-1.0, minSat=1.2, maxSat=1.2,
            minLight=0., maxLight=1., gamma=0.5)
except ImportError:
    cmap = 'spectral'
cmap.set_bad('k',1.)

# STScI Python
from pyraf import iraf


def convIso2Ell(ellTab, xpad=0.0, ypad=0.0):
    """TODO: Docstring for convIso2Ell
    :returns: TODO

    """

    x  = ellTab['x0'] - xpad
    y  = ellTab['y0'] - ypad
    pa = ellTab['pa']
    a  = ellTab['sma'] * 2.0
    b  = ellTab['sma'] * 2.0 * (1.0 - ellTab['ell'])

    ells = [Ellipse(xy=np.array([x[i], y[i]]),
                    width=np.array(b[i]),
                    height=np.array(a[i]),
                    angle=np.array(pa[i]))
            for i in range(x.shape[0])]

    return ells


def maskFits2Pl(inputFits):

    """TODO: Docstring for maskFits2Pl.
    :returns: TODO

    """
    if not os.path.isfile(inputFits):
        raise Exception("Can not find the FITS mask: %s" % inputFits)

    # Name of the .pl mask file for IRAF
    mskIraf = inputFits.replace('.fits', '.pl')

    if os.path.isfile(mskIraf):
        os.remove(mskIraf)

    # Convert the fits format mask into pl format.
    iraf.imcopy.input  = "'" + inputFits + "'"
    iraf.imcopy.output = "'" + mskIraf + "'"
    iraf.imcopy()

    return mskIraf


def ellipseConfig(galX, galY, maxSMA, galEll=0.05, galPA=0.0,
                  iniSMA=6.0, minSMA=0.0, ellStep=0.1,
                  recenter=True, inEllip=None, harmonic=False):
    """TODO: Docstring for ellipseConfig
    :returns: TODO

    """

    # Define parameters for the ellipse run
    # 1. Initial guess of the central X, Y
    if (isinstance(galX, (float, int))) and (isinstance(galY, (float, int))):
        iraf.ellipse.geompar.x0 = galX
        iraf.ellipse.geompar.y0 = galY
    else:
        raise "Make sure that the input X0 and Y0 are meaningful !", galX, galY

    # 2. Initial guess of the ellipticity and PA of the first ISOPHOTE
    if (isinstance(galEll, (float, int))) and (isinstance(galPA, (float, int))):
        iraf.ellipse.geompar.ellip0 = galEll
        iraf.ellipse.geompar.pa0 = galPA
    else:
        raise "Make sure that the input Ellipticity and PA are meaningful !", galEll, galPA

    # 3. Initial radius for ellipse fitting
    iraf.ellipse.geompar.sma0       = iniSMA
    # 4. The minimum and maximum radius for the ellipse fitting
    iraf.ellipse.geompar.minsma     = minSMA
    iraf.ellipse.geompar.maxsma     = maxSMA
    # 5. Parameters about the stepsize during the fitting.
    iraf.ellipse.geompar.linear     = "no"
    iraf.ellipse.geompar.step       = ellStep
    # 6. Do you want to allow the ellipse to decide the galaxy center during the
    iraf.ellipse.geompar.recenter   = "yes"
    # 7. The next three parameters control the behavior of the fit
    iraf.ellipse.controlpar.conver  = 2
    iraf.ellipse.controlpar.hcenter = "no"
    iraf.ellipse.controlpar.hellip  = "no"
    iraf.ellipse.controlpar.hpa     = "no"
    # 8. Parameters about the iterations
    # minit/maxit: minimun and maximum number of the iterations
    iraf.ellipse.controlpar.minit   = 10
    iraf.ellipse.controlpar.maxit   = 200
    # 9. Threshold for the object locator algorithm
    iraf.ellipse.controlpar.olthresh = 1.00000
    # 10. Make sure the Interactive Mode is turned off
    iraf.ellipse.interactive         = "no"
    # 11. Magnitude Zeropoint
    iraf.ellipse.magpar.mag0         = zpPhoto
    # 12. Sampler
    iraf.ellipse.samplepar.integrmode  = intMode1  # "bi-linear" "mean" "median"
    iraf.ellipse.samplepar.usclip      = uppClip1
    iraf.ellipse.samplepar.lsclip      = lowClip1
    iraf.ellipse.samplepar.nclip       = nClip1
    iraf.ellipse.samplepar.fflag       = fracBad1
    if harmonic:
        # 13. Optional Harmonics
        iraf.ellipse.samplepar.harmonics   = '1,2'


def ellipRemoveIndef(outTabName):

    """TODO: Docstring for ellipRemoveIndef
    :returns: TODO

    """

    if os.path.exists(outTabName):
        subprocess.call(['sed', '-i_back', 's/INDEF/NaN/g', outTabName])
    else:
        raise Exception('Can not find the input catalog!')

    return outTabName


def readEllipseOut(outTabName, pix=1.0, zp=27.0, exptime=1.0, bkg=0.0,
                   harmonic=False):

    """TODO: Docstring for readEllipseOut
    :returns: TODO

    """
    name = ellipRemoveIndef(outTabName)

    ellipseOut = Table.read(outTabName, format='ascii.no_header')

    # Rename all the columns
    ellipseOut.rename_column('col1',  'sma')
    ellipseOut.rename_column('col2',  'intens')
    ellipseOut.rename_column('col3',  'int_err')
    ellipseOut.rename_column('col4',  'pix_var')
    ellipseOut.rename_column('col5',  'rms')
    ellipseOut.rename_column('col6',  'ell')
    ellipseOut.rename_column('col7',  'ell_err')
    ellipseOut.rename_column('col8',  'pa')
    ellipseOut.rename_column('col9',  'pa_err')
    ellipseOut.rename_column('col10', 'x0')
    ellipseOut.rename_column('col11', 'x0_err')
    ellipseOut.rename_column('col12', 'y0')
    ellipseOut.rename_column('col13', 'y0_err')
    ellipseOut.rename_column('col14', 'grad')
    ellipseOut.rename_column('col15', 'grad_err')
    ellipseOut.rename_column('col16', 'grad_r_err')
    ellipseOut.rename_column('col17', 'rsma')
    ellipseOut.rename_column('col18', 'mag')
    ellipseOut.rename_column('col19', 'mag_lerr')
    ellipseOut.rename_column('col20', 'mag_uerr')
    ellipseOut.rename_column('col21', 'tflux_e')
    ellipseOut.rename_column('col22', 'tflux_c')
    ellipseOut.rename_column('col23', 'tmag_e')
    ellipseOut.rename_column('col24', 'tmag_c')
    ellipseOut.rename_column('col25', 'npix_e')
    ellipseOut.rename_column('col26', 'npix_c')
    ellipseOut.rename_column('col27', 'a3')
    ellipseOut.rename_column('col28', 'a3_err')
    ellipseOut.rename_column('col29', 'b3')
    ellipseOut.rename_column('col30', 'b3_err')
    ellipseOut.rename_column('col31', 'a4')
    ellipseOut.rename_column('col32', 'a4_err')
    ellipseOut.rename_column('col33', 'b4')
    ellipseOut.rename_column('col34', 'b4_err')
    ellipseOut.rename_column('col35', 'ndata')
    ellipseOut.rename_column('col36', 'nflag')
    ellipseOut.rename_column('col37', 'niter')
    ellipseOut.rename_column('col38', 'stop')
    ellipseOut.rename_column('col39', 'a_big')
    ellipseOut.rename_column('col40', 'sarea')
    if harmonic:
        ellipseOut.rename_column('col41', 'a1')
        ellipseOut.rename_column('col42', 'a1_err')
        ellipseOut.rename_column('col43', 'a2')
        ellipseOut.rename_column('col44', 'a2_err')

    # Normalize the PA
    import angles
    ellipseOut.add_column(Column(name='pa_norm',
                          data=np.array([angles.normalize(pa, 0.0, 180.0)
                              for pa in ellipseOut['pa']])))

    # Convert the intensity into surface brightness
    parea = (pix ** 2.0)
    ellipseOut.add_column(Column(name='sbp',
                          data=(zp-2.5*np.log10((ellipseOut['intens']-bkg)/(parea*exptime)))))
    ellipseOut.add_column(Column(name='sbp_lerr',
                          data=(zp-2.5*np.log10((ellipseOut['intens']-ellipseOut['int_err']-bkg)/(parea*exptime)))))
    ellipseOut.add_column(Column(name='sbp_uerr',
                          data=(zp-2.5*np.log10((ellipseOut['intens']+ellipseOut['int_err']-bkg)/(parea*exptime)))))

    # Convert the unit of radius into arcsecs
    ellipseOut.add_column(Column(name='sma_asec',
                                 data=(ellipseOut['sma'] * pix)))
    ellipseOut.add_column(Column(name='rsma_asec',
                                 data=(ellipseOut['sma'] * pix) ** 0.25))

    return ellipseOut


def zscale(img, contrast=0.25, samples=500):

    """TODO: Docstring for zscale
    :returns: TODO

    """

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


def ellipseGetGrowthCurve(ellipOut, relThreshold=(5E-3)):

    """TODO: Docstring for ellipseGetGrowthCurve
    :returns: TODO

    """

    # The area in unit of pixels covered by an elliptical isophote
    ellArea = np.pi * (ellipOut['sma'] ** 2.0 * (1.0-ellipOut['ell']))
    # The area in unit covered by the "ring"
    isoArea = np.append(ellArea[0], [ellArea[1:] - ellArea[:-1]])
    # The total flux inside the "ring"
    isoFlux = np.append(ellArea[0], [ellArea[1:] - ellArea[:-1]]) * ellipOut['intens']
    #
    isoTFlux = np.asarray(map(lambda x: np.nansum(isoFlux[0:x+1]), range(isoFlux.shape[0])))
    # Flux difference between each "ring"
    diffIsoTFlux = np.append(isoTFlux[0], [isoTFlux[1:] - isoTFlux[:-1]])
    # Relative change of flux in each "ring"
    relDiffIso = (diffIsoTFlux / isoTFlux)

    # TODO: Using more examples to show whether this really works well

    indFluxDecrease = np.where(relDiffIso < relThreshold)
    if (indFluxDecrease[0]).size is 0:
        import warnings
        warnings.warn("WARNING!! the flux increasement is never smaller than the threshold!")
        indMaxIso = (relDiffIso.shape[0] - 1)
    else:
        indMaxIso = (np.where(diffIsoTFlux < 0))[0][0]

    maxIsoSma  = ellipOut['sma'][indMaxIso]

    maxIsoFlux = np.nanmax(isoTFlux[0:indMaxIso])

    # Get the growth curve
    isoGrowthCurve = np.asarray((isoTFlux / maxIsoFlux) * 100.0)

    return isoGrowthCurve, maxIsoSma


def ellipseGetR50(ellipseRsma, isoGrowthCurve, simple=True):

    if len(ellipseRsma) != len(isoGrowthCurve):
        raise "The x and y should have the same size!", len(ellipseRsma), len(isoGrowthCurve)
    else:
        if simple:
            isoRsma50 = ellipseRsma[np.nanargmin(np.abs(isoGrowthCurve - 50.0))]
        else:
            isoRsma50 = (numpy.interp([50.0], isoGrowthCurve, ellipseRsma))[0]

    return isoRsma50


def galEllipse(image, mask, step, inEllip=None, psf=None, galX=None, galY=None,
              galR=None, galQ=None, galPA=None, pixelScale=None,
              gain=None, expTime=None, zpPhoto=None,
              galInnerRatio=0.4, galOuterRatio=20.0, minSmaSys=6.0, maxSmaSysFac=1.2,
              ellStep=0.08, uppClip1=2.0, lowClip=2.0, nClip=3, fracBad=0.5,
              intModel="median", prefix=prefix, plMask=True):

    """TODO: Docstring for galEllipse.

    step  = 1: All Free
            2: Center Fixed
            3: All geometry fixd
            4: Force Photometry, must have inEllip
    :returns: TODO

    """
    # Call the STSDAS.ANALYSIS.ISOPHOTE package
    iraf.stsdas()
    iraf.analysis()
    iraf.isophote()

    """ Check input files """
