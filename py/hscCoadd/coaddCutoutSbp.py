#!/usr/bin/env python
# encoding: utf-8

from __future__ import division

import os
import gc
import copy
import argparse

try:
    import psutil
    psutilOk = True
except Exception:
    psutilOk =False

import numpy as np
# Astropy
from astropy.io import fits

# Cubehelix color scheme from https://github.com/jradavenport/cubehelix
import cubehelix
# For high-contrast image
cmap = cubehelix.cmap(start=0.5, rot=-0.8, gamma=1.0,
                      minSat=1.2, maxSat=1.2,
                      minLight=0.0, maxLight=1.0)
cmap.set_bad('k', 1.)

# Matplotlib related
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['figure.figsize'] = 12, 10
mpl.rcParams['xtick.major.size'] = 10.0
mpl.rcParams['xtick.major.width'] = 2.5
mpl.rcParams['xtick.minor.size'] = 5.0
mpl.rcParams['xtick.minor.width'] = 2.5
mpl.rcParams['ytick.major.size'] = 10.0
mpl.rcParams['ytick.major.width'] = 2.5
mpl.rcParams['ytick.minor.size'] = 5.0
mpl.rcParams['ytick.minor.width'] = 2.5
mpl.rc('axes', linewidth=3.0)
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import MaxNLocator

# Personal
import hscUtils as hUtil
import galSBP

COM = '#' * 100
SEP = '-' * 100
WAR = '!' * 100


def readSbpInput(prefix, root=None, exMask=None):
    """
    Read in necessary files for Ellipse run.

    Parameters:
    """
    # Get the names of necessary input images
    imgFile = prefix + '_img.fits'
    if root is not None:
        imgFile = os.path.join(root, imgFile)

    if exMask is None:
        mskFile = prefix + '_mskfin.fits'
        if root is not None:
            mskFile = os.path.join(root, mskFile)
    else:
        mskFile = exMask

    if os.path.islink(imgFile):
        imgOri = os.readlink(imgFile)
    else:
        imgOri = imgFile

    if os.path.islink(mskFile):
        mskOri = os.readlink(mskFile)
    else:
        mskOri = mskFile

    if not os.path.isfile(imgOri):
        raise Exception("### Can not find the input \
                cutout image : %s !" % imgOri)
    if not os.path.isfile(mskOri):
        raise Exception("### Can not find the input \
                mask image : %s !" % mskOri)
    # Image
    imgHdu = fits.open(imgOri)
    imgArr = imgHdu[0].data
    imgHead = imgHdu[0].header
    # Mask
    mskHdu = fits.open(mskOri)
    mskArr = mskHdu[0].data
    mskHead = mskHdu[0].header

    return imgFile, imgArr, imgHead, mskFile, mskArr, mskHead


def readPsfInput(prefix, root=None):
    """
    Read in the PSF file.

    Parameters:
    """
    psfFile = prefix + '_psf.fits'

    if root is not None:
        psfFile = os.path.join(root, psfFile)

    if os.path.islink(psfFile):
        psfOri = os.readlink(psfFile)

    if not os.path.isfile(psfOri):
        raise Exception("### Can not find the input psf image : %s !" % psfOri)

    # PSF
    psfHdu = fits.open(psfOri)
    psfArr = psfHdu[0].data
    psfDimX, psfDimY = psfArr.shape

    return psfFile, psfDimX, psfDimY


def readInputSky(prefix, root=None, rebin='rebin6'):
    """
    Read in the Input Sky Result.

    Parameters:
    """
    skyFile = prefix + '_' + rebin + '_sky.dat'
    if root is not None:
        skyFile = os.path.join(root, skyFile)

    if not os.path.isfile(skyFile):
        raise Exception("### Can not find the input \
                sky summary : %s !" % skyFile)

    skySum = open(skyFile, 'r').readlines()
    skyMed = float(skySum[3].split(':')[1].strip())
    skyAvg = float(skySum[4].split(':')[1].strip())
    skyStd = float(skySum[5].split(':')[1].strip())

    return skyMed, skyAvg, skyStd


def imgSameSize(img1, img2):
    """Check whether images have the same size."""
    dimX1, dimY1 = img1.shape
    dimX2, dimY2 = img2.shape
    if (dimX1 == dimX2) and (dimY1 == dimY2):
        return True
    else:
        return False


def ellipSummary(ellipOut1, ellipOut2, ellipOut3, image,
                 maxRad=None, mask=None, radMode='rsma',
                 outPng='ellipse_summary.png', zp=27.0, threshold=None,
                 psfOut=None, useKpc=None, pix=0.168,
                 showZoom=True, exptime=1.0, bkg=0.0, outRatio=1.2):
    """
    Make a summary plot for the ellipse run.

    Parameters:
    """
    """ Left side: SBP """
    reg1 = [0.07, 0.05, 0.455, 0.35]
    reg2 = [0.07, 0.40, 0.455, 0.15]
    reg3 = [0.07, 0.55, 0.455, 0.15]
    reg4 = [0.07, 0.70, 0.455, 0.15]
    reg5 = [0.07, 0.85, 0.455, 0.14]

    """ Right side: Curve of growth & IsoMap """
    reg6 = [0.59, 0.05, 0.39, 0.30]
    reg7 = [0.59, 0.35, 0.39, 0.16]
    reg8 = [0.59, 0.55, 0.39, 0.39]

    fig = plt.figure(figsize=(20, 20))
    """ Left """
    ax1 = fig.add_axes(reg1)
    ax2 = fig.add_axes(reg2)
    ax3 = fig.add_axes(reg3)
    ax4 = fig.add_axes(reg4)
    ax5 = fig.add_axes(reg5)
    """ Right """
    ax6 = fig.add_axes(reg6)
    ax7 = fig.add_axes(reg7)
    ax8 = fig.add_axes(reg8)

    """ Image """
    img = fits.open(image)[0].data
    imgX, imgY = img.shape
    imgMsk = copy.deepcopy(img)
    imin, imax = hUtil.zscale(imgMsk, contrast=0.6, samples=500)
    if mask is not None:
        msk = fits.open(mask)[0].data
        imgMsk[msk > 0] = np.nan

    """ Find the proper outer boundary """
    sma = ellipOut3['sma']
    radOuter = galSBP.ellipseGetOuterBoundary(ellipOut3, ratio=outRatio,
                                              threshold=threshold,
                                              polyOrder=12)

    if not np.isfinite(radOuter):
        print " XXX radOuter is NaN, use 0.80 * max(SMA) instead !"
        radOuter = np.nanmax(sma) * 0.80

    indexUse1 = np.where(ellipOut1['sma'] <= (radOuter*1.2))
    indexUse2 = np.where(ellipOut2['sma'] <= (radOuter*1.2))
    indexUse3 = np.where(ellipOut3['sma'] <= (radOuter*1.2))
    print "###     OutRadius", radOuter
    curveOri = ellipOut3['growth_ori']
    curveOri[curveOri < 0.0] = 0.0
    curveCor = ellipOut3['growth_cor']
    curveCor[curveCor < 0.0] = 0.0
    growthCurveOri = -2.5 * np.log10(curveOri) + zp
    growthCurveNew = -2.5 * np.log10(curveCor) + zp

    maxIsoFluxOri = np.nanmax(ellipOut3['growth_ori'][indexUse3])
    # TODO: magFluxOri50  = -2.5 * np.log10(maxIsoFluxOri * 0.50) + zp
    magFluxOri100 = -2.5 * np.log10(maxIsoFluxOri) + zp
    print "###     MagTot OLD : ", magFluxOri100
    ax1.text(0.6, 0.85, 'mag$_{tot,old}=%5.2f$' % magFluxOri100, fontsize=24,
             transform=ax1.transAxes)

    maxIsoFlux = np.nanmax(ellipOut3['growth_cor'][indexUse3])
    magFlux50 = -2.5 * np.log10(maxIsoFlux * 0.50) + zp
    magFlux100 = -2.5 * np.log10(maxIsoFlux) + zp
    print "###     MagTot NEW : ", magFlux100
    ax1.text(0.6, 0.78, 'mag$_{tot,new}=%5.2f$' % magFlux100, fontsize=24,
             transform=ax1.transAxes)

    indMaxFlux = np.nanargmax(ellipOut3['growth_cor'][indexUse3])
    maxIsoSbp = ellipOut3['sbp_upp'][indMaxFlux]
    print "###     MaxIsoSbp : ", maxIsoSbp

    """ Type of Radius """
    if radMode is 'rsma':
        if useKpc is None:
            radStr = 'RSMA (arcsec$^{1/4}$)'
            rad1 = ellipOut1['rsma_asec']
            rad2 = ellipOut2['rsma_asec']
            rad3 = ellipOut3['rsma_asec']
        else:
            radStr = 'RSMA (kpc$^{1/4}$)'
            rad1 = (ellipOut1['sma_asec']*useKpc)**0.25
            rad2 = (ellipOut2['sma_asec']*useKpc)**0.25
            rad3 = (ellipOut3['sma_asec']*useKpc)**0.25

        minRad = 0.41 if 0.41 >= np.nanmin(rad3) else np.nanmin(rad3)
        if useKpc is None:
            imgR50 = (imgX * pix / 2.0) ** 0.25
            radOut = (radOuter * pix * 1.2) ** 0.25
        else:
            imgR50 = (imgX * pix * useKpc / 2.0) ** 0.25
            radOut = (radOuter * pix * useKpc * 1.2) ** 0.25
        if maxRad is None:
            maxSma = np.nanmax(ellipOut3['sma'])
            maxRad = np.nanmax(rad3)
        else:
            maxSma = maxRad
            if useKpc is None:
                maxRad = (maxRad * pix) ** 0.25
            else:
                maxRad = (maxRad * pix * useKpc) ** 0.25
    elif radMode is 'sma':
        if useKpc is None:
            radStr = 'SMA (arcsec)'
            rad1 = ellipOut1['sma_asec']
            rad2 = ellipOut2['sma_asec']
            rad3 = ellipOut3['sma_asec']
        else:
            radStr = 'SMA (kpc)'
            rad1 = ellipOut1['sma_asec'] * useKpc
            rad2 = ellipOut2['sma_asec'] * useKpc
            rad3 = ellipOut3['sma_asec'] * useKpc

        minRad = 0.05 if 0.05 >= np.nanmin(rad3) else np.nanmin(rad3)
        if useKpc is None:
            imgR50 = (imgX * pix / 2.0)
            radOut = radOuter * pix * 1.2
        else:
            imgR50 = (imgX * pix * useKpc / 2.0)
            radOut = radOuter * pix * useKpc * 1.2
        if maxRad is None:
            maxSma = np.nanmax(ellipOut3['sma'])
            maxRad = np.nanmax(rad3)
        else:
            maxSma = maxRad
            if useKpc is None:
                maxRad = maxRad * pix
            else:
                maxRad = maxRad * pix * useKpc
    elif radMode is 'log':
        if useKpc is None:
            radStr = 'log (SMA/arcsec)'
            rad1 = np.log10(ellipOut1['sma_asec'])
            rad2 = np.log10(ellipOut2['sma_asec'])
            rad3 = np.log10(ellipOut3['sma_asec'])
        else:
            radStr = 'log (SMA/kpc)'
            rad1 = np.log10(ellipOut1['sma_asec']*useKpc)
            rad2 = np.log10(ellipOut2['sma_asec']*useKpc)
            rad3 = np.log10(ellipOut3['sma_asec']*useKpc)

        minRad = -1.2 if -1.2 >= np.nanmin(rad3) else np.nanmin(rad3)
        if useKpc is None:
            imgR50 = np.log10(imgX * pix / 2.0)
            radOut = np.log10(radOuter * pix * 1.2)
        else:
            imgR50 = np.log10(imgX * pix * useKpc / 2.0)
            radOut = np.log10(radOuter * pix * useKpc * 1.2)
        if maxRad is None:
            maxSma = np.nanmax(ellipOut3['sma'])
            maxRad = np.nanmax(rad3)
        else:
            maxSma = maxRad
            if useKpc is None:
                maxRad = np.log10(maxRad * pix)
            else:
                maxRad = np.log10(maxRad * pix * useKpc)
    else:
        raise Exception('### Wrong type of Radius: sma, rsma, log')

    """ ax1 SBP """
    ax1.minorticks_on()
    ax1.invert_yaxis()
    ax1.tick_params(axis='both', which='major', labelsize=22, pad=8)

    ax1.set_xlabel(radStr, fontsize=23)
    ax1.set_ylabel('${\mu}$ (mag/arcsec$^2$)', fontsize=28)

    ax1.plot(rad3[indexUse3], ellipOut3['sbp'][indexUse3], '--', color='k',
             linewidth=4.0)

    parea = (pix ** 2.0)
    sbp = zp - 2.5 * np.log10((ellipOut3['intens'] -
                              bkg) / (parea * exptime))
    sbp_low = zp - 2.5 * np.log10((ellipOut3['intens'] +
                                  ellipOut3['int_err'] - bkg) /
                                  (parea * exptime))
    sbp_err = (sbp - sbp_low)
    sbp_upp = (sbp + sbp_err)

    ax1.fill_between(rad3[indexUse3], sbp_upp[indexUse3], sbp_low[indexUse3],
                     facecolor='r', alpha=0.2)
    ax1.plot(rad3[indexUse3], sbp[indexUse3], '-', color='r', linewidth=3.5)

    ax1.set_xlim(minRad, radOut)
    sbpBuffer = 0.5
    minSbp = np.nanmin(ellipOut3['sbp_low'][indexUse3]) - sbpBuffer
    maxSbp = maxIsoSbp + 1.1
    maxSbp = maxSbp if maxSbp >= 29.0 else 28.9
    maxSbp = maxSbp if maxSbp <= 32.0 else 31.9

    if psfOut is not None:
        if radMode is 'rsma':
            psfRad = psfOut['rsma_asec']
        elif radMode is 'sma':
            psfRad = psfOut['sma_asec']
        elif radMode is 'log':
            psfRad = np.log10(psfOut['sma_asec'])
        else:
            raise Exception('### Wrong type of Radius: sma, rsma, log')
        psfSbp = psfOut['sbp'] - (np.nanmin(psfOut['sbp']) -
                                  np.nanmin(ellipOut3['sbp']))
        ax1.plot(psfRad, psfSbp, '-', color='k', linewidth=3.0, alpha=0.7)

    ax1.set_ylim(maxSbp, minSbp)

    """ ax2 Ellipticity """
    ax2.minorticks_on()
    ax2.tick_params(axis='both', which='major', labelsize=20, pad=8)
    ax2.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax2.yaxis.set_major_locator(MaxNLocator(prune='upper'))
    ax2.locator_params(axis='y', tight=True, nbins=4)

    ax2.set_ylabel('$e$', fontsize=30)

    print "###     AvgEll", (1.0 - ellipOut2['avg_q'][0])
    ax2.axhline((1.0 - ellipOut2['avg_q'][0]), color='k', linestyle='--',
                linewidth=3.0)
    ax2.fill_between(rad2[indexUse2],
                     (ellipOut2['ell'][indexUse2] +
                      ellipOut2['ell_err'][indexUse2]),
                     (ellipOut2['ell'][indexUse2] -
                      ellipOut2['ell_err'][indexUse2]),
                     facecolor='r', alpha=0.2)
    ax2.plot(rad2[indexUse2], ellipOut2['ell'][indexUse2], '-', color='r',
             linewidth=3.0)

    ax2.xaxis.set_major_formatter(NullFormatter())
    ax2.set_xlim(minRad, radOut)
    ellBuffer = 0.06
    minEll = np.nanmin(ellipOut2['ell'][indexUse2]) - ellBuffer
    maxEll = np.nanmax(ellipOut2['ell'][indexUse2]) + ellBuffer
    ax2.set_ylim(minEll, maxEll)

    """ ax3 PA """
    ax3.minorticks_on()
    ax3.tick_params(axis='both', which='major', labelsize=20, pad=8)
    ax3.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax3.yaxis.set_major_locator(MaxNLocator(prune='upper'))
    ax3.locator_params(axis='y', tight=True, nbins=4)

    ax3.set_ylabel('PA (degree)',  fontsize=23)

    medPA = np.nanmedian(ellipOut2['pa_norm'][indexUse2])
    avgPA = ellipOut2['avg_pa'][0]
    if (avgPA-medPA >= 85.0) and (avgPA <= 92.0):
        avgPA -= 180.0
    elif (avgPA-medPA <= -85.0) and (avgPA >= -92.0):
        avgPA += 180.0

    print "###     AvgPA", avgPA
    ax3.axhline(avgPA, color='k', linestyle='--', linewidth=3.0)

    ax3.fill_between(rad2[indexUse2],
                     (ellipOut2['pa_norm'][indexUse2] +
                      ellipOut2['pa_err'][indexUse2]),
                     (ellipOut2['pa_norm'][indexUse2] -
                      ellipOut2['pa_err'][indexUse2]),
                     facecolor='r', alpha=0.2)
    ax3.plot(rad2[indexUse2], ellipOut2['pa_norm'][indexUse2], '-', color='r',
             linewidth=3.0)

    ax3.xaxis.set_major_formatter(NullFormatter())
    ax3.set_xlim(minRad, radOut)

    paBuffer = 10.0
    minPA = np.nanmin(ellipOut2['pa_norm'][indexUse2])-paBuffer
    maxPA = np.nanmax(ellipOut2['pa_norm'][indexUse2])+paBuffer
    minPA = minPA if minPA >= -110.0 else -100.0
    maxPA = maxPA if maxPA <= 110.0 else 100.0
    ax3.set_ylim(minPA, maxPA)

    """ ax4 X0/Y0 """
    ax4.minorticks_on()
    ax4.tick_params(axis='both', which='major', labelsize=20, pad=8)
    ax4.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax4.yaxis.set_major_locator(MaxNLocator(prune='upper'))
    ax4.locator_params(axis='y', tight=True, nbins=4)

    ax4.set_ylabel('X0 or Y0 (pixel)', fontsize=23)

    print "###     AvgX0", ellipOut1['avg_x0'][0]
    print "###     AvgY0", ellipOut1['avg_y0'][0]
    ax4.axhline(ellipOut1['avg_x0'][0], linestyle='--', color='r', alpha=0.6,
                linewidth=3.0)
    ax4.fill_between(rad1[indexUse1],
                     (ellipOut1['x0'][indexUse1] +
                      ellipOut1['x0_err'][indexUse1]),
                     (ellipOut1['x0'][indexUse1] -
                      ellipOut1['x0_err'][indexUse1]),
                     facecolor='r', alpha=0.20)
    ax4.plot(rad1[indexUse1], ellipOut1['x0'][indexUse1], '-', color='r',
             linewidth=3.0, label='X0')

    ax4.axhline(ellipOut1['avg_y0'][0], linestyle='-.', color='b', alpha=0.6,
                linewidth=3.0)
    ax4.fill_between(rad1[indexUse1],
                     (ellipOut1['y0'][indexUse1] +
                     ellipOut1['y0_err'][indexUse1]),
                     (ellipOut1['y0'][indexUse1] -
                     ellipOut1['y0_err'][indexUse1]),
                     facecolor='b', alpha=0.25)
    ax4.plot(rad1[indexUse1], ellipOut1['y0'][indexUse1], '-', color='b',
             linewidth=3.0, label='Y0')

    ax4.xaxis.set_major_formatter(NullFormatter())
    ax4.set_xlim(minRad, radOut)
    xBuffer = 3.0
    minX0 = np.nanmin(ellipOut1['x0'][indexUse1])
    maxX0 = np.nanmax(ellipOut1['x0'][indexUse1])
    minY0 = np.nanmin(ellipOut1['y0'][indexUse1])
    maxY0 = np.nanmax(ellipOut1['y0'][indexUse1])
    minCen = minX0 if minX0 <= minY0 else minY0
    maxCen = maxX0 if maxX0 >= maxY0 else maxY0
    ax4.set_ylim(minCen-xBuffer, maxCen+xBuffer)

    """ ax5 A4/B4 """
    ax5.minorticks_on()
    ax5.tick_params(axis='both', which='major', labelsize=20, pad=8)
    ax5.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax5.yaxis.set_major_locator(MaxNLocator(prune='upper'))
    ax5.locator_params(axis='y', tight=True, nbins=4)

    ax5.set_ylabel('A4 or B4',  fontsize=23)

    ax5.axhline(0.0, linestyle='-', color='k', alpha=0.4)
    ax5.fill_between(rad2[indexUse2],
                     (ellipOut2['a4'][indexUse2] +
                     ellipOut2['a4_err'][indexUse2]),
                     (ellipOut2['a4'][indexUse2] -
                     ellipOut2['a4_err'][indexUse2]),
                     facecolor='r', alpha=0.20)
    ax5.plot(rad2[indexUse2], ellipOut2['a4'][indexUse2], '-',
             color='r', linewidth=3.0, label='A4')

    ax5.fill_between(rad2[indexUse2],
                     (ellipOut2['b4'][indexUse2] +
                     ellipOut2['b4_err'][indexUse2]),
                     (ellipOut2['b4'][indexUse2] -
                     ellipOut2['b4_err'][indexUse2]),
                     facecolor='b', alpha=0.20)
    ax5.plot(rad2[indexUse2], ellipOut2['b4'][indexUse2], '-',
             color='b', linewidth=3.0, label='B4')

    ax5.xaxis.set_major_formatter(NullFormatter())
    ax5.set_xlim(minRad, radOut)

    abBuffer = 0.02
    minA4 = np.nanmin(ellipOut2['a4'][indexUse2])
    minB4 = np.nanmin(ellipOut2['b4'][indexUse2])
    maxA4 = np.nanmax(ellipOut2['a4'][indexUse2])
    maxB4 = np.nanmax(ellipOut2['b4'][indexUse2])
    minAB = minA4 if minA4 <= minB4 else minB4
    maxAB = maxA4 if maxA4 >= maxB4 else maxB4
    ax5.set_ylim(minAB-abBuffer, maxAB+abBuffer)

    """ ax6 Growth Curve """
    ax6.minorticks_on()
    ax6.tick_params(axis='both', which='major', labelsize=22, pad=8)
    ax6.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax6.yaxis.set_major_locator(MaxNLocator(prune='upper'))

    ax6.set_xlabel(radStr, fontsize=23)
    ax6.set_ylabel('Curve of Growth (mag)', fontsize=24)

    ax6.axhline(magFlux100, linestyle='-', color='k',
                alpha=0.6, linewidth=3.0, label='mag$_{100}$')
    ax6.axhline(magFlux50,  linestyle='--', color='k',
                alpha=0.6, linewidth=3.0, label='mag$_{50}$')
    ax6.axvline(imgR50, linestyle='-', color='g', alpha=0.4, linewidth=3.0)
    ax6.axvline(radOut, linestyle='--', color='b', alpha=0.8, linewidth=3.0)

    ax6.plot(rad3, growthCurveOri, '--', color='g', linewidth=3.5,
             label='curve$_{old}$')
    ax6.plot(rad3, growthCurveNew, '-', color='r', linewidth=3.5,
             label='curve$_{new}$')

    ax6.legend(loc=[0.34, 0.59], fontsize=24)

    ax6.set_xlim(minRad, maxRad)

    """ ax7 Intensity Curve """
    ax7.minorticks_on()
    ax7.tick_params(axis='both', which='major', labelsize=22, pad=10)
    ax7.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax7.yaxis.set_major_locator(MaxNLocator(prune='upper'))
    ax7.locator_params(axis='y', tight=True, nbins=4)

    ax7.axhline(0.0, linestyle='-', color='k', alpha=0.6, linewidth=3.0)
    ax7.plot(rad3, ellipOut3['intens']+ellipOut3['avg_bkg'], '--', color='g',
             linewidth=2.5)
    ax7.fill_between(rad3, ellipOut3['intens']+ellipOut3['int_err'],
                     ellipOut3['intens']-ellipOut3['int_err'],
                     facecolor='r', alpha=0.3)
    ax7.plot(rad3, ellipOut3['intens'], '-', color='r', linewidth=3.0)
    ax7.axvline(imgR50, linestyle='-', color='g', alpha=0.4, linewidth=3.0)
    ax7.axvline(radOut, linestyle='--', color='b', alpha=0.8, linewidth=3.0)

    indexOut = np.where(ellipOut3['intens'] <= (0.002 *
                        np.nanmax(ellipOut3['intens'])))

    ax7.xaxis.set_major_formatter(NullFormatter())
    ax7.set_xlim(minRad, maxRad)
    minOut = np.nanmin(ellipOut3['intens'][indexOut] -
                       ellipOut3['int_err'][indexOut])
    maxOut = np.nanmax(ellipOut3['intens'][indexOut] +
                       ellipOut3['int_err'][indexOut])
    sepOut = (maxOut - minOut) / 10.0
    ax7.set_ylim(minOut - sepOut, maxOut)

    """ ax8 IsoPlot """

    imgFile = os.path.basename(image)
    imgTitle = imgFile.replace('.fits', '')
    imgTitle = imgTitle.replace('_img', '')

    ax8.tick_params(axis='both', which='major', labelsize=20)
    ax8.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax8.yaxis.set_major_locator(MaxNLocator(prune='upper'))
    ax8.set_title(imgTitle, fontsize=25, fontweight='bold')
    ax8.title.set_position((0.5, 1.05))

    galX0 = ellipOut3['avg_x0'][0]
    galY0 = ellipOut3['avg_y0'][0]
    imgSizeX, imgSizeY = img.shape

    if (galX0 > maxSma) and (galY0 > maxSma) and showZoom:
        zoomReg = imgMsk[np.int(galX0-maxSma):np.int(galX0+maxSma),
                         np.int(galY0-maxSma):np.int(galY0+maxSma)]
        # Define the new center of the cropped images
        xPad = (imgSizeX / 2.0 - maxSma)
        yPad = (imgSizeY / 2.0 - maxSma)
    else:
        zoomReg = imgMsk
        xPad = 0
        yPad = 0

    # Show the image
    ax8.imshow(np.arcsinh(zoomReg), interpolation="none",
               vmin=imin, vmax=imax, cmap=cmap, origin='lower')
    # Get the Shapes
    ellipIso = galSBP.convIso2Ell(ellipOut3, xpad=xPad, ypad=yPad)

    # Overlay the ellipses on the image
    for e in ellipIso:
        ax8.add_artist(e)
        e.set_clip_box(ax8.bbox)
        e.set_alpha(0.85)
        e.set_edgecolor('r')
        e.set_facecolor('none')
        e.set_linewidth(1.8)

    ellOuter = Ellipse(xy=(ellipOut3['avg_x0'][0],
                       ellipOut3['avg_y0'][0]),
                       height=(radOuter * 2.4),
                       width=(radOuter * 2.4 * ellipOut3['avg_q'][0]),
                       angle=(ellipOut3['avg_pa'][0]))
    ax8.add_artist(ellOuter)
    ellOuter.set_clip_box(ax8.bbox)
    ellOuter.set_alpha(0.9)
    ellOuter.set_edgecolor('b')
    ellOuter.set_facecolor('none')
    ellOuter.set_linewidth(3.5)

    """ Save Figure """
    fig.savefig(outPng)
    plt.close(fig)

    return


def coaddCutoutSbp(prefix, root=None, verbose=True, psf=True, inEllip=None,
                   zp=27.0, step=0.12, pix=0.168, exptime=1.0, bkgCor=False,
                   plot=True, galX0=None, galY0=None, galQ0=None, galPA0=None,
                   maxTry=4, galRe=None, redshift=None, psfRecenter=True,
                   showZoom=True, checkCenter=True, updateIntens=False,
                   olthresh=0.5, intMode='mean', lowClip=3.0, uppClip=3.0,
                   nClip=2, fracBad=0.5, minIt=20, maxIt=150, outRatio=1.2,
                   exMask=None, suffix='', plMask=False, noMask=False):
    """
    Generate 1-D SBP Plot.

    Parameters:
    """
    if verbose:
        print SEP
        print "#  CoaddCutoutSbp Start "
        print "##       Input Image: ", prefix
        print "##       Root Directory: ", root
    print SEP
    if psutilOk:
        proc = psutil.Process(os.getpid())
        gc.collect()
        mem0 = proc.memory_info().rss
        print "@@@ Initial: %i" % mem0
    else:
        gc.collect()
    print SEP
    """ 0. Organize Input Data """
    # Read in the input image, mask, psf, and their headers
    sbpInput = readSbpInput(prefix, root=root, exMask=exMask)
    imgFile, imgArr, imgHead, mskFile, mskArr, mskHead = sbpInput
    if (root is not None) and (root[-1] != '/'):
        root += '/'
    if not imgSameSize(imgArr, mskArr):
        raise Exception("### The Image and Mask need to have EXACTLY \
                same dimensions!")
    """
    Delete image and mask array

    TODO Test
    """
    if psutilOk:
        print WAR
        memB = proc.memory_info().rss
    del imgArr
    del mskArr
    if psutilOk:
        memA = proc.memory_info().rss
        print "@@@ Reduce : %5.2f" % ((memA - memB) / memB)
        print WAR
    """"""

    if os.path.islink(imgFile):
        imgOri = os.readlink(imgFile)
    else:
        imgOri = imgFile

    if os.path.islink(mskFile):
        mskOri = os.readlink(mskFile)
    else:
        mskOri = mskFile

    if noMask:
        mskFile = None
        mskOri = None

    """ 0a. Redshift """
    if redshift is not None:
        scale = hUtil.cosmoScale(redshift)
        useKpc = scale
        if verbose:
            print "###      Input Redshift : ", redshift
            print "###      Physical Scale : ", scale, " (kpc/arcsec)"
    else:
        useKpc = None

    """ 0b. Background """
    if bkgCor:
        try:
            skyMed, skyAvg, skyStd = readInputSky(prefix, root=root)
            bkg = skyAvg
            if verbose:
                print "###      Average Background : ", bkg
        except Exception:
            print WAR
            print "XXX   CAN NOT FIND THE BACKGROUND DATA !"
            print WAR
            bkg = 0.00
    else:
        bkg = 0.00

    """ 1. Prepare the Input for SBP """
    if (galX0 is None) or (galY0 is None):
        galX, galY = mskHead['GAL_CENX'], mskHead['GAL_CENY']
    else:
        galX, galY = galX0, galY0
    if (galQ0 is None) or (galPA0 is None):
        galQ, galPA = mskHead['GAL_Q'], mskHead['GAL_PA']
    else:
        galQ, galPA = galQ0, galPA0
    galQ = galQ if galQ <= 0.95 else 0.95
    galPA = hUtil.normAngle(galPA, lower=-90.0, upper=90.0, b=True)
    if galRe is None:
        galR50 = mskHead['GAL_R50']
    else:
        galR50 = galRe
    if verbose:
        print "###      Image : ", imgFile
        print "###      Mask  : ", mskFile
        print "###      galX, galY : ", galX, galY
        print "###      galQ, galPA : ", galQ, galPA
        print "###      galR50 : ", galR50
        print "###      intMode : ", intMode

    maxR = mskHead['NAXIS1'] if (mskHead['NAXIS1'] >=
                                 mskHead['NAXIS2']) else mskHead['NAXIS2']
    maxR = (maxR / 2.0) * np.sqrt(2.0)
    if verbose:
        print "###      maxR : ", maxR

    """
    Actually start to run
    """
    if psutilOk:
        print SEP
        mem1 = proc.memory_info().rss
        print "@@@ Increase: %0.2f%%" % (100.0 * (mem1 - mem0) / mem0)
        print SEP

    if checkCenter and mskHead['MSK_R20'] == 1:
        print WAR
        print "### The central region is masked out"
        print WAR
    else:
        if psf:
            print SEP
            print "\n##   Ellipse Run on PSF "
            print SEP
            psfFile = prefix + '_psf.fits'
            if root is not None:
                psfFile = os.path.join(root, psfFile)
            if os.path.islink(psfFile):
                psfOri = os.readlink(psfFile)
            else:
                psfOri = psfFile

            if not os.path.isfile(psfOri):
                raise Exception("### Can not find the \
                                PSF image: %s !" % psfFile)
            psfOut = galSBP.galSBP(psfFile, iniSma=5.0,
                                   pix=pix,
                                   galQ=0.95,
                                   galPA=0.0,
                                   stage=3,
                                   zpPhoto=zp,
                                   recenter=psfRecenter,
                                   outerThreshold=1e-6)
        if inEllip is None:
            iniSma = (galR50 * 2.0)
            """#        Start with Stage 1 """
            print SEP
            print "##       Ellipse Run on Image %s- Stage 1 " % imgFile
            print SEP
            galSBP.unlearnEllipse()
            ellOut1 = galSBP.galSBP(imgFile,
                                    mask=mskFile,
                                    galX=galX,
                                    galY=galY,
                                    maxSma=maxR,
                                    iniSma=iniSma,
                                    maxTry=maxTry,
                                    galR=galR50,
                                    ellipStep=step,
                                    pix=pix,
                                    bkg=bkg,
                                    galQ=galQ,
                                    galPA=galPA,
                                    stage=1,
                                    zpPhoto=zp,
                                    updateIntens=updateIntens,
                                    olthresh=olthresh,
                                    lowClip=lowClip,
                                    uppClip=uppClip,
                                    nClip=nClip,
                                    fracBad=fracBad,
                                    intMode=intMode,
                                    minIt=minIt,
                                    maxIt=maxIt,
                                    plMask=plMask,
                                    suffix=suffix)
            if ellOut1 is None:
                print WAR
                raise Exception("!!!!! ELLIPSE RUN FAILED AT STAGE 1 !!!!")
                print WAR
            if (galX0 is None) or (galY0 is None):
                galX0 = ellOut1['avg_x0'][0]
                galY0 = ellOut1['avg_y0'][0]

            if checkCenter:
                if (np.abs(galX0 -
                           mskHead['NAXIS1']/2.0) >= 20.0) or \
                   (np.abs(galX0 -
                           mskHead['NAXIS1']/2.0) >= 20.0):
                    raise Exception("The Center is Off !")

            """ # Start with Stage 2 """
            print SEP
            print "##       Ellipse Run on Image %s- Stage 2 " % imgFile
            print SEP
            ellOut2 = galSBP.galSBP(imgFile, mask=mskFile,
                                    galX=galX0,
                                    galY=galY0,
                                    maxSma=maxR,
                                    iniSma=iniSma,
                                    galR=galR50,
                                    ellipStep=step,
                                    maxTry=maxTry,
                                    pix=pix,
                                    bkg=bkg,
                                    galQ=galQ,
                                    galPA=galPA,
                                    stage=2,
                                    zpPhoto=zp,
                                    updateIntens=updateIntens,
                                    olthresh=olthresh,
                                    lowClip=lowClip,
                                    uppClip=uppClip,
                                    nClip=nClip,
                                    fracBad=fracBad,
                                    intMode=intMode,
                                    minIt=minIt,
                                    maxIt=maxIt,
                                    plMask=plMask,
                                    suffix=suffix)
            if ellOut2 is None:
                print WAR
                raise Exception("!!!!! ELLIPSE RUN FAILED AT STAGE 2 !!!!")
                print WAR

            if (galQ0 is None) or (galPA0 is None):
                if (ellOut2['avg_q'][0] <= 0.95):
                    galQ0 = ellOut2['avg_q'][0]
                else:
                    galQ0 = 0.95
                galPA0 = hUtil.normAngle(ellOut2['avg_pa'][0], lower=-90.0,
                                         upper=90.0, b=True)

            """ # Start with Stage 3 """
            print SEP
            print "##       Ellipse Run on Image %s- Stage 3 " % imgFile
            print SEP
            ellOut3 = galSBP.galSBP(imgFile, mask=mskFile,
                                    galX=galX0,
                                    galY=galY0,
                                    maxSma=maxR,
                                    iniSma=iniSma,
                                    galR=galR50,
                                    ellipStep=step,
                                    maxTry=maxTry,
                                    pix=pix,
                                    bkg=bkg,
                                    galQ=galQ0,
                                    galPA=galPA0,
                                    stage=3,
                                    zpPhoto=zp,
                                    updateIntens=updateIntens,
                                    olthresh=olthresh,
                                    lowClip=lowClip,
                                    uppClip=uppClip,
                                    nClip=nClip,
                                    fracBad=fracBad,
                                    intMode=intMode,
                                    plMask=plMask,
                                    suffix=suffix)
            if ellOut3 is None:
                print WAR
                raise Exception("!!!!! ELLIPSE RUN FAILED AT STAGE 3 !!!!")
                print WAR

            if plot:
                print SEP
                print "##       Ellipse Summary Plot "
                print SEP
                if suffix[-1] != '_':
                    suffix = suffix + '_'
                sumPng = root + prefix + '_ellip_' + suffix + 'sum.png'
                if psf:
                    ellipSummary(ellOut1, ellOut2, ellOut3, imgOri,
                                 psfOut=psfOut,
                                 maxRad=maxR, mask=mskOri, radMode='rsma',
                                 outPng=sumPng, zp=zp, useKpc=useKpc, pix=pix,
                                 showZoom=showZoom, exptime=exptime, bkg=bkg,
                                 outRatio=outRatio)
                else:
                    ellipSummary(ellOut1, ellOut2, ellOut3, imgOri,
                                 psfOut=None,
                                 maxRad=maxR, mask=mskOri, radMode='rsma',
                                 outPng=sumPng, zp=zp, useKpc=useKpc, pix=pix,
                                 showZoom=showZoom, exptime=exptime, bkg=bkg,
                                 outRatio=outRatio)
        else:
            """ # Run Ellipse in Forced Photometry Mode """
            print SEP
            print "##       Ellipse Run on " + \
                  "Image %s - Forced Photometry " % imgFile
            print SEP
            ellOut4 = galSBP.galSBP(imgFile, mask=mskFile,
                                    galX=galX,
                                    galY=galY,
                                    inEllip=inEllip,
                                    pix=pix,
                                    bkg=bkg,
                                    stage=4,
                                    zpPhoto=zp,
                                    maxTry=1,
                                    updateIntens=updateIntens,
                                    intMode=intMode,
                                    plMask=plMask,
                                    suffix=suffix)
            if ellOut4 is None:
                print WAR
                raise Exception("!!!!! FORCED ELLIPSE RUN FAILED !!!!")
                print WAR

            print SEP
            if psutilOk:
                mem1 = proc.memory_info().rss
                gc.collect()
                mem2 = proc.memory_info().rss
                print "@@@ Collect: %0.2f%%" % (100.0 * (mem2 - mem1) / mem0)
            else:
                gc.collect()
            print SEP

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="Prefix of the input image")
    parser.add_argument('-r', "--root", dest='root',
                        help="Directory for the input image",
                        default=None)
    parser.add_argument("--intMode", dest='intMode',
                        help="Method for integration",
                        default='mean')
    parser.add_argument('--inEllip', dest='inEllip',
                        help='Input Ellipse table',
                        default=None)
    "" "Optional """
    parser.add_argument("--suffix",
                        help="Suffix of the output image",
                        default='')
    parser.add_argument('--exMask',
                        help="External file for image mask",
                        default=None)
    parser.add_argument('--pix', dest='pix',
                        help='Pixel Scale',
                        type=float, default=0.168)
    parser.add_argument('--step', dest='step',
                        help='Step size',
                        type=float, default=0.12)
    parser.add_argument('--zp', dest='zp',
                        help='Photometric zeropoint',
                        type=float, default=27.0)
    parser.add_argument('--redshift', dest='redshift',
                        help='Photometric zeropoint',
                        type=float, default=None)
    parser.add_argument('--olthresh', dest='olthresh',
                        help='Central locator threshold',
                        type=float, default=0.50)
    parser.add_argument('--uppClip', dest='uppClip',
                        help='Upper limit for clipping',
                        type=float, default=3.0)
    parser.add_argument('--lowClip', dest='lowClip',
                        help='Upper limit for clipping',
                        type=float, default=3.0)
    parser.add_argument('--nClip', dest='nClip',
                        help='Upper limit for clipping',
                        type=int, default=2)
    parser.add_argument('--fracBad', dest='fracBad',
                        help='Outer threshold',
                        type=float, default=0.5)
    parser.add_argument('--minIt', dest='minIt',
                        help='Minimum number of iterations',
                        type=int, default=20)
    parser.add_argument('--maxIt', dest='maxIt',
                        help='Maximum number of iterations',
                        type=int, default=150)
    parser.add_argument('--maxTry', dest='maxTry',
                        help='Maximum number of attempts of ellipse run',
                        type=int, default=4)
    parser.add_argument('--galX0', dest='galX0',
                        help='Center X0',
                        type=float, default=None)
    parser.add_argument('--galY0', dest='galY0',
                        help='Center Y0',
                        type=float, default=None)
    parser.add_argument('--galQ0', dest='galQ0',
                        help='Input Axis Ratio',
                        type=float, default=None)
    parser.add_argument('--galPA0', dest='galPA0',
                        help='Input Position Angle',
                        type=float, default=None)
    parser.add_argument('--galRe', dest='galRe',
                        help='Input Effective Radius in pixel',
                        type=float, default=None)
    parser.add_argument('--outRatio', dest='outRatio',
                        help='Increase the outer boundary by this ratio',
                        type=float, default=1.2)
    parser.add_argument('--verbose', dest='verbose', action="store_true",
                        default=True)
    parser.add_argument('--psf', dest='psf', action="store_true",
                        help='Ellipse run on PSF', default=True)
    parser.add_argument('--plot', dest='plot', action="store_true",
                        help='Generate summary plot', default=True)
    parser.add_argument('--bkgCor', dest='bkgCor', action="store_true",
                        help='Background correction', default=False)
    parser.add_argument('--noCheckCenter', dest='noCheckCenter',
                        action="store_false",
                        help='Check if the center is off', default=True)
    parser.add_argument('--updateIntens', dest='updateIntens',
                        action="store_true", default=True)
    parser.add_argument('--plmask', dest='plmask', action="store_true",
                        default=True)
    parser.add_argument('--nomask', dest='nomask', action="store_true",
                        default=False)

    args = parser.parse_args()

    coaddCutoutSbp(args.prefix, root=args.root,
                   verbose=args.verbose,
                   suffix=args.suffix,
                   psf=args.psf,
                   inEllip=args.inEllip,
                   bkgCor=args.bkgCor,
                   zp=args.zp,
                   step=args.step,
                   galX0=args.galX0,
                   galY0=args.galY0,
                   galQ0=args.galQ0,
                   galPA0=args.galPA0,
                   galRe=args.galRe,
                   checkCenter=args.noCheckCenter,
                   updateIntens=args.updateIntens,
                   pix=args.pix,
                   plot=args.plot,
                   redshift=args.redshift,
                   olthresh=args.olthresh,
                   fracBad=args.fracBad,
                   lowClip=args.lowClip,
                   uppClip=args.uppClip,
                   nClip=args.nClip,
                   intMode=args.intMode,
                   minIt=args.minIt,
                   maxIt=args.maxIt,
                   maxTry=args.maxTry,
                   outRatio=args.outRatio,
                   plMask=args.plmask,
                   exMask=args.exMask,
                   noMask=args.noMask)
