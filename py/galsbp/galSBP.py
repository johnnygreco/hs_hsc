#!/usr/bin/env python
# encoding: utf-8

from __future__ import division

import os
import copy
import argparse
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
from astropy.io import ascii
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

# Personal
import hscUtils as hUtil

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
    iraf.unlearn('imcopy')
    iraf.imcopy.input  = "'" + inputFits + "'"
    iraf.imcopy.output = "'" + mskIraf + "'"
    iraf.imcopy()

    return mskIraf


def defaultEllipse(x0, y0, maxsma, ellip0=0.05, pa0=0.0, sma0=6.0, minsma=0.0,
        linear=False, step=0.08, recenter=True, conver=0.05, hcenter=True,
        hellip=True, hpa=True, minit=10, maxit=250, olthresh=1.00, mag0=27.0,
        integrmode='median', usclip=3.0, lsclip=3.0, nclip=2, fflag=0.5,
        harmonics='none'):
    """
    doc
    """
    ellipConfig = np.recarray((1,), dtype=[('x0', float), ('y0', float),
                                   ('ellip0', float), ('pa0', float),
                                   ('sma0', float), ('minsma', float), ('maxsma', float),
                                   ('linear', bool), ('step', float),
                                   ('recenter', bool), ('conver', int),
                                   ('hcenter', bool), ('hellip', bool), ('hpa', bool),
                                   ('minit', int), ('maxit', int), ('olthresh', float),
                                   ('mag0', float), ('integrmode', 'a10'),
                                   ('usclip', float), ('lsclip', float), ('nclip', int),
                                   ('fflag', float), ('harmonics', 'a10')])
    # Default setting for Ellipse Run
    ellipConfig['x0'] = x0
    ellipConfig['y0'] = y0
    ellipConfig['ellip0'] = ellip0
    ellipConfig['pa0'] = pa0
    ellipConfig['sma0'] = sma0
    ellipConfig['minsma'] = minsma
    ellipConfig['maxsma'] = maxsma
    ellipConfig['linear'] = linear
    ellipConfig['step'] = step
    ellipConfig['recenter'] = recenter
    ellipConfig['conver'] = conver
    ellipConfig['hcenter'] = hcenter
    ellipConfig['hellip'] = hellip
    ellipConfig['hpa'] = hpa
    ellipConfig['minit'] = minit
    ellipConfig['maxit'] = maxit
    ellipConfig['olthresh'] = olthresh
    ellipConfig['mag0'] = mag0
    ellipConfig['integrmode'] = integrmode
    ellipConfig['usclip'] = usclip
    ellipConfig['lsclip'] = lsclip
    ellipConfig['nclip'] = nclip
    ellipConfig['fflag'] = fflag
    ellipConfig['harmonics'] = harmonics

    return ellipConfig


def easierEllipse(ellipConfig):
    """
    doc
    """
    ellipConfig['maxsma'] *= 0.9
    ellipConfig['step']  += 0.04
    ellipConfig['fflag'] += 0.05

    return ellipConfig


def setupEllipse(ellipConfig):
    """TODO: Docstring for ellipseConfig
    :returns: TODO

    """

    cfg = ellipConfig[0]

    iraf.unlearn(iraf.ellipse)
    iraf.unlearn(iraf.geompar)
    iraf.unlearn(iraf.controlpar)
    iraf.unlearn(iraf.samplepar)
    iraf.unlearn(iraf.magpar)

    # Define parameters for the ellipse run
    # 1. Initial guess of the central X, Y
    if (cfg['x0'] > 0) and (cfg['y0'] > 0):
        iraf.ellipse.x0 = cfg['x0']
        iraf.ellipse.y0 = cfg['y0']
    else:
        raise "Make sure that the input X0 and Y0 are meaningful !", cfg['x0'], cfg['y0']

    # 2. Initial guess of the ellipticity and PA of the first ISOPHOTE
    if (cfg['ellip0'] >= 0.0) and (cfg['ellip0'] < 1.0):
        iraf.ellipse.ellip0 = cfg['ellip0']
    else:
        raise "Make sure that the input Ellipticity is meaningful !", cfg['ellip0']
    if (cfg['pa0'] >= 0.0) and (cfg['pa0'] <= 180.0):
        iraf.ellipse.pa0 = cfg['pa0']
    else:
        raise "Make sure that the input Position Angle is meaningful !", cfg['pa0']

    # 3. Initial radius for ellipse fitting
    iraf.ellipse.sma0       = cfg['sma0']
    # 4. The minimum and maximum radius for the ellipse fitting
    iraf.ellipse.minsma     = cfg['minsma']
    iraf.ellipse.maxsma     = cfg['maxsma']
    # 5. Parameters about the stepsize during the fitting.
    if cfg['linear']:
        iraf.ellipse.linear     = 'yes'
    else:
        iraf.ellipse.linear     = 'no'
    iraf.ellipse.geompar.step       = cfg['step']
    # 6. Do you want to allow the ellipse to decide the galaxy center during the
    if cfg['recenter']:
        iraf.ellipse.recenter   = 'yes'
    else:
        iraf.ellipse.recenter   = 'no'
    # 7. The next three parameters control the behavior of the fit
    iraf.ellipse.conver  = cfg['conver']
    if cfg['hcenter']:
        iraf.ellipse.hcenter = 'yes'
    else:
        iraf.ellipse.hcenter = 'no'
    if cfg['hellip']:
        iraf.ellipse.hellip  = 'yes'
    else:
        iraf.ellipse.hellip  = "no"
    if cfg['hpa']:
        iraf.ellipse.hpa     = 'yes'
    else:
        iraf.ellipse.hpa     = 'no'
    # 8. Parameters about the iterations
    # minit/maxit: minimun and maximum number of the iterations
    iraf.ellipse.minit   = cfg['minit']
    iraf.ellipse.maxit   = cfg['maxit']
    # 9. Threshold for the object locator algorithm
    iraf.ellipse.olthresh = cfg['olthresh']
    # 10. Make sure the Interactive Mode is turned off
    iraf.ellipse.interactive         = 'no'
    # 11. Magnitude Zeropoint
    iraf.ellipse.mag0         = cfg['mag0']
    # 12. Sampler
    intMode = cfg['integrmode']
    intMode = intMode.lower().strip()
    if intMode == 'median':
        iraf.ellipse.integrmode  = 'median'
    elif intMode == 'mean':
        iraf.ellipse.integrmode  = 'mean'
    elif intMode == 'bi-linear':
        iraf.ellipse.integrmode  = 'bi-linear'
    else:
        raise Exception("### Only 'mean', 'median', and 'bi-linear' are available !")
    iraf.ellipse.usclip      = cfg['usclip']
    iraf.ellipse.lsclip      = cfg['lsclip']
    iraf.ellipse.nclip       = cfg['nclip']
    iraf.ellipse.fflag       = cfg['fflag']
    # 13. Optional Harmonics
    iraf.ellipse.harmonics   = cfg['harmonics']


def ellipRemoveIndef(outTabName, replace='NaN'):

    """TODO: Docstring for ellipRemoveIndef
    :returns: TODO

    """

    if os.path.exists(outTabName):
        subprocess.call(['sed', '-i_back', 's/INDEF/' + replace + '/g', outTabName])
    else:
        raise Exception('Can not find the input catalog!')

    return outTabName


def readEllipseOut(outTabName, pix=1.0, zp=27.0, exptime=1.0, bkg=0.0,
                   harmonics='none'):

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
    if harmonics != "none":
        # TODO: Read as many harmonics as necessary
        ellipseOut.rename_column('col41', 'a1')
        ellipseOut.rename_column('col42', 'a1_err')
        ellipseOut.rename_column('col43', 'b1')
        ellipseOut.rename_column('col44', 'b1_err')
        ellipseOut.rename_column('col45', 'a2')
        ellipseOut.rename_column('col46', 'a2_err')
        ellipseOut.rename_column('col47', 'b2')
        ellipseOut.rename_column('col48', 'b2_err')

    # Normalize the PA
    ellipseOut.add_column(Column(name='pa_norm',
                          data=np.array([hUtil.normAngle(pa, lower=0.0, upper=180.0)
                              for pa in ellipseOut['pa']])))

    # Convert the intensity into surface brightness
    parea = (pix ** 2.0)
    ellipseOut.add_column(Column(name='sbp',
                          data=(zp-2.5*np.log10((ellipseOut['intens']-bkg)/
                              (parea*exptime)))))
    ellipseOut.add_column(Column(name='sbp_lerr',
                          data=(zp-2.5*np.log10((ellipseOut['intens']-
                              ellipseOut['int_err']-bkg)/(parea*exptime)))))
    ellipseOut.add_column(Column(name='sbp_uerr',
                          data=(zp-2.5*np.log10((ellipseOut['intens']+
                              ellipseOut['int_err']-bkg)/(parea*exptime)))))

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
        raise "The x and y should have the same size!", (len(ellipseRsma),
                len(isoGrowthCurve))
    else:
        if simple:
            isoRsma50 = ellipseRsma[np.nanargmin(np.abs(isoGrowthCurve - 50.0))]
        else:
            isoRsma50 = (numpy.interp([50.0], isoGrowthCurve, ellipseRsma))[0]

    return isoRsma50


def plotIsoMap(img, ellipOut, mask=None, maxRad=None):

    """
    doc
    """


def plotSbpProf(ellipOut, psf=None, maxRad=None):

    """
    doc
    """

def saveEllipOut(ellipOut, prefix):
    """
    doc
    """


def galSBP(image, mask, galX=None, galY=None, inEllip=None, maxSma=None, iniSma=6.0,
           galR=10.0, galQ=0.9, galPA=0.0, pix=0.168, bkg=0.00, stage=3, minSma=0.0,
           gain=3.0, expTime=1.0, zpPhoto=27.0, maxTry=2, minIt=10, maxIt=200,
           ellipStep=0.08, uppClip=3.0, lowClip=3.0, nClip=3, fracBad=0.5,
           intMode="median", suffix=None, plMask=True, conver=0.05, recenter=True,
           verbose=True, visual=True, linearStep=False,
           olthresh=1.00, harmonics='1 2'):

    """TODO: Docstring for galEllipse.

    stage  = 1: All Free
            2: Center Fixed
            3: All geometry fixd
            4: Force Photometry, must have inEllip
    :returns: TODO

    """

    if verbose:
        verStr = 'yes'
    else:
        verStr = 'no'

    """ Check input files """
    if not os.path.isfile(image):
        raise Exception("### Can not find the input image: %s !" % image)
    if not os.path.isfile(mask):
        raise Exception("### Can not find the input mask: %s !" % mask)

    """ Conver the .fits mask to .pl file if necessary """
    if not plMask:
        plFile = maskFits2Pl(mask)
        if not os.path.isfile(plFile):
            raise Exception("### Can not find the .pl mask: %s !" % plFile)
        else:
            mask = plFile

    """ Estimate the maxSMA if none is provided """
    if (maxSma is None) or (galX is None) or (galY is None):
        # TODO: Option to use any HDU
        data = (fits.open(image))[0].data
        dimX, dimY = data.shape
        imgSize = dimX if (dimX >= dimY) else dimY
        imgR = imgSize / 2.0 * 1.3
        imgX = dimX / 2.0
        imgY = dimY / 2.0
        if maxSma is None:
            maxSma = imgR
        if galX is None:
            galX = imgX
        if galY is None:
            galY = imgY

    """ Check the stage """
    if stage == 1:
        hcenter, hellip, hpa = False, False, False
    elif stage == 2:
        hcenter, hellip, hpa = True, False, False
    elif stage == 3:
        print "Fix Geomtry "
        hcenter, hellip, hpa = True, True, True
    elif stage == 4:
        hcenter, hellip, hpa = True, True, True
        if (inEllip is None ) or (not os.path.isfile(inEllip)):
            raise Exception("### Can not find the input ellip file: %s !" % inEllip)
    else:
        raise Exception("### Available step: 1 , 2 , 3 , 4")

    """ Get the default Ellipse settings """
    galEll = (1.0 - galQ)
    ellipCfg = defaultEllipse(galX, galY, maxSma, ellip0=galEll, pa0=galPA,
                sma0=iniSma, minsma=minSma, linear=linearStep, step=ellipStep,
                recenter=recenter, conver=conver, hcenter=hcenter,
                hellip=hellip, hpa=hpa, minit=minIt, maxit=maxIt, olthresh=olthresh,
                mag0=zpPhoto, integrmode=intMode, usclip=uppClip, lsclip=lowClip,
                nclip=nClip, fflag=fracBad, harmonics=harmonics)

    """ Name of the output files """
    if suffix is None:
        suffix = 'ellip'
    suffix = '_' + suffix + '_' + str(stage).strip()
    # TODO: file extension can also be fit or other
    outBin = image.replace('.fits', suffix + '.bin')
    outTab = image.replace('.fits', suffix + '.tab')
    outCdf = image.replace('.fits', suffix + '.cdf')

    outPkl = image.replace('.fits', suffix + '.pkl')
    outCsv = image.replace('.fits', suffix + '.csv')
    outPng = image.replace('.fits', suffix + '.png')
    outCfg = image.replace('.fits', suffix + '.cfg')

    """ Call the STSDAS.ANALYSIS.ISOPHOTE package """
    iraf.stsdas()
    iraf.analysis()
    iraf.isophote()

    """ Start the Ellipse Run """
    attempts = 0
    while attempts < maxTry:
        try:
            """ Config the parameters for ellipse """
            print ellipCfg
            setupEllipse(ellipCfg)
            """ Ellipse run """
            # Check and remove outputs from the previous Ellipse run
            if os.path.exists(outBin):
                os.remove(outBin)
            # Start the Ellipse fitting
            if stage != 4:
                iraf.ellipse(input=image, output=outBin, verbose=verStr)
            else:
                iraf.ellipse(input=image, output=outBin, inellip=inEllip,
                        verbose=verStr)
            break
        except Exception, err:
            print err
            attempts +=1
            ellipCfg = easierEllipse(ellipCfg)
            print "###  !!! Make the Ellipse Run A Little Bit Easier !"

    # Check if the Ellipse run is finished
    if not os.path.isfile(outBin):
        ellipOut = None
        print "###    XXX ELLIPSE RUN FAILED AFTER %3d ATTEMPTS!!!" % maxTry
    else:
        # Remove the existed .tab and .cdf file
        if os.path.isfile(outTab):
            os.remove(outTab)
        if os.path.isfile(outCdf):
            os.remove(outCdf)
        # Tdump the .bin table into a .tab file
        iraf.unlearn('tdump')
        iraf.tdump.columns=''
        iraf.tdump(outBin, datafil=outTab, cdfile=outCdf)
        # Read in the Ellipse output tab
        ellipOut = readEllipseOut(outTab, zp=zpPhoto, pix=pix, exptime=expTime,
                bkg=bkg, harmonics=harmonics)
        nIso = len(ellipOut)
        maxRad = np.nanmax(ellipOut['sma'])
        if verbose:
            print "###   %d elliptical isophotes have been extracted" % nIso
            print "###   The maximum radius is %7.2f pixels" % maxRad

        """ Save a Pickle file """
        hUtil.saveToPickle(ellipOut, outPkl)
        if os.path.isfile(outPkl):
            if verbose:
                print "###     Save Ellipse output to .pkl file: %s" % outPkl
        else:
            raise Exception("### Something is wrong with the .pkl file")

        """ Save the current configuration to a .pkl file """
        hUtil.saveToPickle(ellipCfg, outCfg)
        if os.path.isfile(outCfg):
            if verbose:
                print "###     Save Ellipse configuration to a .cfg file: %s" % outCfg
        else:
            raise Exception("### Something is wrong with the .pkl file")

        """ Save a .CSV file """
        ascii.write(ellipOut, outCsv, format='csv')
        if os.path.isfile(outCsv):
            if verbose:
                print "###     Save Ellipse output to .csv file: %s" % outCsv
        else:
            raise Exception("### Something is wrong with the .csv file")

    return ellipOut


#def galSBP(image, mask, galX=None, galY=None, inEllip=None, maxSma=None, iniSma=6.0,
           #galR=10.0, galQ=0.9, galPA=0.0, pix=0.168, bkg=0.00, step=0.08,
           #gain=3.0, expTime=1.0, zpPhoto=27.0, maxTry=2, minIt=10, maxIt=200,
           #ellStep=0.08, uppClip1=3.0, lowClip=3.0, nClip=3, fracBad=0.5,
           #intMode="median", suffix=None, plMask=True, conver=2, recenter=True,
           #verbose=True, visual=True, harmonic=False, linearStep=False,
           #olthresh=1.00, harmonics='1,2'):
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("image", help="Name of the input image")
    parser.add_argument("mask", help="Name of the input mask")
    parser.add_argument('-x0', dest='galX', help='Galaxy center in X-dimension',
                       type=float, default=None)
    parser.add_argument('-y0', dest='galY', help='Galaxy center in Y-dimension',
                       type=float, default=None)
    parser.add_argument('--maxsma', dest='maxSma', help='Maximum radius for Ellipse Run',
                       type=float, default=None)
    parser.add_argument('-s', dest='stage', help='Stage of Ellipse Run',
                       type=int, default=3, choices=range(1, 4))

    args = parser.parse_args()

    galSBP(args.image, args.mask, galX=args.galX, galY=args.galY,
           maxSma=args.maxSma, stage=args.stage)
