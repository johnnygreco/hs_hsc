#!/usr/bin/env python
# encoding: utf-8
"""Single Sersic fits to HSC cutouts."""

from __future__ import division

import os
import copy
import argparse
import subprocess
import numpy as np
import scipy
from distutils import spawn

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


def coaddRunGalfit(readFile, root=None, imax=120, galfit=None, updateRead=True):
    """
    Run GALFIT
    """

    """ Find GALFIT """
    if galfit is None:
        galfit = spawn.find_executable('galfit')
        if galfit is None:
            raise Exception("XXX Can not find the GALFIT executable")

    """ Check the Read-in File """
    if not os.path.isfile(readFile):
        raise Exception("XXX Can not find the READIN file: %s", readFile)

    """ IMAX string """
    imaxStr = " -imax %4d" % imax

    """ GALFIT command """
    galfitCommand = galfit + ' ' + imaxStr + ' ' + readFile

    """ Excecute the command """
    if root is not None:
        proc = subprocess.Popen([galfitCommand], cwd=root, shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    else:
        proc = subprocess.Popen([galfitCommand], shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in proc.stdout.readlines():
        print line
    retval = proc.wait()

    return


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


def readInputSky(prefix, root=None, rebin='rebin6'):
    """
    doc
    """

    skyFile = prefix + '_' + rebin + '_sky.dat'
    if not os.path.isfile(skyFile):
        raise Exception("### Can not find the input sky summary : %s !" % skyFile)

    skySum = open(skyFile, 'r').readlines()
    skyMed = float(skySum[3].split(':')[1].strip())
    skyAvg = float(skySum[4].split(':')[1].strip())
    skyStd = float(skySum[5].split(':')[1].strip())
    skySkw = float(skySum[6].split(':')[1].strip())

    return skyMed, skyAvg, skyStd


def getInput1Sersic(config, readinFile='cutout_1ser.in', skyGrad=True, useF1=False,
        useF4=False):
    """
    Generate the readin file for 1 Sersic GALFIT fitting
    """

    f = open(readinFile, 'w')

    f.write('\n')
    f.write('===============================================================================\n')
    f.write('# IMAGE and GALFIT CONTROL PARAMETERS\n')
    f.write('A) %s  # Input data image (FITS file)\n' % config['image'][0])
    f.write('B) %s  # Output data image block\n' % config['output'][0])
    f.write('C) %s  # Sigma image name\n' % config['sig'][0])
    f.write('D) %s  # Input PSF image \n' % config['psf'][0])
    f.write('E) 1                   # PSF fine sampling factor relative to data \n')
    f.write('F) %s  # Bad pixel mask\n' % config['mask'][0])
    f.write('G) %s  # File with parameter constraints \n' % config['constr'][0])
    f.write('H)  1  %5d  1  %5d  # Image region to fit\n' % (config['dimx'],
        config['dimy']))
    f.write('I) %5d %5d  # Size of the convolution box\n' % (config['convbox'],
        config['convbox']))
    f.write('J) %6.2f  # Magnitude photometric zeropoint \n' % config['zp'])
    f.write('K) %7.3f  %7.3f # Plate scale (dx dy)\n' % (config['pix'], config['pix']))
    f.write('O) regular             # Display type (regular, curses, both)\n')
    f.write('P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n')
    f.write('\n')

    f.write('# INITIAL FITTING PARAMETERS\n')
    f.write('#\n')
    f.write('#   For object type, the allowed functions are: \n')
    f.write('#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat, \n')
    f.write('#       ferrer, powsersic, sky, and isophote. \n')
    f.write('#  \n')
    f.write('#   Hidden parameters will only appear when they''re specified:\n')
    f.write('#       C0 (diskyness/boxyness), \n')
    f.write('#       Fn (n=integer, Azimuthal Fourier Modes),\n')
    f.write('#       R0-R10 (PA rotation, for creating spiral structures).\n')
    f.write('# \n')
    f.write('# -----------------------------------------------------------------------------\n')
    f.write('#   par)    par value(s)    fit toggle(s)    # parameter description \n')
    f.write('# -----------------------------------------------------------------------------\n')
    f.write('\n')
    f.write('# Object number: 1\n')
    f.write(' 0) sersic    \n')
    f.write(' 1) %7.1f %7.1f  1 1 \n' % (config['x'], config['y']))
    f.write(' 3) %7.3f     1 \n' % config['mag'])
    f.write(' 4) %7.3f     1 \n' % config['re'])
    f.write(' 5) %7.3f     1 \n' % config['nser'])
    f.write(' 6) 0.0000      0          #     ----- \n')
    f.write(' 7) 0.0000      0          #     ----- \n')
    f.write(' 8) 0.0000      0          #     ----- \n')
    f.write(' 9) %7.3f     1 \n' % config['ba'])
    f.write('10) %7.3f     1 \n' % config['pa'])
    if useF1:
        f.write('F1) 0.01 10.00 1 1 ')
    if useF4:
        f.write('F4) 0.01 10.00 1 1 ')
    f.write(' Z) 0                      #  output option (0 = resid., 1 = Dont subtract) \n')
    f.write('\n')
    if config['usesky'] == 1:
        f.write('# Object number: 2\n')
        f.write(' 0) sky                    #  object type\n')
        f.write(' 1) %8.3f  1  #  sky background \n' % config['bkg'])
        if skyGrad:
            f.write(' 2) 0.0000      1          #  dsky/dx (sky gradient in x)\n')
            f.write(' 3) 0.0000      1          #  dsky/dy (sky gradient in y)\n')
        else:
            f.write(' 2) 0.0000      0          #  dsky/dx (sky gradient in x)\n')
            f.write(' 3) 0.0000      0          #  dsky/dy (sky gradient in y)\n')
        f.write(' Z) 0                      #  output option (0 = resid., 1 = Dont subtract) \n')
        f.write('\n')
    f.write('================================================================================\n')

    f.close()


def getInput2Sersic(config, readinFile='cutout_2ser.in', constr=False, skyGrad=True):
    """
    Generate the readin file for 2 Sersic GALFIT fitting
    """

    f = open(readinFile, 'w')

    f.write('\n')
    f.write('===============================================================================\n')
    f.write('# IMAGE and GALFIT CONTROL PARAMETERS\n')
    f.write('A) %s  # Input data image (FITS file)\n' % config['image'][0])
    f.write('B) %s  # Output data image block\n' % config['output'][0])
    f.write('C) %s  # Sigma image name\n' % config['sig'][0])
    f.write('D) %s  # Input PSF image \n' % config['psf'][0])
    f.write('E) 1                   # PSF fine sampling factor relative to data \n')
    f.write('F) %s  # Bad pixel mask\n' % config['mask'][0])
    f.write('G) %s  # File with parameter constraints \n' % config['constr'][0])
    f.write('H)  1  %5d  1  %5d  # Image region to fit\n' % (config['dimx'],
        config['dimy']))
    f.write('I) %5d %5d  # Size of the convolution box\n' % (config['convbox'],
        config['convbox']))
    f.write('J) %6.2f  # Magnitude photometric zeropoint \n' % config['zp'])
    f.write('K) %7.3f  %7.3f # Plate scale (dx dy)\n' % (config['pix'], config['pix']))
    f.write('O) regular             # Display type (regular, curses, both)\n')
    f.write('P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n')
    f.write('\n')
    f.write('# INITIAL FITTING PARAMETERS\n')
    f.write('#\n')
    f.write('#   For object type, the allowed functions are: \n')
    f.write('#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat, \n')
    f.write('#       ferrer, powsersic, sky, and isophote. \n')
    f.write('#  \n')
    f.write('#   Hidden parameters will only appear when they''re specified:\n')
    f.write('#       C0 (diskyness/boxyness), \n')
    f.write('#       Fn (n=integer, Azimuthal Fourier Modes),\n')
    f.write('#       R0-R10 (PA rotation, for creating spiral structures).\n')
    f.write('# \n')
    f.write('# -----------------------------------------------------------------------------\n')
    f.write('#   par)    par value(s)    fit toggle(s)    # parameter description \n')
    f.write('# -----------------------------------------------------------------------------\n')
    f.write('\n')

    f.write('# Object number: 1\n')
    f.write(' 0) sersic    \n')
    f.write(' 1) %7.1f %7.1f  1 1 \n' % (config['x'], config['y']))
    f.write(' 3) %7.3f     1 \n' % (config['mag']+0.6))
    f.write(' 4) %7.3f     1 \n' % (config['re']*0.25))
    f.write(' 5) %7.3f     1 \n' % config['nser'])
    f.write(' 6) 0.0000      0          #     ----- \n')
    f.write(' 7) 0.0000      0          #     ----- \n')
    f.write(' 8) 0.0000      0          #     ----- \n')
    f.write(' 9) %7.3f     1 \n' % config['ba'])
    f.write('10) %7.3f     1 \n' % config['pa'])
    if useF1:
        f.write('F1) 0.01 10.00 1 1 ')
    if useF4:
        f.write('F4) 0.01 10.00 1 1 ')
    f.write(' Z) 0                      #  output option (0 = resid., 1 = Dont subtract) \n')

    f.write('# Object number: 2\n')
    f.write(' 0) sersic    \n')
    f.write(' 1) %7.1f %7.1f  1 1 \n' % (config['x'], config['y']))
    f.write(' 3) %7.3f     1 \n' % (config['mag']))
    f.write(' 4) %7.3f     1 \n' % (config['re']*1.5))
    f.write(' 5) 0.9    1 \n')
    f.write(' 6) 0.0000      0          #     ----- \n')
    f.write(' 7) 0.0000      0          #     ----- \n')
    f.write(' 8) 0.0000      0          #     ----- \n')
    f.write(' 9) %7.3f     1 \n' % config['ba'])
    f.write('10) %7.3f     1 \n' % config['pa'])
    if useF1:
        f.write('F1) 0.01 10.00 1 1 ')
    if useF4:
        f.write('F4) 0.01 10.00 1 1 ')
    f.write(' Z) 0                      #  output option (0 = resid., 1 = Dont subtract) \n')
    f.write('\n')

    if config['usesky'] == 1:
        f.write('# Object number: 3\n')
        f.write(' 0) sky                    #  object type\n')
        f.write(' 1) %8.3f  1  #  sky background \n' % config['bkg'])
        if skyGrad:
            f.write(' 2) 0.0000      1          #  dsky/dx (sky gradient in x)\n')
            f.write(' 3) 0.0000      1          #  dsky/dy (sky gradient in y)\n')
        else:
            f.write(' 2) 0.0000      0          #  dsky/dx (sky gradient in x)\n')
            f.write(' 3) 0.0000      0          #  dsky/dy (sky gradient in y)\n')
        f.write(' Z) 0                      #  output option (0 = resid., 1 = Dont subtract) \n')
        f.write('\n')
    f.write('================================================================================\n')

    f.close()


def getInput3Sersic(config, readinFile='cutout_3ser.in', constr=False, skyGrad=True):
    """
    Generate the readin file for 3 Sersic GALFIT fitting
    """

    f = open(readinFile, 'w')

    f.write('\n')
    f.write('===============================================================================\n')
    f.write('# IMAGE and GALFIT CONTROL PARAMETERS\n')
    f.write('A) %s  # Input data image (FITS file)\n' % config['image'][0])
    f.write('B) %s  # Output data image block\n' % config['output'][0])
    f.write('C) %s  # Sigma image name\n' % config['sig'][0])
    f.write('D) %s  # Input PSF image \n' % config['psf'][0])
    f.write('E) 1                   # PSF fine sampling factor relative to data \n')
    f.write('F) %s  # Bad pixel mask\n' % config['mask'][0])
    f.write('G) %s  # File with parameter constraints \n' % config['constr'][0])
    f.write('H)  1  %5d  1  %5d  # Image region to fit\n' % (config['dimx'],
        config['dimy']))
    f.write('I) %5d %5d  # Size of the convolution box\n' % (config['convbox'],
        config['convbox']))
    f.write('J) %6.2f  # Magnitude photometric zeropoint \n' % config['zp'])
    f.write('K) %7.3f  %7.3f # Plate scale (dx dy)\n' % (config['pix'], config['pix']))
    f.write('O) regular             # Display type (regular, curses, both)\n')
    f.write('P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n')
    f.write('\n')
    f.write('# INITIAL FITTING PARAMETERS\n')
    f.write('#\n')
    f.write('#   For object type, the allowed functions are: \n')
    f.write('#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat, \n')
    f.write('#       ferrer, powsersic, sky, and isophote. \n')
    f.write('#  \n')
    f.write('#   Hidden parameters will only appear when they''re specified:\n')
    f.write('#       C0 (diskyness/boxyness), \n')
    f.write('#       Fn (n=integer, Azimuthal Fourier Modes),\n')
    f.write('#       R0-R10 (PA rotation, for creating spiral structures).\n')
    f.write('# \n')
    f.write('# -----------------------------------------------------------------------------\n')
    f.write('#   par)    par value(s)    fit toggle(s)    # parameter description \n')
    f.write('# -----------------------------------------------------------------------------\n')

    f.write('\n')
    f.write('# Object number: 1\n')
    f.write(' 0) sersic    \n')
    f.write(' 1) %7.1f %7.1f  1 1 \n' % (config['x'], config['y']))
    f.write(' 3) %7.3f     1 \n' % (config['mag']+1.2))
    f.write(' 4) %7.3f     1 \n' % (config['re']*0.25))
    f.write(' 5) %7.3f     1 \n' % config['nser'])
    f.write(' 6) 0.0000      0          #     ----- \n')
    f.write(' 7) 0.0000      0          #     ----- \n')
    f.write(' 8) 0.0000      0          #     ----- \n')
    f.write(' 9) %7.3f     1 \n' % config['ba'])
    f.write('10) %7.3f     1 \n' % config['pa'])
    if useF1:
        f.write('F1) 0.01 10.00 1 1 ')
    if useF4:
        f.write('F4) 0.01 10.00 1 1 ')
    f.write(' Z) 0                      #  output option (0 = resid., 1 = Dont subtract) \n')

    f.write('# Object number: 2\n')
    f.write(' 0) sersic    \n')
    f.write(' 1) %7.1f %7.1f  1 1 \n' % (config['x'], config['y']))
    f.write(' 3) %7.3f     1 \n' % (config['mag']+0.9))
    f.write(' 4) %7.3f     1 \n' % (config['re']*0.9))
    f.write(' 5) 0.9    1 \n')
    f.write(' 6) 0.0000      0          #     ----- \n')
    f.write(' 7) 0.0000      0          #     ----- \n')
    f.write(' 8) 0.0000      0          #     ----- \n')
    f.write(' 9) %7.3f     1 \n' % config['ba'])
    f.write('10) %7.3f     1 \n' % config['pa'])
    if useF1:
        f.write('F1) 0.01 10.00 1 1 ')
    if useF4:
        f.write('F4) 0.01 10.00 1 1 ')
    f.write(' Z) 0                      #  output option (0 = resid., 1 = Dont subtract) \n')
    f.write('\n')

    f.write('# Object number: 3\n')
    f.write(' 0) sersic    \n')
    f.write(' 1) %7.1f %7.1f  1 1 \n' % (config['x'], config['y']))
    f.write(' 3) %7.3f     1 \n' % (config['mag']+0.7))
    f.write(' 4) %7.3f     1 \n' % (config['re']*1.3))
    f.write(' 5) 0.5    1 \n')
    f.write(' 6) 0.0000      0          #     ----- \n')
    f.write(' 7) 0.0000      0          #     ----- \n')
    f.write(' 8) 0.0000      0          #     ----- \n')
    f.write(' 9) %7.3f     1 \n' % config['ba'])
    f.write('10) %7.3f     1 \n' % config['pa'])
    if useF1:
        f.write('F1) 0.01 10.00 1 1 ')
    if useF4:
        f.write('F4) 0.01 10.00 1 1 ')
    f.write(' Z) 0                      #  output option (0 = resid., 1 = Dont subtract) \n')
    f.write('\n')

    if config['usesky'] == 1:
        f.write('# Object number: 3\n')
        f.write(' 0) sky                    #  object type\n')
        f.write(' 1) %8.3f  1  #  sky background \n' % config['bkg'])
        if skyGrad:
            f.write(' 2) 0.0000      1          #  dsky/dx (sky gradient in x)\n')
            f.write(' 3) 0.0000      1          #  dsky/dy (sky gradient in y)\n')
        else:
            f.write(' 2) 0.0000      0          #  dsky/dx (sky gradient in x)\n')
            f.write(' 3) 0.0000      0          #  dsky/dy (sky gradient in y)\n')
        f.write(' Z) 0                      #  output option (0 = resid., 1 = Dont subtract) \n')
        f.write('\n')
    f.write('================================================================================\n')

    f.close()



def coaddCutoutGalfitSimple(prefix, root=None, pix=0.168, useBkg=True, zp=27.0, usePsf=True,
        galX0=None, galY0=None, galQ0=None, galPA0=None, galRe=None, galSer=2.0,
        model=None, inFile=None, outFile=None, useSig=True, mag=18.0, constrFile=None,
        verbose=True, run=True, skyGrad=True, ser2Comp=True, ser3comp=True,
        useF4=False, useF1=False):

    """
    Run 1-Sersic fitting on HSC cutout image
    """

    print "## Input Image: ", prefix
    """ 0. Organize Input Data """
    # Read in the input image, mask, psf, and their headers
    imgFile, imgArr, imgHead, mskFile, mskArr, mskHead = readSbpInput(prefix,
            root=root)
    if not imgSameSize(imgArr, mskArr):
        raise Exception("### The Image and Mask need to have EXACTLY same dimensions!")
    dimX, dimY = imgArr.shape

    if mskHead['MSK_R20'] == 1:
        raise Exception("### The central region is masked out")

    """ 0a. PSF """
    if usePsf:
        psfFile = prefix + '_psf.fits'
        if root is not None:
            psfFile = os.path.join(root, psfFile)
        if not os.path.isfile(psfFile):
            raise Exception(" XXX Can not find the PSF image : %s", psfFile)
    else:
        psfFile = ''

    """ 0b. Sigma Image """
    if useSig:
        sigFile = prefix + '_sig.fits'
        if root is not None:
            sigFile = os.path.join(root, sigFile)
        if not os.path.isfile(sigFile):
            raise Exception(" XXX Can not find the Sigma image : %s", sigFile)
    else:
        sigFile = ''

    """ 0c. Background """
    if useBkg:
        try:
            skyMed, skyAvg, skyStd = readInputSky(prefix, root=root)
            bkg = skyAvg
            if verbose:
                print " ### Average Background : ", bkg
        except Exception:
            print " XXX CAN NOT FIND THE BACKGROUND DATA !"
            bkg = 0.00
    else:
        bkg = 0.00

    """ 0d. Read-in File """
    if model is None:
        suffix = ''
    else:
        suffix = '_' + suffix
    if inFile is None:
        inFile = prefix + '_1ser' + suffix + '.in'
        if root is not None:
            inFile = os.path.join(root, inFile)

    """ 0e. Output File """
    if outFile is None:
        outFile = prefix + '_1ser' + suffix + '.fits'
        if root is not None:
            outFile = os.path.join(root, outFile)

    """ 0f. Prepare the Input for SBP """
    if (galX0 is None) or (galY0 is None):
        galX, galY = mskHead['GAL_CENX'], mskHead['GAL_CENY']
    else:
        galX, galY = galX0, galY0
    if (galQ0 is None) or (galPA0 is None):
        galQ, galPA = mskHead['GAL_Q'], mskHead['GAL_PA']
    else:
        galQ, galPA = galQ0, galPA0
    galQ = galQ if galQ <= 0.95 else 0.95
    galPA = hUtil.normAngle(galPA, lower=0.0, upper=180.0)

    if galRe is None:
        galR50 = mskHead['GAL_R50']
    else:
        galR50 = galRe

    """ 0g. Convolution Box Size """
    convbox = int(galR50 * 26.0)
    convbox = convbox if convbox <= int(dimX*0.9) else int(dimX*0.9)

    if verbose:
        print " ### Image : ", imgFile
        print " ### Mask  : ", mskFile
        print " ### Sigma : ", sigFile
        print " ### PSF   : ", psfFile
        print " ### galX, galY : ", galX, galY
        print " ### galQ, galPA : ", galQ, galPA
        print " ### galR50 : ", galR50
        print " ### galSer : ", galSer
        print " ### convbox : ", convbox

    """ 0h. Generate the configuration file """
    galfitConfig = np.recarray((1,), dtype=[('x', float), ('y', float),
                                   ('ba', float), ('pa', float), ('mag', float),
                                   ('re', float), ('nser', float), ('bkg', float),
                                   ('image', 'a120'), ('psf', 'a120'), ('mask', 'a120'),
                                   ('constr', 'a50'), ('sig', 'a120'), ('pix', float),
                                   ('zp', float), ('convbox', int), ('usesky', int),
                                   ('dimx', int), ('dimy', int), ('output', 'a120')
                                   ])
    if useBkg:
        galfitConfig['usesky'] = 1
    else:
        galfitConfig['usesky'] = 0

    galfitConfig['x']     = galX
    galfitConfig['y']     = galX
    galfitConfig['ba']    = galQ
    galfitConfig['pa']    = galPA
    galfitConfig['mag']   = mag
    galfitConfig['re']    = galR50
    galfitConfig['nser']  = galSer
    galfitConfig['bkg']   = bkg
    galfitConfig['image'] = imgFile
    galfitConfig['psf']   = psfFile
    galfitConfig['sig']   = sigFile
    galfitConfig['mask']  = mskFile
    galfitConfig['pix']   = pix
    galfitConfig['zp']    = zp
    galfitConfig['convbox'] = convbox
    galfitConfig['dimx']  = dimX
    galfitConfig['dimy']  = dimY
    galfitConfig['output']  = outFile
    galfitConfig['constr']  = ''

    """ 1a. Generate the Read-in File for 1Ser model"""
    getInput1Sersic(galfitConfig, readinFile=inFile, skyGrad=skyGrad,
            useF1=useF1, useF4=useF4)
    """ 1b. Execute the GALFIT run """
    if run:
        coaddRunGalfit(inFile, root=root, imax=120)

    """ Optional: 2-Sersic Model """
    if ser2Comp:
        """ 2a. Generate the Read-in File for 2Ser model"""
        inFile2 = inFile.replace('1ser', '2ser')
        """ 2b. Output File """
        outFile2 = outFile.replace('1ser', '2ser')
        """ 2c. Config """
        config2Ser = copy.deepcopy(galfitConfig)
        config2Ser['output'] = outFile2
        getInput2Sersic(config2Ser, readinFile=inFile2, skyGrad=skyGrad,
            useF1=useF1, useF4=useF4)
        """ 2d. Execute the GALFIT run """
        if run:
            coaddRunGalfit(inFile2, root=root, imax=120)

    """ Optional: 3-Sersic Model """
    if ser3Comp:
        """ 3a. Generate the Read-in File for 3Ser model"""
        inFile3 = inFile.replace('1ser', '3ser')
        """ 3b. Output File """
        outFile3 = outFile.replace('1ser', '3ser')
        """ 3c. Config """
        config3Ser = copy.deepcopy(galfitConfig)
        config3Ser['output'] = outFile3
        getInput3Sersic(config3Ser, readinFile=inFile3, skyGrad=skyGrad,
            useF1=useF1, useF4=useF4)
        """ 3d. Execute the GALFIT run """
        if run:
            coaddRunGalfit(inFile3, root=root, imax=150)

    return


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="Prefix of the cutout image files")
    parser.add_argument('-r', '--root', dest='root', help='Path to the image files',
                        default=None)
    parser.add_argument('-model', dest='model', help='Suffix of the model',
                        default=None)
    parser.add_argument('-inFile', dest='inFile', help='Name of the read-in file',
                        default=None)
    parser.add_argument('-outFile', dest='outFile', help='Name of the output model',
                        default=None)
    parser.add_argument('-constrFile', dest='constrFile', help='Name of the constraint',
                        default=None)
    parser.add_argument('--pix', dest='pix', help='Pixel Scale',
                       type=float, default=0.168)
    parser.add_argument('--zp', dest='zp', help='Photometric zeropoint',
                       type=float, default=27.0)
    parser.add_argument('--mag', dest='mag', help='Total magnitude',
                       type=float, default=18.00)
    parser.add_argument('--galX0', dest='galX0', help='Galaxy Center: X',
                       type=float, default=None)
    parser.add_argument('--galY0', dest='galY0', help='Galaxy Center: Y',
                       type=float, default=None)
    parser.add_argument('--galQ0', dest='galQ0', help='Axis ratio',
                       type=float, default=None)
    parser.add_argument('--galPA0', dest='galPA0', help='Position angle',
                       type=float, default=None)
    parser.add_argument('--galRe', dest='galRe', help='Effective Radius',
                       type=float, default=None)
    parser.add_argument('--galSer', dest='galSer', help='Sersic Index',
                       type=float, default=2.0)
    parser.add_argument('--useBkg', dest='useBkg', action="store_true",
                       default=True)
    parser.add_argument('--usePsf', dest='usePsf', action="store_true",
                       default=True)
    parser.add_argument('--useSig', dest='useSig', action="store_true",
                       default=True)
    parser.add_argument('--verbose', dest='verbose', action="store_true",
                       default=True)
    parser.add_argument('--run', dest='run', action="store_true", default=True)
    parser.add_argument('--ser2Comp', dest='ser2Comp', action="store_true",
                       default=False)
    parser.add_argument('--ser3Comp', dest='ser3Comp', action="store_true",
                       default=False)
    parser.add_argument('--skyGrad', dest='skyGrad', action="store_true",
                       default=True)
    parser.add_argument('--useF1', dest='useF1', action="store_true",
                       default=False)
    parser.add_argument('--useF4', dest='useF1', action="store_true",
                       default=False)

    args = parser.parse_args()

    coaddCutoutGalfitSimple(args.prefix, root=args.root, pix=args.pix, useBkg=args.useBkg,
            zp=args.zp, usePsf=args.usePsf, galX0=args.galX0, galY0=args.galY0,
            galQ0=args.galQ0, galPA0=args.galPA0, galRe=args.galRe, galSer=args.galSer,
            model=args.model, inFile=args.inFile, outFile=args.outFile,
            useSig=args.useSig, mag=args.mag, constrFile=args.constrFile,
            verbose=args.verbose, run=args.run, ser2Comp=args.ser2Comp,
            ser3Comp=args.ser3Comp, skyGrad=args.skyGrad, useF1=args.useF1,
            useF4=args.useF4)
