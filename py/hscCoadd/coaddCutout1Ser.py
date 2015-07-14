#!/usr/bin/env python
# encoding: utf-8

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
    imaxStr = " -imax %4d", imax

    """ GALFIT command """
    galfitCommand = galfit + ' ' + imaxStr + ' ' + readFile

    """ Excecute the command """
    proc = subprocess.Popen([galfitCommand], cwd=root, shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in proc.stdout.readlines():
        print line
    retval = proc.wait()

    return


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


def coaddCutout1Ser(prefix, root=None):

    """
    Run 1-Sersic fitting on HSC cutout image
    """
