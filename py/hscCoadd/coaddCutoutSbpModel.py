#!/usr/bin/env python
# encoding: utf-8
"""Wrapper of BMODEL."""

from __future__ import division

import os
import gc
import copy
import string
import random
import warnings
import argparse
import subprocess
import numpy as np

# Matplotlib default settings
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Ellipse
mpl.rcParams['figure.figsize'] = 12, 10
mpl.rcParams['xtick.major.size'] = 10.0
mpl.rcParams['xtick.major.width'] = 2.5
mpl.rcParams['xtick.minor.size'] = 5.0
mpl.rcParams['xtick.minor.width'] = 2.5
mpl.rcParams['ytick.major.size'] = 10.0
mpl.rcParams['ytick.major.width'] = 2.5
mpl.rcParams['ytick.minor.size'] = 5.0
mpl.rcParams['ytick.minor.width'] = 2.5
mpl.rc('axes', linewidth=3.5)

# Astropy related
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, Column
from pyraf import iraf

# Color table
try:
    cmap = plt.get_cmap('viridis')
    cmap.set_bad('k', 1.)
except Exception:
    from palettable.cubehelix import perceptual_rainbow_16
    cmap = perceptual_rainbow_16.mpl_colormap
    cmap.set_bad('k', 1.)

# Personal
import hscUtils as hUtil

COM = '#' * 100
SEP = '-' * 100
WAR = '!' * 100


def writeBmodelPar(binary, image, backgr=0.00, interp='spline',
                   highar=False):
    """Write a parameter file for x_isophote.e bmodel."""
    outPar = binary.replace('.bin', '_bmod.par')
    outFile = binary.replace('.bin', '_bmod.fits')

    """Remove the old"""
    if os.path.isfile(outPar):
        os.remove(outPar)

    """Write the parameter file"""
    f = open(outPar, 'w')
    # ----------------------------------------------------------------- #
    f.write('\n')
    f.write('model.table = %s' % binary)
    f.write('model.output = %s' % outFile)
    f.write('model.parent = %s' % image)
    f.write('model.interp = %s' % interp)
    f.write('model.backgr = %10.5f' % backgr)
    f.write('model.fulltable = yes')
    if highar:
        f.write('model.highar = no')
    else:
        f.write('model.highar = yes')
    f.write('model.verbose = no')
    f.write('model.mode = "al"')
    f.write('# EOF')
    # ----------------------------------------------------------------- #
    f.close()

    if os.path.isfile(outPar):
        return True
    else:
        return False


def run(args):
    """Generate the residual image using Ellipse run."""
    """Call the STSDAS.ANALYSIS.ISOPHOTE package"""
    if args.isophote is None:
        if args.verbose:
            print '\n' + SEP
            print "##       Call STSDAS.ANALYSIS.ISOPHOTE() "
            print SEP
        iraf.stsdas()
        iraf.analysis()
        iraf.isophote()
        iraf.ttools()



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="Prefix of the input image")
    parser.add_argument('-r', "--root", dest='root',
                        help="Directory for the input image",
                        default=None)
    parser.add_argument('--stage', dest='stage',
                        help='Stage of Ellipse Run',
                        type=int, default=3,
                        choices=range(1, 5))
    """Optional parameters"""
    parser.add_argument('--highar', dest='highar', action="store_true",
                        help='Include a3/a4/b3/b4 ?', default=False)
    parser.add_argument("--interp", dest='interp',
                        help="Method for interpolation",
                        default='spline')
    parser.add_argument("--isophote", dest='isophote',
                        help="Location of the x_isophote.e file",
                        default=None)
    parser.add_argument("--xttools", dest='xttools',
                        help="Location of the x_ttools.e file",
                        default=None)
    parser.add_argument('--verbose', dest='verbose', action="store_true",
                        default=False)

    args = parser.parse_args()

    run(args)
