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


def srcToRaDec(src):
    """
    Return a list of (RA, DEC)
    """

    raArr  = np.asarray([coord[0] * 180.0 / np.pi for coord in src['coord']])
    decArr = np.asarray([coord[1] * 180.0 / np.pi for coord in src['coord']])

    return raArr, decArr


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
