#!/usr/bin/env python
# encoding: utf-8

"""
Collection of useful tools to deal with galaxy images from HSC and other surveys

 * Cosmology related procedures
 * Galactic extinction related
 * Angle normalization
 * Coordinate related shortcuts

"""
import os
import numpy as np
import astropy.units as u
from astropy.utils.misc  import isiterable
from astropy.coordinates import SkyCoord

""" Path of the package"""
__path__     = os.path.realpath(__file__)
""" Path of the data directory"""
__dataDir__  = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
if not os.path.exists(__dataDir__):
    raise Exception("Can not find the data directory: %s" % __dataDir__)

"""
Angle related functions:

Strongly based on: https://github.com/phn/angles/blob/master/angles.py
by Prasanth Nair

"""

def rad2deg(rad):
    """ Convert radians into degrees"""
    return (rad * 180.0 / np.pi)

def deg2rad(deg):
    """ Convert degrees into radians"""
    return (deg * np.pi / 180.0)

def hr2deg(deg):
    """ Convert degrees into hours"""
    return (deg *(24.0 / 360.0))

def deg2hr(hr):
    """ Convert hours into degrees"""
    return (hr * 15.0)

def normAngle(num, lower=0, upper=360, b=False):
    """Normalize number to range [lower, upper) or [lower, upper].
    Parameters
    ----------
    num : float
        The number to be normalized.
    lower : int
        Lower limit of range. Default is 0.
    upper : int
        Upper limit of range. Default is 360.
    b : bool
        Type of normalization. Default is False. See notes.
    Returns
    -------
    n : float
        A number in the range [lower, upper) or [lower, upper].
    Raises
    ------
    ValueError
      If lower >= upper.
    Notes
    -----
    If the keyword `b == False`, then the normalization is done in the
    following way. Consider the numbers to be arranged in a circle,
    with the lower and upper ends sitting on top of each other. Moving
    past one limit, takes the number into the beginning of the other
    end. For example, if range is [0 - 360), then 361 becomes 1 and 360
    becomes 0. Negative numbers move from higher to lower numbers. So,
    -1 normalized to [0 - 360) becomes 359.
    If the keyword `b == True`, then the given number is considered to
    "bounce" between the two limits. So, -91 normalized to [-90, 90],
    becomes -89, instead of 89. In this case the range is [lower,
    upper]. This code is based on the function `fmt_delta` of `TPM`.
    Range must be symmetric about 0 or lower == 0.
    Examples
    --------
    >>> normalize(-270,-180,180)
    90.0
    >>> import math
    >>> math.degrees(normalize(-2*math.pi,-math.pi,math.pi))
    0.0
    >>> normalize(-180, -180, 180)
    -180.0
    >>> normalize(180, -180, 180)
    -180.0
    >>> normalize(180, -180, 180, b=True)
    180.0
    >>> normalize(181,-180,180)
    -179.0
    >>> normalize(181, -180, 180, b=True)
    179.0
    >>> normalize(-180,0,360)
    180.0
    >>> normalize(36,0,24)
    12.0
    >>> normalize(368.5,-180,180)
    8.5
    >>> normalize(-100, -90, 90)
    80.0
    >>> normalize(-100, -90, 90, b=True)
    -80.0
    >>> normalize(100, -90, 90, b=True)
    80.0
    >>> normalize(181, -90, 90, b=True)
    -1.0
    >>> normalize(270, -90, 90, b=True)
    -90.0
    >>> normalize(271, -90, 90, b=True)
    -89.0
    """
    from math import floor, ceil
    # abs(num + upper) and abs(num - lower) are needed, instead of
    # abs(num), since the lower and upper limits need not be 0. We need
    # to add half size of the range, so that the final result is lower +
    # <value> or upper - <value>, respectively.
    res = num
    if not b:
        if lower >= upper:
            raise ValueError("Invalid lower and upper limits: (%s, %s)" %
                             (lower, upper))

        res = num
        if num > upper or num == lower:
            num = lower + abs(num + upper) % (abs(lower) + abs(upper))
        if num < lower or num == upper:
            num = upper - abs(num - lower) % (abs(lower) + abs(upper))

        res = lower if num == upper else num
    else:
        total_length = abs(lower) + abs(upper)
        if num < -total_length:
            num += ceil(num / (-2 * total_length)) * 2 * total_length
        if num > total_length:
            num -= floor(num / (2 * total_length)) * 2 * total_length
        if num > upper:
            num = total_length - num
        if num < lower:
            num = -total_length - num

        res = num

    res *= 1.0  # Make all numbers float, to be consistent

    return res


"""
Coordinate related shortcuts

    * Convert from (ra, dec) to (l, b)
    * Conversion between ICRS and FK5
"""

def radec2lb(ra, dec, radian=False, FK5=False):

    """
    Convert (ra, dec) into Galactic coordinate (l, b)

    Parameters
    ----------
    ra : float or list or array
        RA Coordinates in degree
    dec : float or list or array
        DEC Coordinates in degree

    Returns
    -------
    l : float or list or array
    b : float or list or array
    """

    """ See if the input is number or array"""
    if not (isiterable(ra) or isiterable(dec)):
        returnScalar = True
        if not FK5:
            raDec = [SkyCoord(ra, dec, frame='icrs', unit='deg')]
        else:
            raDec = [SkyCoord(ra, dec, frame='fk5', unit='deg')]
    else:
        returnScalar = False
        if not FK5:
            raDec = [SkyCoord(ra, dec, frame='icrs', unit='deg')
                     for rrr, ddd in zip(ra, dec)]
        else:
            raDec = [SkyCoord(ra, dec, frame='fk5', unit='deg')
                     for rrr, ddd in zip(ra, dec)]

    """ Convert to galactic coordinates
        Currently, coordinates do not support arrays; have to loop.
    """
    l = np.empty(len(raDec), dtype=np.float)
    b = np.empty(len(raDec), dtype=np.float)

    for ii, cc in enumerate(raDec):
        gg = cc.galactic
        # Hack to support both astropy v0.2.4 and v0.3.dev
        # TODO: remove this hack once v0.3 is out (and array-ify this
        # whole thing)
        if radian:
            l[ii] = gg.l.radian
            b[ii] = gg.b.radian
        else:
            l[ii] = gg.l.degree
            b[ii] = gg.b.degree

    if returnScalar:
        return l[0], b[0]
    else:
        return l, b


def icrs2fk5(ra, dec, radian=False):

    """ Convert coordinates from ICRS to FK5 frame """
    if not radian:
        raDec = SkyCoord(ra, dec, frame='icrs', unit='deg')
    else:
        raDec = SkyCoord(ra, dec, frame='icrs', unit='radian')

    raDecFK5 = raDec.transform_to('fk5')

    if not radian:
        return raDecFK5.ra.degree, raDecFK5.dec.degree
    else:
        return raDecFK5.ra.radian, raDecFK5.dec.radian

def fk52icrs(ra, dec, radian=False):

    """ Convert coordinates from FK5 to ICRS frame """
    if not radian:
        raDec = SkyCoord(ra, dec, frame='fk5', unit='deg')
    else:
        raDec = SkyCoord(ra, dec, frame='fk5', unit='radian')

    raDecFK5 = raDec.transform_to('icrs')

    if not radian:
        return raDecFK5.ra.degree, raDecFK5.dec.degree
    else:
        return raDecFK5.ra.radian, raDecFK5.dec.radian


"""
Image Visualization Related

    * zScale of image
"""

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





