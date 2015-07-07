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

import scipy.interpolate
import scipy.ndimage

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


"""
Image Manipulation

    * Image resampling, like the rebin function in IDL
    * Image rotation
    * Image shift (sub-pixel accurarcy)
"""
def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''
    From: http://wiki.scipy.org/Cookbook/Rebinning
    Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).

    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [np.float64, np.float32]:
        a = np.cast[float](a)

    m1 = np.cast[int](minusone)
    ofs = np.cast[int](centre) * 0.5
    old = np.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print "[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions."
        return None
    newdims = np.asarray( newdims, dtype=float )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = np.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = np.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = np.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [np.arange(i, dtype = np.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = np.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = np.mgrid[nslices]

        newcoords_dims = range(np.rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (np.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print "Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported."
        return None


"""
File Manipulation

    * Save numpy array to cPickle file format
    * Save numpy array to hickle/HDF5 format
    * Save numpy array to csv file format
"""
def saveToPickle(array, name):
    """
    Save a numpy array to a cPickle/Pickle format binary file
    """
    try:
       import cPickle as pickle
    except:
       import pickle

    output= open(name, 'w')
    pickle.dump(array, output, protocol=2)
    output.close()


def saveToHickle(array, name):
    """
    Save a numpy array to a hickle/HDF5 format binary file
    """
    try:
       import hickle
    except:
        raise Exception("### The Hickle package is required!")

    output= open(name, 'w')
    hickle.dump(array, output, protocol=2)
    output.close()


def saveToCSV(array, name):
    """
    Save a numpy array to a CSV file

    Use the dtype.name as column name if possible
    """
    output= open(name, 'w')
    colNames = array.dtype.names
    output.write("#" + ', '.join(colNames) + '\n')
    for item in array:
        line = ''
        for i in range(0, len(colNames)-1):
            col = colNames[i]
            line += str(item[col]) + ' , '
        line += str(item[colNames[-1]]) + '\n'
        output.write(line)
    output.close()

def parseRegEllipse(regName):
    """
    Parse a DS9 .reg files, convert the Ellipse or Circle regions
    into arrays of parameters for ellipse:
    x, y, a, b, theta
    """
    if os.path.isfile(regName):
        raise Exception("### Can not find the .reg file!")
    # Parse the .reg file into lines
    lines = [line.strip() for line in open(regName, 'r')]
    # Coordinate type of this .reg file: e.g. 'image'
    coordType = lines[2].strip()
    # Parse each region
    regs = [reg.split(" ") for reg in lines[3:]]

    xc = []
    yc = []
    ra = []
    rb = []
    theta = []

    for reg in regs:
        if reg[0].strip() == 'ellipse' and len(reg) is 6:
            xc.append(float(reg[1]))
            yc.append(float(reg[2]))
            ra.append(float(reg[3]))
            rb.append(float(reg[4]))
            theta.append(float(reg[5]) * np.pi / 180.0)
        elif reg[0].strip() == 'circle' and len(reg) is 4:
            xc.append(float(reg[1]))
            yc.append(float(reg[2]))
            ra.append(float(reg[3]))
            rb.append(float(reg[3]))
            theta.append(0.0)

    xc = np.array(xc, dtype=np.float32)
    yc = np.array(yc, dtype=np.float32)
    ra = np.array(ra, dtype=np.float32)
    rb = np.array(rb, dtype=np.float32)
    theta = np.array(theta, dtype=np.float32)

    return xc, yc, ra, rb, theta, coordType


"""
Cosmology Related

    * Get luminosity distance at redshift=z
    * Get angular diameter distance at redshift = z
    * Get pixel scale at redshift = z
    * Get distance module at redshift = z
    * Get comoving volume at redshift = z
    * Get differential volume at redsfhit = z
    * Get the age of the Universe at redshift = z
    * Get the look-back time at redshift = z

"""

def cosmoDL(redshift, WMAP9=True, H0=69.3, Om0=0.287,
             Planck15=True, kpc=False):
    """
    Get the Luminosity Distance at redshift=z

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    dl = cosmo.luminosity_distance(redshift)

    if not kpc:
        return dl.value
    else:
        return dl.to(u.kpc).value


def cosmoDA(redshift, WMAP9=True, H0=69.3, Om0=0.287,
             Planck15=True, kpc=False):
    """
    Get the Angular Diameter Distance at redshift=z

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    da = cosmo.angular_diameter_distance(redshift)

    if not kpc:
        return da.value
    else:
        return da.to(u.kpc).value


def cosmoScale(redshift, WMAP9=True, H0=69.3, Om0=0.287,
                Planck15=True):
    """
    Get the Angular Scale (kpc/") at redshift=z

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    scale = cosmo.kpc_proper_per_arcmin(redshift).to(u.kpc / u.arcsec)

    return scale.value


def cosmoDistMod(redshift, WMAP9=True, H0=69.3, Om0=0.287,
                  Planck15=True):
    """
    Get the Distance Module at redshift=z

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    dm = cosmo.distmod(redshift)

    return dm.value


def cosmoComVol(redshift, WMAP9=True, H0=69.3, Om0=0.287,
                 Planck15=True, Gpc=False):
    """
    Get the Comoving Volume at redshift=z

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    v = cosmo.comoving_volume(redshift)

    if not Gpc:
        return v.value
    else:
        return v.to(u.Gpc).value


def cosmodVol(redshift, WMAP9=True, H0=69.3, Om0=0.287,
              Planck15=True):
    """
    Get the Differential Comoving Volume at redshift=z

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    dv = cosmo.differential_comoving_volume(redshift)

    return dv.value


def cosmoAge(redshift, WMAP9=True, H0=69.3, Om0=0.287,
             Planck15=True, Myr=False):
    """
    Get the Age of the Universe at redshift=z

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    age = cosmo.age(redshift)

    if not Myr:
        return age.value
    else:
        return age.to(u.Myr).value


def cosmoLookBack(redshift, WMAP9=True, H0=69.3, Om0=0.287,
                 Planck15=True, Myr=False):
    """
    Get the Look-back Time at redshift=z

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    lbt = cosmo.lookback_time(redshift)

    if not Myr:
        return lbt.value
    else:
        return lbt.to(u.Myr).value


"""
Galactic Extinction Related

"""


"""
Geometry Related
"""

def ellipDist(x, y, x0, y0, pa=0.0, q=0.9):
    """
    doc
    """
    theta = (pa * np.pi / 180.0)

    distA = ((x - x0) * np.cos(theta) + (y - y0) * np.sin(theta)) ** 2.0
    distB = (((y - y0) * np.cos(theta) - (x - x0) * np.sin(theta)) / q) ** 2.0

    return np.sqrt(distA + distB)

