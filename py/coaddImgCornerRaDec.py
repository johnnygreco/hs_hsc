#!/usr/bin/env python

from astropy.io import fits
from astropy    import wcs
import argparse
import numpy
import os

def coaddImgCornersRaDec(fileName, hdu=0):

    if not os.path.isfile(fileName):
        raise Exception('Can not find the image file: %s' % fileName)
    else:
        # Load the FITS Hdulist
        hduList = fits.open(fileName)
        # Get the header of the image HDU
        header  = hduList[hdu].header

    # Read the WCS information from the header
    w = wcs.WCS(header)
    # Get the dimension of the image
    xSize = header['NAXIS1']
    ySize = header['NAXIS2']

    # Define the Image coordinates of 4 corners
    pixLL = [0, 0]          # Lower left
    pixUL = [0, ySize]      # Upper left
    pixLR = [xSize, 0]      # Lower right
    pixUR = [xSize, ySize]  # Upper right

    # Numpy array
    pixCorners = numpy.array([pixLL, pixUL, pixLR, pixUR], numpy.float_)

    # Get the (RA, Dec) of four corners
    raDec = w.wcs_pix2world(pixCorners, 1)

    # Output results TODO: Need to improve
    print "%s , %12.7f , %12.7f , %12.7f , %12.7f , %12.7f , %12.7f , %12.7f , %12.7f" % (fileName, raDec[0][0], raDec[0][1], raDec[1][0], raDec[1][1],
               raDec[2][0], raDec[2][1], raDec[3][0], raDec[3][1])

    return raDec


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("name", help="Name of the image list file")
    args = parser.parse_args()

    if not os.path.isfile(args.name):
        raise Exception('Can not find the list file: %s' % args.name)
    else:
        imgList = open(args.name, 'r')
        for img in imgList.readlines():
            imgFile = args.root + '/' + img.strip()
            imgRaDec = coaddImgCornersRaDec(imgFile, hdu=1)

