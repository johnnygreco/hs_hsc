#!/usr/bin/env python

import argparse
import numpy
import matplotlib.pyplot as pyplot
import lsst.daf.persistence as dafPersist

def main(rootDir, visit, ccd):

    # make a butler and specify your dataId
    butler = dafPersist.Butler(rootDir)
    dataId = {'visit': visit, 'ccd':ccd}

    # get the exposure from the butler
    exposure = butler.get('calexp', dataId)

    # get the maskedImage from the exposure, and the image from the mimg
    mimg = exposure.getMaskedImage()
    img = mimg.getImage()

    # convert to a numpy ndarray
    nimg = img.getArray()

    # stretch it with arcsinh and make a png with pyplot
    pyplot.imshow(numpy.arcsinh(nimg), cmap='gray')
    pyplot.gcf().savefig("test.png")
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("visit", type=int, help="Visit to show")
    parser.add_argument("ccd", type=int, help="CCD to show")
    args = parser.parse_args()
    
    main(args.root, args.visit, args.ccd)
