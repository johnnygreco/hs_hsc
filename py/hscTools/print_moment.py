#!/usr/bin/env python

import argparse
import numpy
import matplotlib.pyplot as pyplot
import lsst.daf.persistence as dafPersist
import lsst.afw.geom.ellipses as geomEllip

def main(rootDir, visit, ccd, n=5):

    # make a butler and specify your dataId
    butler = dafPersist.Butler(rootDir)
    dataId = {'visit': visit, 'ccd':ccd}

    # get the exposure from the butler
    sources = butler.get('src', dataId)

    for src in sources[0:n]:

        # get the adaptive moments
        m = src.get('shape.sdss')
        ixx, ixy, iyy = m.getIxx(), m.getIxy(), m.getIyy()

        # convert to an ellipse (note that theta is in radians and is not an Angle object)
        e  = geomEllip.Axes(m)
        a, b, theta = e.getA(), e.getB(), e.getTheta()
        print ixx, ixy, iyy, "  ", a, b, theta

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("visit", type=int, help="Visit to show")
    parser.add_argument("ccd", type=int, help="CCD to show")
    parser.add_argument("-n", "--number", type=int, default=5, help="Number of sources to print")
    args = parser.parse_args()
    
    main(args.root, args.visit, args.ccd, n=args.number)
