#!/usr/bin/env python

import argparse
import numpy
import matplotlib.pyplot as pyplot
import lsst.daf.persistence as dafPersist

def main(rootDir, visit, ccd, n=5):

    # make a butler and specify your dataId
    butler = dafPersist.Butler(rootDir)
    dataId = {'visit': visit, 'ccd':ccd}

    # get the exposure from the butler
    sources = butler.get('src', dataId)

    for src in sources[0:n]:
        
        # get the coordinate object (an ICRS coord)
        icrs      = src.get('coord')
        galactic  = icrs.toGalactic()
        fk5       = icrs.toFk5()
        ecliptic  = icrs.toEcliptic()

        # get the RA and Dec.  Both are angles, but they're only defined for ICRS and Fk5
        ra, dec    = icrs.getRa(),         icrs.getDec()
        
        # the relevant type is available only when it makes sense (i.e. you can't ask an fk5 for l, and b)
        l, b       = galactic.getL(),      galactic.getB()
        ra2, dec2  = fk5.getRa(),          fk5.getDec()
        lamb, beta = ecliptic.getLambda(), ecliptic.getBeta()
        
        # you can always ask for Latitude and longitude.
        # For icrs they return RA,Dec; for galactic they return l,b; etc
        lon, lat   = icrs.getLongitude(),  icrs.getLatitude()

        sid = src.getId()

        print "ID: ", sid
        print "   ICRS        RA/Dec    (deg)", ra.asDegrees(), dec.asDegrees()
        print "   FK5         RA/Dec    (rad)", ra2.asRadians(), dec2.asRadians()
        print "   Galactic       l/b (arcmin)", l.asArcminutes(), b.asArcminutes()
        print "   Ecliptic lamb/beta (arcsec)", lamb.asArcseconds(), beta.asArcseconds()
        print "   Generic   Long/Lat    (str)", lon, lat
        
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("visit", type=int, help="Visit to show")
    parser.add_argument("ccd", type=int, help="CCD to show")
    parser.add_argument("-n", "--number", type=int, default=5, help="Number of sources to print")
    args = parser.parse_args()
    
    main(args.root, args.visit, args.ccd, n=args.number)
