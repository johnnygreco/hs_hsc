#!/usr/bin/env python
# encoding: utf-8

import os
import copy
import argparse
import warnings
import collections
import numpy as np

# HSC Pipeline
import lsst.daf.persistence   as dafPersist
import lsst.afw.coord         as afwCoord
import lsst.afw.image         as afwImage
import lsst.afw.geom          as afwGeom
import lsst.afw.table         as afwTable

hscFilters = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y']


def getCoaddPsfImage(ra, dec, filter=None, prefix='hsc_psf', butler=None, root=None):

    """ Get the SkyMap from database """
    if butler is None:
        if root is None:
            raise Exception("### You have to either provide the Butler or Root")
        else:
            butler = dafPersist.Butler(root)
    skyMap = butler.get("deeepCoadd_skyMap", immediate=True)

    """ Define the Ra, Dec pair """
    if isinstance(ra, (int, float, numpy.number)) and isinstance(dec, (int, float,
        numpy.number)):
        point = afwGeom.Point2D(ra, dec)
        coord = afwCoord.IcrsCoord(point)
    if isinstance(ra, collections.Iterable) and isinstance(dec, collections.Iterable):
        points = map(lambda x, y: afwGeom.Point2D(x, y), ra, dec)
        coord  = map(lambda x: afwCoord.IcrsCoord(x), points)
    """ Search for overlapped tracts """
    matches = skyMap.findClosestTractPatchList([coord])

    for tract, patch in matches:

        """ TractId """
        tractId = tract.getId()
        """ PatchId: Even this (RA, DEC) is covered by more than 1 patches
            Only return 1"""
        patchId = "%d,%d" % patch[0].getIndex()
        """ Try to get the data """
        if filter is None
            filtList = hscFilters
        else:
            filtList = [filter]
        for filt in filtList:
            try:
                coadd = butler.get("deepCoadd_psf", tract=tractId,
                        patch=patchId, filter=filt, immediate=True)
            except Exception, errMsg:
                print "### No Useful data for %7.3f, %7.3f at " % (ra, dec)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("ra",   type=float, help="RA  to search")
    parser.add_argument("dec",  type=float, help="Dec to search")
    parser.add_argument("size", type=float, help="Half size of the cutout box")
    parser.add_argument('-f', '--filter', dest='filter', help="Filter",
                       default='HSC-I')
    parser.add_argument('-p', '--prefix', dest='outfile',
                        help='Prefix of the output file',
                        default='hsc_coadd_cutout')
    args = parser.parse_args()

    coaddGetPsfImage(args.root, args.ra, args.dec, args.size,
                     filter=args.filter, prefix=args.outfile)
