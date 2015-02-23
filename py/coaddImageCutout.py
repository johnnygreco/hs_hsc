#!/usr/bin/env python

import argparse
import numpy                  as np
import lsst.daf.persistence   as dafPersist
import lsst.afw.coord         as afwCoord
import lsst.afw.image         as afwImage
import lsst.afw.geom          as afwGeom


def coaddImageCutout(root, ra, dec, filt, prefix):

    butler = dafPersist.Butler(root)
    skyMap = butler.get("deepCoadd_skyMap", immediate=True)

    point = afwGeom.Point2D(ra, dec)
    coord = afwCoord.IcrsCoord(point)

    for tract, patch in skyMap.findClosestTractPatchList([coord]):
        tractId = tract.getId()
        patchId = "%d,%d" % patch.getIndex()
        coadd = butler.get("deepCoadd", tract=tractId,
                           patch=patchId, filter=filt, immediate=True)

        pixel = coadd.getWcs().skyToPixel(coord)
        bbox = afwGeom.Box2I(pixel, pixel)
        bbox.grow(300)
        bbox.clip(coadd.getBBox(afwImage.PARENT))
        if bbox.isEmpty():
            continue
        subImage = afwImage.ExposureF(coadd, bbox, afwImage.PARENT)

        outFits = prefix + str(tractId) + '_' + patchId + '_' + filt + '.fits'
        subImage.writeFits(outFits)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("ra",   type=float, help="RA  to search")
    parser.add_argument("dec",  type=float, help="Dec to search")
    parser.add_argument('-f', '--filter', dest='filt', help="Filter",
                       default='HSC-I')
    parser.add_argument('-p', '--prefix', dest='outfile',
                        help='Prefix of the output file',
                        default='hsc_coadd_cutout')
    args = parser.parse_args()

    coaddImageCutout(args.root, args.ra, args.dec, args.filt, args.outfile)
