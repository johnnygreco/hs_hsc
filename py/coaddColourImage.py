#!/usr/bin/env python

import argparse
import lsst.daf.persistence as dafPersist
import lsst.afw.display.rgb as afwRgb
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.image as afwImage

def coaddColourImage(root, ra, dec, size, filt='gri',
                     prefix='hsc_coadd_cutout', info=None,
                     min=-0.0, max=0.3, Q=8):

    # Get the SkyMap of the database
    butler = dafPersist.Butler(root)
    skyMap = butler.get("deepCoadd_skyMap", immediate=True)

    # [Ra, Dec] pair
    raDec = afwCoord.Coord(ra*afwGeom.degrees, dec*afwGeom.degrees)

    # Check the choice of filters
    hscFilters = ['g', 'r', 'i', 'z', 'y']
    if len(filt) is not 3:
        raise Exception("Have to be three filters!")
    elif filt[0] not in hscFilters:
        raise Exception("%s is not a valid HSC filter!" % filt[0])
    elif filt[1] not in hscFilters:
        raise Exception("%s is not a valid HSC filter!" % filt[1])
    elif filt[2] not in hscFilters:
        raise Exception("%s is not a valid HSC filter!" % filt[2])

    # Define the name of the output file
    outRgb = prefix + '_' + filt + '_color.png'

    #Figure out the area we want, and read the data.
    #For coadds the WCS is the same in all bands, but the code handles the general case
    #Start by finding the tract and patch
    tractInfo, patchInfo = skyMap.findTractPatchList([raDec])[0]
    tract = tractInfo.getId()
    patch = "%d,%d" % patchInfo[0].getIndex()

    # Check if the coordinate is available in all three bands.
    try:

        # Get the correct HSC filter name
        filter1 = "HSC-%s" % filt[0].upper()
        filter2 = "HSC-%s" % filt[1].upper()
        filter3 = "HSC-%s" % filt[2].upper()

        # Get the metadata
        md1 = butler.get("deepCoadd_md", immediate=True,
                        tract=tract, patch=patch, filter=filter1)
        md2 = butler.get("deepCoadd_md", immediate=True,
                        tract=tract, patch=patch, filter=filter2)
        md3 = butler.get("deepCoadd_md", immediate=True,
                        tract=tract, patch=patch, filter=filter3)

        filtArr = [filter1, filter2, filter3]

    except Exception, errMsg:

        colorFound = False
        print "#############################################"
        print " The desired coordinate is not available !!! "
        print "#############################################"
        print errMsg

    else:

        #Then we can read the desired pixels
        images = {}
        cutoutSize = int(size)

        for i in range(3):

            # Find the file of the coadd image
            coadd = butler.get("deepCoadd", immediate=True,
                                  tract=tract, patch=patch, filter=filtArr[i])

            # Get the WCS information
            wcs = coadd.getWcs()

            # Convert the central coordinate from Ra,Dec to pixel unit
            pixel = wcs.skyToPixel(raDec)
            pixel = afwGeom.Point2I(pixel)

            # Define the bounding box for the central pixel
            bbox = afwGeom.Box2I(pixel, pixel)

            # Grow the bounding box to the desired size
            bbox.grow(int(cutoutSize))
            bbox.clip(coadd.getBBox(afwImage.PARENT))

            # Get the masked image
            subImage  = afwImage.ExposureF(coadd, bbox, afwImage.PARENT)
            images[i] = subImage.getMaskedImage()

        # Define the Blue, Green, and Red channels
        B, G, R = images[0].getImage(), images[1].getImage(), images[2].getImage()

        # Generate the RGB image
        imgRgb = afwRgb.makeRGB(R, G, B, min=min, range=(max - min), Q=Q,
                               saturatedPixelValue=None)

        #afwRgb.displayRGB(rgb)

        # Save it to a .PNG file
        afwRgb.writeRGB(outRgb, imgRgb)

        colorFound = True

    return colorFound


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("ra",   type=float, help="RA  to search")
    parser.add_argument("dec",  type=float, help="Dec to search")
    parser.add_argument("size", type=float, help="Half size of the cutout box")
    parser.add_argument('-f', '--filters', dest='filt',
                        help="Combination of three filters for the colour image",
                        default='gri')
    parser.add_argument('-p', '--prefix', dest='outfile',
                        help='Prefix of the output file',
                        default='hsc_coadd_cutout')
    parser.add_argument('-i', '--info', dest='info',
                        help='Information to show on the image',
                        default=None)
    args = parser.parse_args()

    coaddColourImage(args.root, args.ra, args.dec, args.size,
                     filt=args.filt, prefix=args.outfile,
                     info=args.info)
