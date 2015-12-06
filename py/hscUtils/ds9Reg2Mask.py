#!/usr/bin/env python
# encoding: utf-8

from __future__ import division

import os
import argparse

from astropy.io import fits
import cubehelix
# For high-contrast image
cmap1 = cubehelix.cmap(start=0.5, rot=-0.8, gamma=1.0,
                       minSat=1.2, maxSat=1.2,
                       minLight=0.0, maxLight=1.0)
cmap1.set_bad('k', 1.)

from coaddCutoutPrepare import saveFits, showSEPImage


def run(imgFile, regFile, mskFile=None, hdu=0, show=False):
    """
    Mask out the regions in a DS9 region file.

    Parameters:
    """
    try:
        import pyregion
        # import pyregion._region_filter as regFilter
    except Exception:
        raise Exception("### Please have pyregion installed first")

    if not os.path.isfile(imgFile):
        raise Exception("### Can not find the Image: %s" % imgFile)
    else:
        print imgFile
        img = fits.open(imgFile)[0].data
        head = fits.open(imgFile)[0].header

    if not os.path.isfile(regFile):
        raise Exception("### Can not find the Region file: %s" % regFile)
    else:
        print "## Load in the region file : %s" % regFile
        reg = pyregion.open(regFile).as_imagecoord(head)

    print "## Get the region mask"
    imgX, imgY = img.shape
    regMask = reg.get_mask(shape=(imgX, imgY))
    print "## Convert the mask into an int array"
    intMask = regMask.astype(int)

    if mskFile is None:
        mskFile = regFile.replace('.reg', '_msk.fits')
    print "## Save the regMask to %s" % mskFile
    saveFits(intMask, mskFile, head=head, clobber=True)

    if show:
        pngName = mskFile.replace('.fits', '.png')
        showSEPImage(img, pngName=pngName, mask=intMask, cmap=cmap1,
                     title=mskFile.replace('.fits', ''))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("imgFile", help="Name of the image")
    parser.add_argument("regFile", help="Name of the DS9 region")
    parser.add_argument("-m", "--mskFile", dest='mskFile', default=None,
                        help="Name of the DS9 region")
    parser.add_argument("-s", "--show", dest='show', action="store_true",
                        default=False,
                        help="Whether to show the mask in PNG figure")
    parser.add_argument("--hdu", dest='hdu', default=0,
                        help="The HDU to be used", type=int)
    args = parser.parse_args()

    run(args.imgFile, args.regFile, mskFile=args.mskFile, hdu=args.hdu,
        show=args.show)
