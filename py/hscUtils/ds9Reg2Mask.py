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


def reg2Mask(imgFile, regFile, mskFile=None, hdu=0, show=False,
             save=True, imgHead=None):
    """
    Mask out the regions in a DS9 region file.

    Parameters:
    """
    try:
        import pyregion
    except Exception:
        raise Exception("### Please have pyregion installed first")

    if imgHead is None:
        if not os.path.isfile(imgFile):
            raise Exception("### Can not find the Image: %s" % imgFile)
        else:
            img = fits.open(imgFile)[0].data
            head = fits.open(imgFile)[0].header
    else:
        img = imgFile
        head = imgHead

    if not os.path.isfile(regFile):
        raise Exception("### Can not find the Region file: %s" % regFile)
    else:
        reg = pyregion.open(regFile).as_imagecoord(head)
    imgX, imgY = img.shape
    regMask = reg.get_mask(shape=(imgX, imgY))
    intMask = regMask.astype(int)

    if save:
        if mskFile is None:
            mskFile = regFile.replace('.reg', '_msk.fits')
        saveFits(intMask, mskFile, head=head, clobber=True)

    if show:
        pngName = mskFile.replace('.fits', '.png')
        showSEPImage(img, pngName=pngName, mask=intMask, cmap=cmap1,
                     title=mskFile.replace('.fits', ''))

    return intMask


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

    reg2Mask(args.imgFile, args.regFile, mskFile=args.mskFile, hdu=args.hdu,
             show=args.show)
