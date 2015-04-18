#!/usr/bin/env python

import argparse
import numpy
import matplotlib.pyplot      as pyplot
import lsst.daf.persistence   as dafPersist
import lsst.afw.geom          as afwGeom
import lsst.afw.image         as afwImage
import lsst.afw.display.utils as dispUtils
import lsst.afw.display.ds9   as ds9

def zscale(img, contrast=0.25, samples=500):

    ravel = img.ravel()
    if len(ravel) > samples:
        imsort = numpy.sort(numpy.random.choice(ravel, size=samples))
    else:
        imsort = numpy.sort(ravel)

    n = len(imsort)
    idx = numpy.arange(n)

    med = imsort[n/2]
    w = 0.25
    i_lo, i_hi = int((0.5-w)*n), int((0.5+w)*n)
    p = numpy.polyfit(idx[i_lo:i_hi], imsort[i_lo:i_hi], 1)
    slope, intercept = p

    z1 = med - (slope/contrast)*(n/2-n*w)
    z2 = med + (slope/contrast)*(n/2-n*w)

    return z1, z2
    

def main(rootDir, visit, ccd, n=4):

    # make a butler and specify your dataId
    butler = dafPersist.Butler(rootDir)
    dataId = {'visit': visit, 'ccd':ccd}

    # get the souces and the exposure from the butler
    sources = butler.get('src', dataId)
    exposure = butler.get("calexp", dataId)

    # make a list of the ones with the calib.psf.used flag set
    psfSources = [s for s in sources if s.get("calib.psf.used")]

    wx, wy = 25, 25         # the width to use for the bounding box

    subImages = []
    labels    = []
    for psfSrc in psfSources[0:n]:
        
        # make a bounding box for this source
        x, y = int(psfSrc.getX()), int(psfSrc.getY())
        bbox = afwGeom.Box2I(afwGeom.Point2I(x-wx/2, y-wy/2), afwGeom.Extent2I(wx, wy))

        # create a new image, and add it to our list
        subimg = afwImage.ImageF(exposure.getMaskedImage().getImage(), bbox)
        subImages.append(subimg)

        # use the ID for a label
        label = "ID="+str(psfSrc.getId())
        labels.append(label)


    # create a Mosaic object and set the specifications for your mosaic
    m = dispUtils.Mosaic()
    m.setGutter(2)          # width of the space between each subimage
    m.setBackground(0)
    m.setMode("square")

    # create the mosaic
    mosaic = m.makeMosaic(subImages)

    # display it with labels in ds9
    ds9.mtv(mosaic)                
    m.drawLabels(labels)

    # make a pyplot png (note: y-axis flipped for pyplot)
    img = mosaic.getArray()[::-1,:]
    vmin, vmax = zscale(img)
    pyplot.imshow(img, interpolation="none", cmap="gray", vmin=vmin, vmax=vmax)
    pyplot.gcf().savefig("mosaic.png")

    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("visit", type=int, help="Visit to show")
    parser.add_argument("ccd", type=int, help="CCD to show")
    parser.add_argument("-n", "--number", type=int, default=4, help="Number of sources to show")
    args = parser.parse_args()
    
    main(args.root, args.visit, args.ccd, n=args.number)
