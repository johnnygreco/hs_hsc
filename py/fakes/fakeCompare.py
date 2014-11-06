#!/usr/bin/env python

import argparse
import numpy
import matplotlib.pyplot as pyplot
import lsst.daf.persistence as dafPersist
from matchFakeGalaxy import getFakeSources

def getExpArray(root, visit, ccd, filter=None):

    # make a butler and specify your dataId
    butler = dafPersist.Butler(root)
    if filter:
        dataId = {'tract': visit, 'patch':ccd, 'filter':filter}
        prefix = "deepCoadd_"
    else:
        dataId = {'visit': visit, 'ccd':int(ccd)}
        prefix = ""

    # get the exposure from the butler
    exposure = butler.get(prefix+'calexp', dataId)

    # get the maskedImage from the exposure, and the image from the mimg
    #mimg = exposure.getMaskedImage()
    #img = mimg.getImage()
    img = exposure.getImage()

    # convert to a numpy ndarray
    return img.getArray()

def getNoMatchXY(rootDir, visit, ccd):

    print visit, ccd

    # TODO: Need to be organized
    (ind, fakeXY, matchX, matchY, psfMag, psfMerr) = getFakeSources(rootDir,
                                                        visit, ccd)
    nFakes   = len(fakeXY)
    nNoMatch = (nFakes - len(numpy.argwhere(psfMag)))
    noMatchX = []
    noMatchY = []
    for i in range(nFakes):
        injectXY = fakeXY[i]
        if matchX[i] > 0:
            pass
        else:
            noMatchX.append(injectXY[0])
            noMatchY.append(injectXY[1])

    print nNoMatch, len(noMatchX)
    if len(noMatchX) is not nNoMatch:
        raise Exception("Something is wrong about the number of noMatch stars!")

    return noMatchX, noMatchY

def main(root1, root2, visit, ccd):

    # get the image array before the fake objects are added
    imgBefore = getExpArray(root1, visit, ccd)
    imgAfter  = getExpArray(root2, visit, ccd)

    # get the difference between the two image
    imgDiff = (imgAfter - imgBefore)

    # get the X, Y lists of noMatch stars
    noMatchX, noMatchY = getNoMatchXY(root2, visit, ccd)

    # stretch it with arcsinh and make a png with pyplot
    fig, axes = pyplot.subplots(1, 3, sharex=True, sharey=True, figsize=(15,10))
    pyplot.subplots_adjust(left=0.04, bottom=0.03, right=0.99, top=0.97,
                           wspace=0.01, hspace = 0.01)

    imgs   = imgBefore, imgAfter, imgDiff
    titles = "Before", "After", "Diff"
    for i in range(3):
        axes[i].imshow(numpy.arcsinh(imgs[i]), cmap='gray')
        axes[i].set_title(titles[i])
        # highlight the noMatch stars with a circle
        if i < 2:
            area = numpy.pi * 4 ** 2
            axes[i].scatter(noMatchX, noMatchY, s=area, c='r', alpha=0.5)

    pyplot.gcf().savefig("fakeCompare-%d-%s.png"%(visit,str(ccd)))


if __name__ == '__main__':

    root = '/lustre/Subaru/SSP/rerun/song/'

    parser = argparse.ArgumentParser()
    parser.add_argument("root1", help="Root directory of data before adding fake objects")
    parser.add_argument("root2", help="Root directory of data after adding fake objects")
    parser.add_argument("visit", type=int, help="Visit to show")
    parser.add_argument("ccd", type=int, help="CCD to show")
    args = parser.parse_args()

    root1 = root + args.root1
    root2 = root + args.root2

    main(root1, root2, args.visit, args.ccd)
