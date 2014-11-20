#!/usr/bin/env python

import argparse
import numpy
import matplotlib.pyplot as pyplot
import lsst.daf.persistence as dafPersist

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
    mimg = exposure.getMaskedImage()
    img = mimg.getImage()

    # convert to a numpy ndarray
    return img.getArray()


def showFakeObjects(root1, root2, visit, ccd, root="", matchObjs=None,
                    noMatch=None, badMatch=None):

    # get the image array before the fake objects are added
    imgBefore = getExpArray(root + root1, visit, ccd)
    imgAfter  = getExpArray(root + root2, visit, ccd)

    # get the difference between the two image
    imgDiff = (imgAfter - imgBefore)

    # stretch it with arcsinh and make a png with pyplot
    fig, axes = pyplot.subplots(1, 3, sharex=True, sharey=True, figsize=(15,10))
    pyplot.subplots_adjust(left=0.04, bottom=0.03, right=0.99, top=0.97,
                           wspace=0.01, hspace = 0.01)

    imgs   = imgBefore, imgAfter, imgDiff
    titles = "Before", "After", "Diff"
    for i in range(3):
        axes[i].imshow(numpy.arcsinh(imgs[i]), cmap='gray')
        axes[i].set_title(titles[i])

        area1 = numpy.pi * 6 ** 2
        area2 = numpy.pi * 4 ** 2

        if matchObjs is not None:
            axes[i].scatter(matchObjs['X'], matchObjs['Y'], s=area1,
                            edgecolors='g', alpha=0.9)
        if noMatch is not None:
            axes[i].scatter(noMatch['X'], noMatch['Y'], s=area2, c='r',
                            alpha=0.3)
        if badMatch is not None:
            axes[i].scatter(badMatch['X'], badMatch['Y'], s=area2, c='b',
                            alpha=0.4)

    pyplot.gcf().savefig("%s-%d-%s.png"%(root2, visit, str(ccd)))


def plotParamDiff(catMatch, paramIn, paramOut):

    # TODO


    return None


def main(root1, root2, visit, ccd, cat, root=""):

    # TODO


    return None


if __name__ == '__main__':

    root = '/lustre/Subaru/SSP/rerun/song/'

    parser = argparse.ArgumentParser()
    parser.add_argument("root1", help="Root directory of data before adding fake objects")
    parser.add_argument("root2", help="Root directory of data after adding fake objects")
    parser.add_argument("visit", type=int, help="Visit to show")
    parser.add_argument("ccd", type=int, help="CCD to show")
    parser.add_argument("cat", type=int, help="input FITS catalog")
    args = parser.parse_args()

    main(args.root1, args.root2, args.visit, args.ccd, args.cat, root=root)
