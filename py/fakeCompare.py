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


def main(root1, root2, visit, ccd, filter=None):

    # get the image array before the fake objects are added
    imgBefore = getExpArray(root1, visit, ccd, filter=None)
    imgAfter  = getExpArray(root2, visit, ccd, filter=None)

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
    pyplot.gcf().savefig("fakeCompare-%d-%s.png"%(visit,str(ccd)))


if __name__ == '__main__':

    root = '/lustre/Subaru/SSP/rerun/song/'

    parser = argparse.ArgumentParser()
    parser.add_argument("root1", help="Root directory of data before adding fake objects")
    parser.add_argument("root2", help="Root directory of data after adding fake objects")
    parser.add_argument("visit", type=int, help="Visit to show")
    parser.add_argument("ccd", type=str, help="CCD to show")
    parser.add_argument("-f", "--filter", help="Filter (will cause visit/ccd to be read as tract/patch)")
    args = parser.parse_args()

    root1 = root + args.root1
    root2 = root + args.root2

    main(root1, root2, args.visit, args.ccd, filter=args.filter)
