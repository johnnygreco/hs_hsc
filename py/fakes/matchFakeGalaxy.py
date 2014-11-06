#!/usr/bin/env python
"""
matchFakes.py
matches fakes based on position stored in the calibrated exposure image header
"""

import lsst.daf.persistence as dafPersist
from lsst.afw.table.tableLib import SourceCatalog
import numpy as np
import argparse
import re
import collections

def getFakeSources(rootdir, visit, ccd, tol=2.0):
    """Get list of sources which agree in position with fake ones with tol
    """
    butler = dafPersist.Butler(rootdir)
    visit  = int(visit)
    ccd    = int(ccd)
    dataId = {'visit':visit, 'ccd':ccd}

    sources = butler.get('src', dataId)
    cal_md  = butler.get('calexp_md', dataId)

    # Get the X, Y locations of objects on the CCD
    srcX, srcY = sources.getX(), sources.getY()
    # Get the zeropoint
    zeropoint = 2.5*np.log10(cal_md.get("FLUXMAG0"))
    # Get the PSF flux and its error
    flux, ferr = sources.get('cmodel.flux'), sources.get('cmodel.flux.err')
    # Convert them into magnitude and its error
    mag,  merr = 2.5*np.log10(flux), 2.5/np.log(10)*(ferr/flux)
    mag = zeropoint - mag

    # X, Y locations of the fake stars
    fakeXY   = collections.defaultdict(tuple)
    fakeCount = 0
    # Regular Expression
    #fakename = re.compile('FAKE([0-9]+)')
    fakename = re.compile(r"FAKE*")
    for card in cal_md.names():
        m = fakename.match(card)
        if m is not None:
            x,y = map(float, (cal_md.get(card)).split(','))
            #fakeXY[int(m.group(1))] = (x,y)
            fakeXY[int(fakeCount)] = (x,y)
            fakeCount += 1

    srcIndex = collections.defaultdict(list)
    for fid, fcoord  in fakeXY.items():
        matched = ((np.abs(srcX-fcoord[0]) < tol) &
                   (np.abs(srcY-fcoord[1]) < tol))
        s1 = sources.subset(matched)
        srcIndex[fid] = np.where(matched)[0]

    #srcList    = None
    srcPsfMag  = []
    srcPsfMerr = []
    matchX     = []
    matchY     = []
    for s in srcIndex.values():
        #for ss in s:
        #if srcList is None:
        #   srcList = SourceCatalog(sources.getSchema())
        #   srcList.append(sources[ss])
        #
        if len(s) > 0:
            ss = s[0]
            srcPsfMag.append(mag[ss])
            srcPsfMerr.append(merr[ss])
            matchX.append(srcX[ss])
            matchY.append(srcY[ss])
        else:
            srcPsfMag.append(0)
            srcPsfMerr.append(0)
            matchX.append(0)
            matchY.append(0)

    return srcIndex, fakeXY, matchX, matchY, srcPsfMag, srcPsfMerr


def main():

    #TODO: this should use the LSST/HSC conventions
    parser = argparse.ArgumentParser()
    parser.add_argument('rootDir', help='root dir of data repo')
    parser.add_argument('visit', help='id of visit', type=int)
    parser.add_argument('ccd', help='id of ccd', type=int)
    args = parser.parse_args()

    visit = int(args.visit)
    ccd   = int(args.ccd)

    #(starIndex,starList) = getFakeSources(args.rootDir, {'visit':args.visit, 'ccd':args.ccd})
    (starIndex, fakeXY, matchX, matchY, starPsfMag, starPsfMerr) = getFakeSources(args.rootDir,
                                                                    visit, ccd)

    nInject = len(fakeXY)
    nMatch  = len(np.argwhere(starPsfMag))
    print "# Number of Injected Stars : %d" % nInject
    print "# Number of Matched  Stars : %d" % nMatch
    print "# Visit = %d   CCD = %d" % (args.visit, args.ccd)
    print "# FakeX  FakeY  PSFMag  PSFMagErr  Deblend "

    for i in range(nInject):
       if len(starIndex[i]) > 1:
           deblend = "blended"
       elif matchX[i] > 0:
           deblend = "isolate"
       else:
           deblend = "nomatch"

       injectXY = fakeXY[i]

       print "%6.1d   %6.1d   %7.3f  %6.3f  %s" % (injectXY[0], injectXY[1],
                                            starPsfMag[i], starPsfMerr[i], deblend)


if __name__=='__main__':
    main()
