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

def getFakeSources(rootdir, dataId, tol=0.1):
    """Get list of sources which agree in position with fake ones with tol
    """
    butler = dafPersist.Butler(rootdir)

    sources = butler.get('src', dataId)
    cal_md = butler.get('calexp_md', dataId)


    fakeXY = collections.defaultdict(tuple)
    fakename = re.compile('FAKE([0-9]+)')
    for card in cal_md.names():
        m = fakename.match(card)
        if m is not None:
            x,y = map(float, (cal_md.get(card)).split(','))
            fakeXY[int(m.group(1))] = (x,y)


    srcX, srcY = sources.getX(), sources.getY()
    srcIndex = collections.defaultdict(list)
    for fid, fcoord  in fakeXY.items():
        matched = ((np.abs(srcX-fcoord[0]) < tol) &
                   (np.abs(srcY-fcoord[1]) < tol))
        s1 = sources.subset(matched)
        srcIndex[fid] = np.where(matched)[0]

    srcList = None
    for s in srcIndex.values():
        for ss in s:
            if srcList is None:
                srcList = SourceCatalog(sources.getSchema())
            srcList.append(sources[ss])
    return srcIndex, srcList


def main():

    #TODO: this should use the LSST/HSC conventions
    parser = argparse.ArgumentParser()
    parser.add_argument('rootDir', help='root dir of data repo')
    parser.add_argument('visit', help='id of visit', type=int)
    parser.add_argument('ccd', help='id of ccd', type=int)
    args = parser.parse_args()

    print getFakeSources(args.rootDir, {'visit':args.visit, 'ccd':args.ccd})



if __name__=='__main__':
    main()
