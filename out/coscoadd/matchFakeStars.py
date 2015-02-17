#!/usr/bin/env python
"""
matchFakes.py
matches fakes based on position stored in the calibrated exposure image header
"""

import lsst.daf.persistence as dafPersist
from lsst.afw.table import SourceCatalog, SchemaMapper
import numpy as np
import argparse
import re
import collections

def getStars(rootdir, visit, ccd, tol):
    """Get list of sources which agree in position with fake ones with tol
    """
    # Call the butler
    butler = dafPersist.Butler(rootdir)
    dataId = {'visit':visit, 'ccd':ccd}
    tol = float(tol)

    # Get the source catalog and metadata
    sources = butler.get('src', dataId)
    cal_md  = butler.get('calexp_md', dataId)

    # Get the X, Y locations of objects on the CCD
    srcX, srcY = sources.getX(), sources.getY()
    # Get the zeropoint
    zeropoint = (2.5 * np.log10(cal_md.get("FLUXMAG0")))
    # Get the parent ID
    parentID = sources.get('parent')
    # Check the star/galaxy separation
    extendClass = sources.get('classification.extendedness')
    # Get the nChild
    nChild = sources.get('deblend.nchild')
    # Get the aperture corrections
    # apcorr = sources.get('correctfluxes.apcorr')
    apcorr = sources.get('flux.sinc')

    # For Stars: Get these parameters
    # Get the PSF flux and its error
    flux, ferr = sources.getPsfFlux(), sources.getPsfFluxErr()
    # Convert them into magnitude and its error
    mag,  merr = 2.5*np.log10(flux), 2.5/np.log(10)*(ferr/flux)
    mag = zeropoint - mag

    apcorr = zeropoint - 2.5*np.log10(apcorr)

    # X, Y locations of the fake stars
    fakeList = collections.defaultdict(tuple)
    # Regular Expression
    # Search for keywords like FAKE12
    fakename = re.compile('FAKE([0-9]+)')
    # Go through all the keywords
    counts = 0
    for card in cal_md.names():
        # To see if the card matches the pattern
        m = fakename.match(card)
        if m is not None:
            # Get the X,Y location for fake object
            x,y    = map(float, (cal_md.get(card)).split(','))
            # Get the ID or index of the fake object
            fakeID = int(m.group(1))
            fakeList[counts] = [fakeID, x, y]
            counts += 1

    # Match the fake object to the source list
    srcIndex = collections.defaultdict(list)
    for fid, fcoord  in fakeList.items():
        separation = np.sqrt(np.abs(srcX-fcoord[1])**2 +
                             np.abs(srcY-fcoord[2])**2)
        matched = (separation <= tol)
        matchId = np.where(matched)[0]
        matchSp = separation[matchId]
        sortId = [matchId for (matchSp, matchId) in sorted(zip(matchSp, matchId))]
        # DEBUG:
        # print fid, fcoord, matchId
        # print sortId, sorted(matchSp), matchId
        # Select the index of all matched object
        srcIndex[fid] = sortId

    # Return the source list
    mapper = SchemaMapper(sources.schema)
    mapper.addMinimalSchema(sources.schema)
    newSchema = mapper.getOutputSchema()
    newSchema.addField('fakeId', type=int,
                       doc='id of fake source matched to position')
    srcList = SourceCatalog(newSchema)
    srcList.reserve(sum([len(s) for s in srcIndex.values()]))

    # Return a list of interesting parameters
    #srcParam = collections.defaultdict(list)
    srcParam = []
    nFake = 0
    for matchIndex in srcIndex.values():
        # Check if there is a match
        if len(matchIndex) > 0:
            # Only select the one with the smallest separation
            ss = matchIndex[0]
            fakeObj = fakeList[nFake]
            diffX = srcX[ss] - fakeObj[1]
            diffY = srcY[ss] - fakeObj[2]
            paramList = (fakeObj[0], fakeObj[1], fakeObj[2],
                         mag[ss], merr[ss], apcorr[ss], diffX, diffY,
                         parentID[ss], nChild[ss], extendClass[ss])
            srcParam.append(paramList)
        else:
            fakeObj = fakeList[nFake]
            paramList = (fakeObj[0], fakeObj[1], fakeObj[2],
                         0, 0, -1, -1, -1, -1, -1, -1)
            srcParam.append(paramList)
        # Go to another fake object
        nFake += 1

    # Make a numpy record array
    srcParam = np.array(srcParam, dtype=[('fakeID', int),
                                         ('fakeX', float),
                                         ('fakeY', float),
                                         ('psfMag', float),
                                         ('psfMagErr', float),
                                         ('apCorr', float),
                                         ('diffX', float),
                                         ('diffY', float),
                                         ('parentID', int),
                                         ('nChild', int),
                                         ('extendClass', float)])

    return srcIndex, srcParam, srcList, zeropoint


def main():

    #TODO: this should use the LSST/HSC conventions
    parser = argparse.ArgumentParser()
    parser.add_argument('rootDir', help='root dir of data repo')
    parser.add_argument('visit',   help='id of visit', type=int)
    parser.add_argument('ccd',     help='id of ccd',   type=int)
    parser.add_argument('tol',     help='tolerence in matching', type=float)
    args = parser.parse_args()

    # Get the information of the fake objects from the output source catalog
    (fakeIndex, fakeParam, fakeList, zp) = getStars(args.rootDir,
                                                    args.visit,
                                                    args.ccd,
                                                    args.tol)

    root = args.rootDir
    temp = root.split("/")
    if root[-1] is '/':
        rerun = temp[-2]
    else:
        rerun = temp[-1]

    outTxt = rerun + '_' + str(args.visit).strip() + '_' + \
             str(args.ccd).strip() + '_match.txt'

    fakeID = fakeParam['fakeID']
    psfMag = fakeParam['psfMag']
    psfErr = fakeParam['psfMagErr']
    apCorr = fakeParam['apCorr']
    parent = fakeParam['parentID']
    nchild = fakeParam['nChild']
    fakeX  = fakeParam['fakeX']
    fakeY  = fakeParam['fakeY']
    diffX  = fakeParam['diffX']
    diffY  = fakeParam['diffY']
    extend = fakeParam['extendClass']

    # Number of injected fake objects, and the number of the objects recovered
    # by the pipeline (using the selected tol during matching)
    nInject = len(fakeID)
    nMatch  = len(np.argwhere(psfMag))

    output = open(outTxt, 'w')

    # Print out some information
    print '###################################################################'
    print "# Number of Injected Objects : %d" % nInject
    print "# Number of Matched  Objects : %d" % nMatch
    print "# The zeropoint of this CCD is %6.3f" % zp
    print "# Visit = %d   CCD = %d" % (args.visit, args.ccd)
    print '###################################################################'
    print "# FakeX    FakeY   DiffX   DiffY   PSFMag   PSFMagErr  SincMag  " + \
          "Matched   Deblend "

    header = "# FakeX    FakeY   DiffX   DiffY   PSFMag   PSFMagErr  SincMag  " + \
             "Matched   Deblend   Extended   nChild  Visit  CCD\n"
    output.write(header)

    for i in range(nInject):
       if len(fakeIndex[i]) > 1:
           matched = "multiple"
       elif psfMag[i] > 0:
           matched = "  single"
       else:
           matched = " nomatch"

       if (parent[i] > 0):
           deblend = "deblend"
       else:
           deblend = "isolate"

       print "%7.2f   %7.2f   %6.2f   %6.2f  %7.3f  %6.3f  %7.3f  %s  %s" % (
             fakeX[i], fakeY[i], diffX[i], diffY[i], psfMag[i], psfErr[i],
             apCorr[i], matched, deblend)

       line = "%7.2f   %7.2f   %6.2f   %6.2f  %7.3f  %6.3f  %7.3f  %s  %s  %d  %d  %d  %d \n" % (
              fakeX[i], fakeY[i], diffX[i], diffY[i], psfMag[i], psfErr[i],
              apCorr[i], matched, deblend, extend[i], nchild[i],
              args.visit, args.ccd)
       output.write(line)

    output.close()

if __name__=='__main__':
    main()
