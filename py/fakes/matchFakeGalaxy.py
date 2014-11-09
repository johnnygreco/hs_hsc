#!/usr/bin/env python
"""
matchFakes.py
matches fakes based on position stored in the calibrated exposure image header
"""

import lsst.daf.persistence as dafPersist
from lsst.afw.table import SourceCatalog, SchemaMapper
import lsst.afw.geom.ellipses as geomEllip
import numpy as np
import argparse
import re
import collections
import pyfits as fits

def getSizeAndShape(m):
    """
    Get the major axis radius, axis ratio, position angle from moments
    m = Moments
    """
    ixx, ixy, iyy = m.getIxx(), m.getIxy(), m.getIyy()
    # convert to an ellipse (note that theta is in radians and is not an
    # Angle object)
    ellipse = geomEllip.Axes(m)
    a, b, theta = ellipse.getA(), ellipse.getB(), ellipse.getTheta()

    pa = theta * 180.0 / np.pi

    return a, (b/a), pa


def getGalaxy(rootdir, visit, ccd, tol):
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
    nSources = srcX.shape[0]
    # Get the zeropoint
    zeropoint = (2.5 * np.log10(cal_md.get("FLUXMAG0")))
    # Get the parent ID
    parentID = sources.get('parent')
    # Check the star/galaxy separation
    extendClass = sources.get('classification.extendedness')

    # For Galaxies: Get these parameters
    # 1. Get the Kron flux and its error
    fluxKron, ferrKron = sources.get('flux.kron'), sources.get('flux.kron.err')
    magKron, merrKron = (zeropoint - 2.5*np.log10(fluxKron)), (2.5/np.log(10)*
                                                            (ferrKron/fluxKron))
    # 2. Get the CModel flux and its error
    fluxCmod, ferrCmod = sources.get('cmodel.flux'), sources.get('cmodel.flux.err')
    magCmod, merrCmod = (zeropoint - 2.5*np.log10(fluxCmod)), (2.5/np.log(10)*
                                                            (ferrCmod/fluxCmod))
    # 3. Get the Exponential flux and its error
    fluxExp, ferrExp = sources.get('cmodel.exp.flux'), sources.get('cmodel.exp.flux.err')
    magExp, merrExp = (zeropoint - 2.5*np.log10(fluxExp)), (2.5/np.log(10)*
                                                            (ferrExp/fluxExp))
    # 4. Get the de Vacouleurs flux and its error
    fluxDev, ferrDev = sources.get('cmodel.dev.flux'), sources.get('cmodel.dev.flux.err')
    magDev, merrDev = (zeropoint - 2.5*np.log10(fluxDev)), (2.5/np.log(10)*
                                                            (ferrDev/fluxDev))
    # 5. Get the SDSS shapes (Re, b/a, PA)
    sdssR  = []
    sdssBa = []
    sdssPa = []
    expR  = []
    expBa = []
    expPa = []
    devR  = []
    devBa = []
    devPa = []
    for ss in sources:
        # 5. Get the SDSS shapes (Re, b/a, PA)
        sdssMoment = ss.get('shape.sdss')
        r, ba, pa = getSizeAndShape(sdssMoment)
        sdssR.append(r)
        sdssBa.append(ba)
        sdssPa.append(pa)
        # 6. Get the Exponential shapes (Re, b/a, PA)
        expMoment = ss.get('cmodel.exp.ellipse')
        r, ba, pa = getSizeAndShape(expMoment)
        expR.append(r)
        expBa.append(ba)
        expPa.append(pa)
        # 7. Get the de Vaucouleurs shapes (Re, b/a, PA)
        devMoment = ss.get('cmodel.dev.ellipse')
        r, ba, pa = getSizeAndShape(devMoment)
        devR.append(r)
        devBa.append(ba)
        devPa.append(pa)
    # 8. Get the fracDev
    fracDev = sources.get('cmodel.fracDev')

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
            # TODO: actually get the one with minimum separation
            ss = matchIndex[0]
            fakeObj = fakeList[nFake]
            diffX = srcX[ss] - fakeObj[1]
            diffY = srcY[ss] - fakeObj[2]
            paramList = (fakeObj[0], fakeObj[1], fakeObj[2],
                         magKron[ss], merrKron[ss], magCmod[ss], merrCmod[ss],
                         magExp[ss], merrExp[ss], magDev[ss], merrDev[ss],
                         sdssR[ss], sdssBa[ss], sdssPa[ss],
                         expR[ss], expBa[ss], expPa[ss],
                         devR[ss], devBa[ss], devPa[ss],
                         diffX, diffY, fracDev[ss],
                         parentID[ss], extendClass[ss])
            srcParam.append(paramList)
        else:
            paramList = (fakeObj[0], fakeObj[1], fakeObj[2],
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         -1, -1, 0, -1, -1)
            srcParam.append(paramList)
        # Go to another fake object
        nFake += 1

    # Make a numpy record array
    srcParam = np.array(srcParam, dtype=[('fakeID', int),
                                         ('fakeX', float),
                                         ('fakeY', float),
                                         ('magKron', float),
                                         ('errKron', float),
                                         ('magCmod', float),
                                         ('errCmod', float),
                                         ('magExp', float),
                                         ('errExp', float),
                                         ('magDev', float),
                                         ('errDev', float),
                                         ('sdssR', float),
                                         ('sdssBa', float),
                                         ('sdssPa', float),
                                         ('expR', float),
                                         ('expBa', float),
                                         ('expPa', float),
                                         ('devR', float),
                                         ('devBa', float),
                                         ('devPa', float),
                                         ('diffX', float),
                                         ('diffY', float),
                                         ('fracDev', float),
                                         ('parentID', int),
                                         ('extendClass', float)])

    return srcIndex, srcParam, srcList, zeropoint

def save_to_fits(params, root, visit, ccd):

    temp = root.split("/")
    if root[-1] is '/':
        rerun = temp[-2]
    else:
        rerun = temp[-1]

    outFits = rerun + '_' + str(visit).strip() + '_' + str(ccd).strip() + '_' +\
              '_match.fits'

    tabHdu = fits.BinTableHDU.from_columns([
        fits.Column(name='fakeID', format='A', array=params['fakeID']),
        fits.Column(name='fakeX',  format='E', array=params['fakeX']),
        fits.Column(name='fakeY',  format='E', array=params['fakeY']),
        fits.Column(name='diffX',  format='E', array=params['diffX']),
        fits.Column(name='diffY',  format='E', array=params['diffY']),
        fits.Column(name='magKron',  format='E', array=params['magKron']),
        fits.Column(name='errKron',  format='E', array=params['errKron']),
        fits.Column(name='magCmod',  format='E', array=params['magCmod']),
        fits.Column(name='errCmod',  format='E', array=params['errCmod']),
        fits.Column(name='magExp',  format='E', array=params['magExp']),
        fits.Column(name='errExp',  format='E', array=params['errExp']),
        fits.Column(name='magDev',  format='E', array=params['magDev']),
        fits.Column(name='errDev',  format='E', array=params['errDev']),
        fits.Column(name='sdssR',  format='E', array=params['sdssR']),
        fits.Column(name='sdssBa',  format='E', array=params['sdssBa']),
        fits.Column(name='sdssPa',  format='E', array=params['sdssPa']),
        fits.Column(name='expR',  format='E', array=params['expR']),
        fits.Column(name='expBa',  format='E', array=params['expBa']),
        fits.Column(name='expPa',  format='E', array=params['expPa']),
        fits.Column(name='devR',  format='E', array=params['devR']),
        fits.Column(name='devBa',  format='E', array=params['devBa']),
        fits.Column(name='devPa',  format='E', array=params['devPa']),
        fits.Column(name='fracDev',  format='E', array=params['fracDev']),
        fits.Column(name='parentID',  format='I', array=params['parentID']),
        fits.Column(name='extendClass', format='E', array=params['extendClass'])
    ])

    tabHdu.writeto(outFits)

    return outFits

def main():

    #TODO: this should use the LSST/HSC conventions
    parser = argparse.ArgumentParser()
    parser.add_argument('rootDir', help='root dir of data repo')
    parser.add_argument('visit',   help='id of visit', type=int)
    parser.add_argument('ccd',     help='id of ccd',   type=int)
    parser.add_argument('tol',     help='tolerence in matching', type=float)
    args = parser.parse_args()

    # Get the information of the fake objects from the output source catalog
    (fakeIndex, fakeParam, fakeList, zp) = getGalaxy(args.rootDir,
                                                     args.visit,
                                                     args.ccd,
                                                     args.tol)

    #outFits = save_to_fits(fakeParam, args.rootDir, args.visit, args.ccd)

    fakeID = fakeParam['fakeID']
    magKron = fakeParam['magKron']
    magCmod = fakeParam['magCmod']
    errCmod = fakeParam['errCmod']
    magExp = fakeParam['magExp']
    magDev = fakeParam['magDev']
    parent = fakeParam['parentID']
    fakeX  = fakeParam['fakeX']
    fakeY  = fakeParam['fakeY']
    diffX  = fakeParam['diffX']
    diffY  = fakeParam['diffY']

    # Number of injected fake objects, and the number of the objects recovered
    # by the pipeline (using the selected tol during matching)
    nInject = len(fakeID)
    nMatch  = len(np.argwhere(magCmod))

    # Print out some information
    print '###################################################################'
    print "# Number of Injected Objects : %d" % nInject
    print "# Number of Matched  Objects : %d" % nMatch
    print "# The zeropoint of this CCD is %6.3f" % zp
    print "# Visit = %d   CCD = %d" % (args.visit, args.ccd)
    print '###################################################################'
    print "# FakeX    FakeY    DiffX  DiffY  CmodMag  CmodErr  " + \
          "ExpMag  DevMag  KronMag  Matched  Deblend "

    for i in range(nInject):
       if len(fakeIndex[i]) > 1:
           matched = "multiple"
       elif magCmod[i] > 0:
           matched = "  single"
       else:
           matched = " nomatch"

       if (parent[i] > 0):
           deblend = "deblend"
       else:
           deblend = "isolate"

       print "%7.2f   %7.2f   %6.2f   %6.2f  %7.3f  %6.3f  %7.3f  %7.3f  %7.3f  %s  %s" % (
             fakeX[i], fakeY[i], diffX[i], diffY[i], magCmod[i], errCmod[i],
             magExp[i], magDev[i], magKron[0], matched, deblend)

if __name__=='__main__':
    main()
