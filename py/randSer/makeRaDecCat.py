#!/usr/bin/env python
"""
Create N random selected (RA, Dec) pairs within a desired region.
The inputs can be:
   1. RA, Dec range: (minRa, maxRa, minDec, maxDec)
   2. DataId for single frame image: (visit, ccd)
   3. DataId for coadd image: (tract, patch, filter)
"""
import numpy as np
from numpy.random import uniform

def getImageRaDecRange(rootDir, dataId, dataType='calexp'):
    """
    Get the Ra,Dec range for certain single frame or coadded image using the WCS
    """

    import lsst.daf.persistence as dafPersist
    import lsst.afw.image as afwImage

    butler = dafPersist.Butler(rootDir)
    exposure = butler.get(dataType, dataId)

    bboxI = exposure.getBBox(afwImage.PARENT)
    wcs   = exposure.getWcs()

    minPix = wcs.pixelToSky(bboxI.getMinX(), bboxI.getMinY())
    maxPix = wcs.pixelToSky(bboxI.getMaxX(), bboxI.getMaxY())

    ra1  = minPix.getLongitude().asDegrees()
    ra2  = maxPix.getLongitude().asDegrees()
    dec1 = minPix.getLatitude().asDegrees()
    dec2 = maxPix.getLatitude().asDegrees()

    minRa, maxRa   = min(ra1, ra2), max(ra1, ra2)
    minDec, maxDec = min(dec1, dec2), max(dec1, dec2)

    return [minRa, maxRa, minDec, maxDec]


def getRandomRaDec(nRand, minRa, maxRa, minDec, maxDec, rad=None):
    """
    Randomly select Ra,Dec pairs from the input Ra,Dec range
    """
    if minRa > maxRa or minDec > maxDec:
        raise Exception('Please provide appropriate Ra,Dec range !')

    if rad is None:
        raArr  = uniform(low=minRa,  high=maxRa,  size=nRand)
        decArr = uniform(low=minDec, high=maxDec, size=nRand)
        return zip(raArr, decArr)
    else:
        import lsst.afw.coord as afwCoord
        import lsst.afw.geom  as afwGeom

        minSep = float(rad)
        raArr  = []
        decArr = []
        numTry = 0
        while len(raArr) < nRand:
            if numTry == 0:
                raArr.append(uniform(low=minRa,   high=maxRa))
                decArr.append(uniform(low=minDec, high=maxDec))
                numTry += 1
            else:
                raTry  = uniform(low=minRa,  high=maxRa)
                decTry = uniform(low=minDec, high=maxDec)
                coordTry = afwCoord.Coord(afwGeom.Point2D(raTry, decTry))
                nExist = len(raArr)
                sepGood = True
                for ii in range(nExist):
                    coordTest = afwCoord.Coord(afwGeom.Point2D(raArr[ii],
                                                               decArr[ii]))
                    sep = coordTry.angularSeparation(coordTest).asArcseconds()
                    if sep <= minSep:
                        print "## BAD ONE %8d : %8d -- %f10 <= %f10 !!" % (len(raArr),
                                                                     ii, sep, minSep)
                        sepGood = False
                        break
                if sepGood:
                    raArr.append(raTry)
                    decArr.append(decTry)
        return zip(raArr, decArr)


def plotRandomRaDec(randomRaDec, rangeRaDec=None):
    """
    Plot the distribution of radom Ra, Dec for examination
    """
    import matplotlib.pyplot as plt

    plt.scatter(*zip(*randomRaDec))
    plt.xlabel(r'RA (J2000)',  fontsize=20, labelpad=20)
    plt.ylabel(r'DEC (J2000)', fontsize=20, labelpad=20)

    if rangeRaDec is not None:
        if type(rangeRaDec) is dict:
            raMin,  raMax  = rangeRaDec['raMin'], rangeRaDec['raMax']
            decMin, decMax = rangeRaDec['raMin'], rangeRaDec['raMax']
        else:
            raMin,  raMax  = rangeRaDec[0], rangeRaDec[1]
            decMin, decMax = rangeRaDec[0], rangeRaDec[1]
        raDec0   = (raMin, decMin)
        raRange  = (raMax  - raMin)
        decRange = (decMax - decMin)

    plt.gcf().savefig('randomRaDec.png')

    return None

def makeRaDecCat(nRand, dataId=None, rangeRaDec=None, rad=None,
                 rootDir='/lustre/Subaru/SSP/rerun/song/cosmos-i2',
                 inputCat=None, plot=False):
    """
    Generate nRand random RA,Dec pairs in a desired region of sky
    The region can be defined by:
       1) Single frame image
       2) Coadd image
       3) RA, Dec range
    """

    if dataId is not None:
        if 'visit' in dataId.keys() and 'ccd' in dataId.keys():
            # Input should be a single frame image
            rangeRaDec  = getImageRaDecRange(rootDir, dataId)
            randomRaDec = getRandomRaDec(nRand, rangeRaDec[0], rangeRaDec[1],
                                         rangeRaDec[2], rangeRaDec[3], rad=rad)
        elif 'tract' in dataId.keys() and 'patch' in dataId.keys() and 'filter' in dataId.keys:
            # Input should be a coadd image
            rangeRaDec  = getImageRaDecRange(rootDir, dataId,
                                            dataType='deepCoadd_calexp')
            randomRaDec = getRandomRaDec(nRand, rangeRaDec[0], rangeRaDec[1],
                                         rangeRaDec[2], rangeRaDec[3], rad=rad)
        else:
            raise KeyError('Please provide the correct dataId !')

    elif rangeRaDec is not None:

        if type(rangeRaDec) is dict:
            rKeys = rangeRaDec.keys()
            if 'raMin' in rKeys and 'raMax' in rKeys and 'decMin' in rKeys and 'decMax' in rKeys:
                randomRaDec = getRandomRaDec(nRand, rangeRaDec['minRa'],
                                            rangeRaDec['maxRa'],
                                            rangeRaDec['minDec'],
                                            rangeRaDec['maxDec'], rad=rad)
            else:
                raise KeyError('Please provide the correct rangeRaDec!')
        elif type(rangeRaDec) is list or type(rangeRaDec).__module__ == 'numpy':
            if len(rangeRaDec) >= 4:
                randomRaDec = getRandomRaDec(nRand, rangeRaDec[0],
                                             rangeRaDec[1], rangeRaDec[2],
                                             rangeRaDec[3], rad=rad)
            else:
                raise Exception('randomRaDec should have at least 4 elements!')
        else:
            raise Exception('randomRaDec need to be Dict/List/Numpy.Array')

    else:
        raise Exception("Need to provide either dataId or rangeRaDec")

    if plot:
        plotRandomRaDec(randomRaDec, rangeRaDec=rangeRaDec)

    if inputCat is not None:

        import os
        import astropy.table

        if os.path.exists(inputCat):

            galCat = astropy.table.Table.read(inputCat, format='fits')
            outCat = inputCat.strip().replace('.fits', '_radec.fits')

            nGal   = len(galCat)

            if nGal == nRand:
                raArr, decArr = np.array(zip(*randomRaDec))
                raCol  = astropy.table.Column(name='RA',  data=raArr)
                decCol = astropy.table.Column(name='Dec', data=decArr)
                galCat.add_columns([raCol, decCol])
            elif nGal < nRand:
                import random
                raArr, decArr = np.array(zip(*random.sample(randomRaDec, nGal)))
                raCol  = astropy.table.Column(name='RA',  data=raArr)
                decCol = astropy.table.Column(name='Dec', data=decArr)
                galCat.add_columns([raCol, decCol])
            else:
                raise Exception('There are not enough random RA, Dec!')

            galCat.write(outCat, format='fits', overwrite=True)

        else:
            raise Exception('Can not find input catalog %s!' % inputCat)

    return None #randomRaDec

