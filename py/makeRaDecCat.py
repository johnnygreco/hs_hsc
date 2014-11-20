import numpy as np
import numpy.random.uniform as uniform

def getRaDecRangePatch(dataId):

    rangeRaDec = ''

    return rangeRaDec

def getRaDecRangeCcd(dataId):

    rangeRaDec = ''

    return rangeRaDec

def getRandomRaDec(nRand, minRa, minDec, maxRa, maxDec, rad=None):
    """
    Get random RA and Dec array
    """
    raArr  = uniform(low=minRa,  high=maxRa,  size=nRand)
    decArr = uniform(low=minDec, high=maxDec, size=nRand)

def main(nRand, dataId=None, rangeRaDec=None):

    return randomRaDec

