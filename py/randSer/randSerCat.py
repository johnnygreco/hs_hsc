import numpy as np
from numpy.random import uniform
from astropy.io import fits

def serGen(nModels, minRe=2.0, maxRe=10.0, minMag=20.0, maxMag=23.5,
           minSer=0.5, maxSer=3.5, minBa=0.4, maxBa=0.9,
           minPa=0.0, maxPa=45.0, nSer=None, outFits=None):
    """
    Generate nModels sing Serisc galaxy models
    The parameters are drawn from a flat distribution
    TODO: Drawn from other distribution?
    """
    indArr = (np.arange(nModels) + 1)
    magArr = uniform(low=minMag, high=maxMag, size=nModels)
    reArr  = uniform(low=minRe,  high=maxRe,  size=nModels)
    baArr  = uniform(low=minBa,  high=maxBa,  size=nModels)
    paArr  = uniform(low=minPa,  high=maxPa,  size=nModels)
    if nSer is None:
        serArr = uniform(low=minSer, high=maxSer, size=nModels)
    else:
        serArr = np.empty(nModels)
        serArr.fill(nSer)

    modelArr = []
    for i in range(nModels):
        modelRec = (i, magArr[i], reArr[i], serArr[i], baArr[i], paArr[i])
        modelArr.append(modelRec)

    serModels = np.array(modelArr, dtype=[('ID','int'), ('mag','float'),
                                          ('sersic_n','float'),
                                          ('reff_pix','float'),
                                          ('b_a','float'), ('theta','float')])

    if outFits is not None:
        col1 = fits.Column(name='ID',  format='I', array=indArr)
        col2 = fits.Column(name='mag', format='D', array=magArr)
        col3 = fits.Column(name='sersic_n', format='D', array=serArr)
        col4 = fits.Column(name='reff_pix', format='D', array=reArr)
        col5 = fits.Column(name='b_a',   format='D', array=baArr)
        col6 = fits.Column(name='theta', format='D', array=paArr)

        cols   = fits.ColDefs([col1, col2, col3, col4, col5, col6])
        tbhdu  = fits.BinTableHDU.from_columns(cols)
        prihdr = fits.Header()
        prihdr['NMODELS'] = nModels
        prihdu = fits.PrimaryHDU(header=prihdr)

        thdulist = fits.HDUList([prihdu, tbhdu])
        thdulist.writeto(outFits)

    return serModels

