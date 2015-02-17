#!/usr/bin/env python
"""
matchFakes.py
matches fakes based on position stored in the calibrated exposure image header
"""

import lsst.daf.persistence as dafPersist
from lsst.afw.table import SourceCatalog, SchemaMapper
import lsst.afw.geom
import lsst.pex.exceptions
import numpy as np
import argparse
import re
import collections
import astropy.table
import lsst.afw.geom.ellipses

def getMag(flux, fluxerr, zeropoint):
    """
    return the magnitude and error
    """
    mag, magerr = -2.5 * np.log10(flux), 2.5/np.log(10.0)*fluxerr/flux
    return (mag.T + zeropoint).T, magerr


def getEllipse(quad):
    """
    returns the semi-major axis, axes ratio and PA for a given quadrupole moment
    """
    e = lsst.afw.geom.ellipses.Axes(quad)
    return e.getA(), e.getB()/e.getA(), e.getTheta() * 180.0/np.pi


def matchToFakeCatalog(sources, fakeCatalog):
    """
    match to the fake catalog and append those columns to the source table
    
    this assumes the sources are an astropy table
    or it will throw a TypeError
    """
    if not isinstance(sources, astropy.table.Table):
        raise TypeError("expects and astropy table for sources, use getAstroTable to convert")
    
    
    fakes = astropy.table.Table().read(fakeCatalog)
    fakes.rename_column('ID', 'fakeId')
    return astropy.table.join(sources, fakes, keys='fakeId', join_type='left')


def getFakeMatchesHeader(cal_md, sources, tol=1.0):
    """
    returns the fake matches based on the information in the header
    
    returns a tuple with:
    the positions in pixels of the fake sources added to the chip
    the match is in a dictionary of the form: {fakeid:[ind_of_match_in_sources,...],...}

    look within a tolerance of 1 pixel in each direction
    """
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
        distX = srcX - fcoord[0]
        distY = srcY - fcoord[1]
        matched = (np.abs(distX) < tol) & (np.abs(distY) < tol)
        srcIndex[fid] = np.where(matched)[0]

    return fakeXY, srcIndex


def getFakeMatchesRaDec( sources, radecCatFile, bbox, wcs, tol=1.0):
    """
    Returns the fake matches based on an radec match to the catalog of fake sources inputed.
    
    Args:
    sources: source table object
    radecCatFile: filename for fits file of fake object table, including ra/dec
    bbox: Bounding Box of exposure (ccd or patch) within which to match fakes
    wcs: Wcs of source image
    
    KeywordArgs:
    tol: tolerance within which to match sources, given in PIXELS

    Returns:
    fakeXY: set of pixel positions of fake sources
    srcIdx: matches to fake sources in a dictionary of the form {fakeid:[ind_of_match_in_sources,...],...}

    Raise:
    IOError: couldn't open radecCatFile
    """
    
    fakeXY = collections.defaultdict(tuple)
    try:
        fakeCat = astropy.table.Table().read(radecCatFile)
    except IOError:
        raise
    
    for fakeSrc in fakeCat:
        fakeCoord = wcs.skyToPixel(lsst.afw.geom.Angle(fakeSrc['RA'], lsst.afw.geom.degrees),
                                   lsst.afw.geom.Angle(fakeSrc['Dec'], lsst.afw.geom.degrees))
        if bbox.contains(fakeCoord):
            fakeXY[int(fakeSrc['ID'])] = (fakeCoord.getX(), fakeCoord.getY())
        
    srcX, srcY = sources.getX(), sources.getY()
    srcIndex = collections.defaultdict(list)
    for fid, fcoord in fakeXY.items():
        distX = srcX - fcoord[0]
        distY = srcY - fcoord[1]
        matched = (np.abs(distX) < tol) & (np.abs(distY) < tol)
        srcIndex[fid] = np.where(matched)[0]
    
    return fakeXY, srcIndex



def getFakeSources(butler, dataId, tol=1.0, extraCols=('zeropoint', 'visit', 'ccd'),
                   includeMissing=False, footprints=False, radecMatch=None):
    """Get list of sources which agree in pixel position with fake ones with tol
    
    this returns a sourceCatalog of all the matched fake objects,
    note, there will be duplicates in this list, since I haven't checked deblend.nchild,
    and I'm only doing a tolerance match, which could include extra sources
    
    the outputs can include extraCols as long as they are one of:
      zeropoint, visit, ccd, thetaNorth, pixelScale

    if includeMissing is true, then the pipeline looks at the fake sources
    added in the header and includes an entry in the table for sources without
    any measurements, specifically the 'id' column will be 0

    radecMatch is the fakes table. if it's not None(default), then do an ra/dec 
    match with the input catalog instead of looking in the header for where the 
    sources where added
    """
    
    availExtras = {'zeropoint':{'type':float, 'doc':'zeropoint'}, 
                   'visit':{'type':int, 'doc':'visit id'}, 
                   'ccd':{'type':int, 'doc':'ccd id'},
                   'thetaNorth':{'type':lsst.afw.geom.Angle, 'doc':'angle to north'},
                   'pixelScale':{'type':float, 'doc':'pixelscale in arcsec/pixel'}}
    
    if not np.in1d(extraCols, availExtras.keys()).all():
        print "extraCols must be in ",availExtras

    if not 'filter' in dataId:
        sources = butler.get('src', dataId, 
                             flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS)
        cal = butler.get('calexp', dataId)
        cal_md = butler.get('calexp_md', dataId)
    else:
        sources = butler.get('deepCoadd_src', dataId, flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS)
        cal = butler.get('deepCoadd', dataId)
        cal_md = butler.get('deepCoadd_md', dataId)

    if ('pixelScale' in extraCols) or ('thetaNorth' in extraCols):
        wcs = cal.getWcs()
        availExtras['pixelScale']['value'] =  wcs.pixelScale().asArcseconds()
        availExtras['thetaNorth']['value'] = lsst.afw.geom.Angle(
            np.arctan2(*tuple(wcs.getLinearTransform().invert()
                              (lsst.afw.geom.Point2D(1.0,0.0)))))
    if 'visit' in extraCols:
        availExtras['visit']['value'] = dataId['visit']
    if 'ccd' in extraCols:
        availExtras['ccd']['value'] = dataId['ccd']
    if 'zeropoint' in extraCols:
        availExtras['zeropoint']['value'] = 2.5*np.log10(cal_md.get('FLUXMAG0'))

        
    if radecMatch is None:
        fakeXY, srcIndex = getFakeMatchesHeader(cal_md, sources, tol=tol)
    else:
        fakeXY, srcIndex = getFakeMatchesRaDec(sources, radecMatch, 
                                               lsst.afw.geom.Box2D(cal.getBBox(lsst.afw.image.PARENT)),
                                               cal.getWcs(), 
                                               tol=tol)

    mapper = SchemaMapper(sources.schema)
    mapper.addMinimalSchema(sources.schema)
    newSchema = mapper.getOutputSchema()
    newSchema.addField('fakeId', type=int, doc='id of fake source matched to position')
    newSchema.addField('fakeOffset', type=lsst.afw.geom.Point2D,
                       doc='offset from input fake position (pixels)')
    for extraName in set(extraCols).intersection(availExtras):
        newSchema.addField(extraName, type=availExtras[extraName]['type'],
                           doc=availExtras[extraName]['doc'])

    srcList = SourceCatalog(newSchema)
    srcList.reserve(sum([len(s) for s in srcIndex.values()]) + 
                    (0 if not includeMissing else srcIndex.values().count([])))

    for ident, sindlist in srcIndex.items():
        if includeMissing and len(sindlist)==0:
            newRec = srcList.addNew()
            newRec.set('fakeId', ident)
            newRec.set('id', 0)
        for ss in sindlist:
            newRec = srcList.addNew()
            newRec.assign(sources[ss], mapper)
            newRec.set('fakeId', ident)
            newRec.set('fakeOffset', 
                       lsst.afw.geom.Point2D(sources[ss].get('centroid.sdss').getX() - fakeXY[ident][0],
                                             sources[ss].get('centroid.sdss').getY() - fakeXY[ident][1]))

    if includeMissing:
        srcList = srcList.copy(deep=True)

    for extraName in set(extraCols).intersection(availExtras):
        tempCol = srcList.get(extraName)
        tempCol.fill(availExtras[extraName]['value'])

    return srcList


def getAstroTable(src, mags=True):
    """
    returns an astropy table with all the src entries
    if the entries are complex objects, it breaks them down:
      ellipse entries are broken into 
           ellipse_a = semi-major axis
           ellipse_q = axis ratio (always < 1)
           ellipse_theta = rotation of semi-major axis from chip x-axis in degrees 
    if mags is True, returns the magnitudes for all the flux columns
    """
    
    tab = astropy.table.Table()
    for name in src.schema.getNames():
        try: 
            tab.add_column(astropy.table.Column(name=name,
                                                data=src.get(name)))
        except lsst.pex.exceptions.LsstException:
            if type(src[0].get(name)) is lsst.afw.geom.ellipses.ellipsesLib.Quadrupole:
                reff, q, theta = zip(*[getEllipse(s.get(name)) for s in src])
                tab.add_column(astropy.table.Column(name=name+'_a', data=reff))
                tab.add_column(astropy.table.Column(name=name+'_q', data=q))
                tab.add_column(astropy.table.Column(name=name+'_theta', data=theta))
            elif type(src[0].get(name)) is lsst.afw.coord.coordLib.IcrsCoord:
                x, y= zip(*[(s.get(name).getRa().asDegrees(), 
                             s.get(name).getDec().asDegrees()) for s in src])
                tab.add_column(astropy.table.Column(name=name+'_ra', data=x))
                tab.add_column(astropy.table.Column(name=name+'_dec', data=y))
            else:
                tab.add_column(astropy.table.Column(name=name, 
                                                    data=np.array([s.get(name) for s in src])))
            #report angles in degrees
        if isinstance(src[0].get(name), lsst.afw.geom.Angle):
            tab.remove_column(name)
            tab.add_column(astropy.table.Column(data=[s.get(name).asDegrees()
                                                      for s in src],
                                                dtype=float, name=name))

    if mags:
        #this is a horrible hack, but I don't think we can use the slots, since 
        #not all the fluxes end up in the slots
        for col in tab.colnames:
            if (re.match('^flux\.[a-z]+$', col) or 
                re.match('^flux\.[a-z]+.apcorr$', col) or
                re.match('^cmodel.+flux$', col) or 
                re.match('^cmodel.+flux.apcorr$', col)):
                mag, magerr = getMag(tab[col], tab[col+'.err'], 
                                     tab['zeropoint'] if not re.search('apcorr', col) else 0.0)
            
                tab.add_column(astropy.table.Column(name=re.sub('flux', 'mag', col),
                                                    data=mag))
                tab.add_column(astropy.table.Column(name=re.sub('flux', 'mag', col+'.err'),
                                                    data=magerr))
                
    return tab



def returnMatchTable(rootDir, visit, ccdList, outfile=None, fakeCat=None,
                     overwrite=False, filt=None):
    """
    driver (main function) for return match to fakes
    INPUT: rootDir = rerun directory
           visit = visit id (int) (or tracts)
           ccdList = list of ccds to look at (or patches)
           outdir = output directory for matched file, None means no output written
           fakeCat = fake catalog to match to, None means the fake sources are just
                     extracted from the header of the CCDs based on position but no matching is done
           overwrite = whether to overwrite the existing output file, default is False
    OUTPUT: returns an astropy.table.Table with all the entries from the source catalog for 
            objects which match in pixel position to the fake sources
    """
    
    butler = dafPersist.Butler(rootDir)
    slist = None

    for ccd in ccdList:
        if filt is None:
            print 'doing ccd %d'%int(ccd)
            temp = getFakeSources(butler,
                                  {'visit':visit, 'ccd':int(ccd)}, 
                                  includeMissing=True,
                                  extraCols=('visit', 'ccd', 
                                             'zeropoint', 'pixelScale', 
                                             'thetaNorth'), radecMatch=fakeCat)
        else:
            print 'doing patch %s'%ccd
            temp = getFakeSources(butler,
                                  {'tract':visit, 'patch':ccd,
                                   'filter':filt},
                                  includeMissing=True, extraCols=('thetaNorth', 'pixelScale',
                                                                  'zeropoint'),
                                  radecMatch=fakeCat)
        if slist is None:
            slist = temp.copy(True)
        else:
            slist.extend(temp, True)
        del temp

    astroTable = getAstroTable(slist, mags=True)
    
    if fakeCat is not None:
        astroTable = matchToFakeCatalog(astroTable, fakeCat)

    if outfile is not None:        
        astroTable.write(outfile+'.fits', 
                         format='fits',
                         overwrite=overwrite)

    return astroTable

    
if __name__=='__main__':
        #TODO: this should use the LSST/HSC conventions
    parser = argparse.ArgumentParser()
    parser.add_argument('rootDir', help='root dir of data repo')
    parser.add_argument('visit', 
                        help='id of visit (or tract, if filter is specified)', type=int)
    parser.add_argument('-f', '--filter', dest='filt',
                        help='name of filter, if none assume single visit',
                        default=None)
    parser.add_argument('--ccd', nargs='+', help='id of ccd(s) or patches')
    parser.add_argument('-o', help='outputfilename', default=None, dest='outfile')
    parser.add_argument('-c', help='fake catalog', default=None, dest='fakeCat')
    parser.add_argument('-w', '--overwrite', help='over write output file', 
                        dest='ow', default=False, action='store_true')
    args = parser.parse_args()

    
    returnMatchTable(args.rootDir, args.visit, args.ccd, args.outfile, args.fakeCat,
                     overwrite=args.ow, filt=args.filt)
