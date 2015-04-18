#!/usr/bin/env python

import numpy
import lsst.afw.geom  as afwGeom
import lsst.afw.table as afwTable

def writeFakeData():
    """Write two text files with columns: id ra dec"""
    
    n = 20
    ra  = [10.0 + numpy.random.uniform(0.0, 1.0, n)]
    dec = [20.0 + numpy.random.uniform(0.0, 1.0, n)]
    ra.append(ra[0] + numpy.random.normal(0.0, (0.2*afwGeom.arcseconds).asDegrees(), n))
    dec.append(dec[0] + numpy.random.normal(0.0, (0.2*afwGeom.arcseconds).asDegrees(), n))

    id_offsets = 1000, 2000
    for i in 0, 1:
        filename = "coord{id}.dat".format(id=i+1)
        with open(filename, 'w') as fp:
            for i_p in range(n):
                out = "{id} {ra} {dec}\n".format(id=id_offsets[i]+i_p, ra=ra[i][i_p], dec=dec[i][i_p])
                fp.write(out)

                
def loadSourceCatalog(filename):
    """Load file into a SourceCatalog.  Assuming id, ra, and dec are columns 1, 2, and 3"""
    
    schema = afwTable.SourceTable.makeMinimalSchema()
    table  = afwTable.SourceTable.make(schema)
    scat   = afwTable.SourceCatalog(table)

    with open(filename) as fp:
        for line in fp.readlines():
            id, ra, dec = line.split()
            s = scat.addNew()
            s.setId(int(id))
            s.setRa(float(ra)*afwGeom.degrees)
            s.setDec(float(dec)*afwGeom.degrees)
            
    return scat
    

    
def main():

    # write some fake data
    writeFakeData()

    # read in the fake data
    scat1 = loadSourceCatalog("coord1.dat")
    scat2 = loadSourceCatalog("coord2.dat")
    
    # do the match
    distance = 1.0*afwGeom.arcseconds
    matches = afwTable.matchRaDec(scat1, scat2, distance)    

    # print the IDs
    for match in matches:
        s1, s2, dist = match
        print s1.getId(), s2.getId(), (dist*afwGeom.radians).asArcseconds()
        
    
    
if __name__ == '__main__':
    main()
