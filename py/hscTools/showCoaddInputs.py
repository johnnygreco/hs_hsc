import math
import lsst.afw.geom as afwGeom

def showCoaddInputs(butler, dataId, pos=None, coaddType="deepCoadd"):
    """Show the inputs for the specified dataId, optionally at the specified position
    @param butler    Butler to provide inputs
    @param dataId    Desired inputs
    @param pos       Position we want to know about (PointD with pixel coords in coadd or afwCoord.Coord)
    @param coaddType Type of coadd to examine
    """
    coadd = butler.get(coaddType, **dataId)
    visitInputs = coadd.getInfo().getCoaddInputs().visits
    ccdInputs = coadd.getInfo().getCoaddInputs().ccds
    #
    # If we have a position to evaluate things at we'll need to know it on the input images,
    # so we need it in (ra, dec)
    #
    if pos:
        if hasattr(pos, "toFk5"):
            posSky = pos
        else:
            posSky = coadd.getWcs().pixelToSky(pos)
    else:
        posSky = None

    print "%6s %3s %7s %5s %5s" % ("visit", "ccd", "exptime", "FWHM", "weight")

    totalExpTime = 0.0
    expTimeVisits = set()
    for i in range(len(ccdInputs) + 1):
        try:
            input = ccdInputs[i]
        except:
            v = "coadd"
            calexp = coadd
        else:
            ccd = input.get("ccd")
            v = input.get("visit")
            bbox = input.getBBox()
            # It's quicker to not read all the pixels, so just read 1
            calexp = butler.get("calexp_sub", bbox=afwGeom.BoxI(afwGeom.PointI(0, 0), afwGeom.ExtentI(1, 1)),
                                visit=int(v), ccd=ccd)

        calib = calexp.getCalib()
        psf = calexp.getPsf()

        if posSky:
            pos = calexp.getWcs().skyToPixel(posSky)

            if not (v == "coadd" or bbox.contains(afwGeom.PointI(pos))):
                continue
        else:
            pos = psf.getAveragePosition()

        sigma = psf.computeShape(pos).getDeterminantRadius()

        if v == "coadd":
            exptime = totalExpTime
            ccd = ""
            weight = ""
        else:
            exptime = calib.getExptime()
            weight = "%5.2f" % (input.get("weight"))

            if v not in expTimeVisits:
                totalExpTime += exptime
                expTimeVisits.add(v)

        print  "%6s %3s %7.0f %5.2f %5s" % (v, ccd, exptime, sigma, weight)
