#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import argparse
import numpy  as np
from collections import defaultdict, Iterable

import lsst.daf.persistence as dafPer
import lsst.afw.coord       as afwCoord
import lsst.pex.exceptions


def getVisits(tpDict, butler, filter='HSC-I'):
    """
    get visits that cover the tract/patch list

    From Claire Lackner's getTractPatches.py in gxystack
    """
    visitDict = defaultdict(set)

    for tract in tpDict:
        for patch in tpDict[tract]:
            dataId =  {'tract':tract,
                       'patch':"{0!s},{1!s}".format(*patch),
                       'filter':filter}
            coadd = butler.get("deepCoadd", dataId)
            try:
                ccds = coadd.getInfo().getCoaddInputs().ccds
            except lsst.pex.exceptions.LsstCppException:
                continue
            print dataId
            for ccd in ccds:
                visitDict[ccd.get('visit')].add(ccd.get('ccd'))

    visCcds = []
    for visit, ccds in visitDict.items():
        for ccd in ccds:
            visitCcds.append((visit, ccd))

    # Return a list of (visit, ccd) tuples
    return visCcds


def raDec2VisitCcd(ra, dec, root, filter='HSC-I', butler=None):
    """TODO: Docstring for raDec2VisitCcd.

    Find all the (visit, ccd) of HSC single images that cover given (RA, DEC)

    Based on Claire Lackner's getTractPatches.py in gxystack

    """

    # Get the butler if not provided
    if butler is None:
        butler = dafPer.Butler(root)
    skyMap = butler.get("deepCoadd_skyMap", immediate=True)

    # Fake a very simple source table
    schema = afwTable.SourceTable.makeMinimalSchema()
    table  = afwTable.SourceTable.make(schema)
    scat   = afwTable.SourceCatalog(table)

    # See if the input is number or array
    if isinstance(ra, (int, float)) and isinstance(dec, (int, float)):
        s = scat.addNew()
        s.setId(1)
        s.setRa(float(ra).afwGeom.degrees)
        s.setDec(float(dec).afwGeom.degrees)
    elif isinstance(ra, Iterable) and isinstance(dec, Iterable) and (len(ra) == len(dec)):
        raArr  = np.asarray(ra)
        decArr = np.asarray(dec)
        for i, (rr, dd) in enumerate(zip(raArr, decArr)):
            s = scat.addNew()
            s.setId(int(i))
            s.setRa(float(rr)*afwGeom.degrees)
            s.setDec(float(dd)*afwGeom.degrees)

    # Get the list of tract and patch
    tpList = [skyMap.findTractPatchList([scat[i].get("coord"),])
              for i in range(len(scat))]
    tpDict = defaultdict(set)

    # Build a dictionary of matched tract and patch
    for tp in tpList:
        for tract in tp:
            # tract is tuple: t[0] is TractInfo(id=9216)
            # t[1] is PatchInfo(index=(2,4), ....)
            for patch in tract[1]:
                # tpDict[9216].add((2,4))
                tpDict[tract[0].getId()].add(patch.getIndex())
    print "number of patches ", sum([len(td) for td in tpDict.values()])

    visitDict = getVisits(tpDict[0], butler, filter=filter)

    return visitDict

