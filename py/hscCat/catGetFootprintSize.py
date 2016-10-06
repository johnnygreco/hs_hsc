#!/usr/bin/env python
# encoding: utf-8
"""Getting the footprint size for object in a catalog."""

import os
import argparse
import warnings
import numpy as np

# HSC Pipeline
import lsst.daf.persistence as dafPersist

# Astropy
from astropy.table import Table, Column


def getFootprintSize(catalog, root, filt='HSC-I', idCol='object_id',
                     tractCol='tract', patchCol='patch_s',
                     patchNCol='patch', catType='ref',
                     butler=None, verbose=True,
                     catFormat='fits', getParentSize=True):
    """Extract information about the footprint size."""
    # Read in the catalog
    data = Table.read(catalog, format=catFormat)
    # Sort the catalog based on Tract and Patch Number
    data.sort([tractCol, patchNCol])

    # Prepare the Butler
    if butler is None:
        try:
            butler = dafPersist.Butler(root)
            if verbose:
                print "## Load in the Butler"
        except:
            raise Exception("### Can not load the Butler")

    # Setup necessary data structure
    trPast, paPast = None, None
    footprintSize = []
    if getParentSize:
        footprintSizeParent = []

    # Loop through the catalog
    for obj in data:
        idUse = obj[idCol]
        trUse, paUse = int(obj[tractCol]), str(obj[patchCol]).strip()
        if (trUse != trPast) or (paUse != paPast):
            try:
                cat = butler.get('deepCoadd_' + catType,
                                 tract=trUse,
                                 patch=paUse,
                                 filter=filt,
                                 immediate=True)
                if verbose:
                    print("## Load in catalog for %d-%s" % (trUse, paUse))
            except Exception:
                if verbose:
                    print("!! Can not load catalog for %d-%s" % (trUse, paUse))
                cat = None

        if cat is not None:
            try:
                objFound = (cat[cat['id'] == idUse])[0]
                nPix = objFound.getFootprint().getArea()
                footprintSize.append(nPix)
            except Exception:
                if verbose:
                    print("!! Can not find object: %d-%s-%d" % (trUse, paUse,
                                                                idUse))
                objFound = None
                footprintSize.append(np.nan)
        else:
            objFound = None
            footprintSize.append(np.nan)

        if getParentSize:
            if (cat is not None) and (objFound is not None):
                parentID = objFound.getParent()
                if parentID != 0L:
                    try:
                        parentFound = (cat[cat['id'] == parentID])[0]
                        nPixParent = parentFound.getFootprint().getArea()
                        footprintSizeParent.append(nPixParent)
                    except Exception:
                        if verbose:
                            print("!! Can not find parent: %d-%s-%d" % (trUse,
                                                                        paUse,
                                                                        parentID))
                        footprintSizeParent.append(np.nan)
                else:
                    footprintSizeParent.append(np.nan)
            else:
                footprintSizeParent.append(np.nan)

        # Pass the current Tract, Patch to the trPast and paPast parameters
        trPast, paPast = trUse, paUse

    # Add new columns to the table
    data.add_column(Column(np.asarray(footprintSize),
                           name='footprintSize'))
    if getParentSize:
        data.add_column(Column(np.asarray(footprintSizeParent),
                                name='footprintSizeParent'))

    # Save the new catalog
    newCatalog = catalog.replace('.' + catFormat,
                                 '_fpSize.' + catFormat)
    data.write(newCatalog, format=catFormat, overwrite=True)


def testDavid():
    """Sanity check."""
    cat = 'hsc_cosmosfield_i_david_ID.fits'
    #cat = 'hsc_cosmosfield_i_david_ID_test.fits'
    root = '/lustre2/Subaru/SSP/rerun/S15Bi/udeep'
    #root = '/lustre2/HSC_DR/dr1/s15b/data/s15b_udeep'

    getFootprintSize(cat, root)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("catalog", help="Location and name of the catalogs")
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument('-f', '--filter', dest='filt', help="Filter",
                        default='HSC-I')
    parser.add_argument('-p', '--parent', action="store_true",
                        dest='getParentSize', default=False)
    parser.add_argument('-t', '--type', dest='catType',
                        help='Type of the HSC output catalog to use',
                        default='ref')

    args = parser.parse_args()

    getFootprintSize(args.catalog, args.root,
                     filt=args.filt, catType=args.catType,
                     getParentSize=args.getParentSize)
