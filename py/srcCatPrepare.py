#!/usr/bin/env python

import argparse
import numpy
import matplotlib.pyplot as pyplot
import lsst.daf.persistence as dafPersist

def main(rootDir, tract, patch, filter):

    # Make a butler and specify the dataId
    butler = dafPersist.Butler(rootDir)
    dataId = {'tract': tract, 'patch': patch, 'filter': filter}

    # Get the exposure from the butler
    exposure = butler.get('deepCoadd_calexp')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("tract", type=int, help="Visit to show")
    parser.add_argument("patch",  help="CCD to show")
    parser.add_argument("filter", help="Filter")
    args = parser.parse_args()

    main(args.root, args.tract, args.patch, args.filter)
