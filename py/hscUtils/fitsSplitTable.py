#!/usr/bin/env python
# encoding: utf-8

from __future__ import division, print_function

import os
import argparse
import warnings

import numpy as np
from astropy.table import Table


def run(args):
    """
    Split a FITS table into roughly equal chunks.

    Parameters:
        args.file :  FITS table
        args.num  :  Number of output tables
    """
    file = str(args.file).strip()
    num = int(args.num)

    if not os.path.isfile(file):
        raise Exception("# Can not find the input FITS table: %s" % file)
    tab = Table.read(file, format='fits')
    prefix = file.replace('.fits', '')
    if not num >= len(tab):
        start, end = 0, num
        for index in range(np.ceil(len(tab)/num).astype(int)):
            suffix = str(index+1).strip()
            output = prefix + '_' + suffix + '.fits'
            if args.verbose:
                print("# Saving the %d table: %s" % ((index+1), output))
            if end >= len(tab):
                end = (len(tab) - 1)
            chunk = tab[start:end]
            chunk.write(output, format='fits', overwrite=True)
            start = end
            end = (end + num)
    else:
        warnings.warn("# Chunk size is larger or equal ",
                      "to the size of the input. Do Nothing!")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='Fits file')
    parser.add_argument('num', help='Number of output catalog', type=int)
    parser.add_argument('-v', '--verbose', help='Blah..Blah..',
                        default=False, action='store_true', dest='verbose')
    args = parser.parse_args()

    run(args)
