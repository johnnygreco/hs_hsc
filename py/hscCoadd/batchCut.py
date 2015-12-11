#!/usr/bin/env python
# encoding: utf-8

import os
import argparse
import coaddBatchCutout as cbc

COM = '#' * 100
SEP = '-' * 100
WAR = '!' * 100


def run(args):
    """Run cutout in batch mode."""
    if os.path.isfile(args.incat):
        cbc.coaddBatchCutFull(args.root, args.incat,
                              filter=args.filter,
                              idField=args.idField,
                              prefix=args.prefix,
                              zCutoutSize=args.zCutout,
                              zField=args.zField,
                              onlyColor=args.onlyColor,
                              noColor=args.noColor,
                              saveSrc=args.saveSrc,
                              makeDir=args.makeDir,
                              raField=args.raField,
                              decField=args.decField)
    else:
        raise Exception("### Can not find the input catalog: %s" % args.incat)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("incat", help="The input catalog for cutout")
    parser.add_argument("-s", '--size', dest='size', type=int,
                        help="Half size of the cutout box", default=200)
    parser.add_argument('-f', '--filter', dest='filter', help="Filter",
                        default='HSC-I')
    parser.add_argument('-cf', '--color-filters', dest='colorFilters',
                        help="Choice of filters for color images",
                        default='riz')
    parser.add_argument('-sf', '--size-field', dest='sizeField',
                        help="Column name for cutout size",
                        default='cutout_size')
    parser.add_argument('-info1', '--infoField1', dest='infoField1',
                        help="Column name for first extra information",
                        default='MSTAR')
    parser.add_argument('-info2', '--infoField2', dest='infoField2',
                        help="Column name for second extra information",
                        default='VDISP')
    parser.add_argument('-oc', '--onlyColor', action="store_true",
                        dest='onlyColor', default=False)
    parser.add_argument('-safe', '--safe', action="store_true", dest='safe',
                        default=False)
    parser.add_argument('-clean', '--clean', action="store_true", dest='clean',
                        default=False)
    parser.add_argument('-v', '--verbose', action="store_true", dest='verbose',
                        default=False)
    parser.add_argument('-src', '--src', action="store_true", dest='saveSrc',
                        default=True)
    parser.add_argument('-makeDir', '--makeDir', action="store_true",
                        dest='makeDir', default=True)
    parser.add_argument('-zc', '--zCutoutSize', action="store_true",
                        dest='zCutout', default=True)
    parser.add_argument('-nc', '--noColor', action="store_true",
                        dest='noColor', default=True)
    parser.add_argument('-p', '--prefix', dest='prefix',
                        help='Prefix of the output file',
                        default='sl')
    parser.add_argument('-id', '--id', dest='idField',
                        help="Column name for ID", default='ID')
    parser.add_argument('-ra', '--ra', dest='raField',
                        help="Column name for RA", default='RA')
    parser.add_argument('-dec', '--dec', dest='decField',
                        help="Column name for DEC", default='DEC')
    parser.add_argument('-z', '--redshift', dest='zField',
                        help="Column name for z", default='Z')
    args = parser.parse_args()

    run(args)
