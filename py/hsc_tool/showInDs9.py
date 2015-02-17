#!/usr/bin/env python

import argparse
import lsst.daf.persistence  as dafPersist
import lsst.afw.image        as afwImage
import lsst.afw.display.ds9  as ds9

def main(rootDir, visit, ccd, frame=1, title="", scale="zscale", zoom="to fit", trans=60, useEllipse=False):

    # make a butler and specify your dataId
    butler = dafPersist.Butler(rootDir)
    dataId = {'visit': visit, 'ccd':ccd}

    # get the exposure from the butler
    exposure = butler.get('calexp', dataId)

    # put the settings in a dict object and call ds9.mtv()
    settings = {'scale':scale, 'zoom': zoom, 'mask' : 'transparency %d' %(trans)}
    ds9.mtv(exposure, frame=frame, title=title, settings=settings)


    # now get the source catalog and overplot the points

    sources = butler.get('src', dataId)

    with ds9.Buffering():
        for i,source in enumerate(sources):
            color = ds9.RED
            size = 5.0
            
            if useEllipse:
                # show an ellipse symbol
                symbol = "@:{ixx},{ixy},{iyy}".format(ixx=source.getIxx(),
                                                      ixy=source.getIxy(),
                                                      iyy=source.getIyy())
            else:
                # just a simple point (symbols +, x, *, o are all accepted)
                symbol = "o"
                
            ds9.dot(symbol, source.getX(), source.getY(), ctype=color,
                    size=size, frame=frame, silent=True)
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("visit", type=int, help="Visit to show")
    parser.add_argument("ccd", type=int, help="CCD to show")
    parser.add_argument("-e", "--useEllipse", action='store_true', help="Overplot ellipses")
    parser.add_argument("-f", "--frame", type=int, default=1, help="Frame")
    parser.add_argument("-s", "--scale", default="zscale", help="Gray-scale")
    parser.add_argument("-t", "--title", default="", help="Figure title")
    parser.add_argument("-T", "--trans", default=60, help="Transparency")
    parser.add_argument("-z", "--zoom",  default="to fit", help="Zoom")
    args = parser.parse_args()

    main(args.root, args.visit, args.ccd,
         frame=args.frame, title=args.title, scale=args.scale, zoom=args.zoom, trans=args.trans,
         useEllipse=args.useEllipse
     )
