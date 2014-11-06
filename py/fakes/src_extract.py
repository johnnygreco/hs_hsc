#!/usr/bin/env python

import argparse
import numpy
import lsst.daf.persistence as dafPersist

def printInfo( butler, dataId ):

    src = butler.get( "deepCoadd_src", dataId )
    num_obj = len( src )

    # Misc
    comma = " , "

    # Magnitude and Error
    # PSF
    mag_psf = ( 27.0 - 2.5 * numpy.log10( src.get( 'flux.psf' )  ) )
    err_psf = ( 2.5 / numpy.log(10)*( src.get( 'flux.psf.err' ) /
                                    src.get( 'flux.psf' ) ) )
    # CModel
    mag_mod = ( 27.0 - 2.5 * numpy.log10( src.get( 'cmodel.flux' )  ) )
    err_mod = ( 2.5 / numpy.log(10)*( src.get( 'cmodel.flux.err' ) /
                                    src.get( 'cmodel.flux' ) ) )
    # Kron
    mag_kro = ( 27.0 - 2.5 * numpy.log10( src.get( 'flux.kron' ) ) )
    err_kro = ( 2.5 / numpy.log(10)*( src.get( 'flux.kron.err' ) /
                                    src.get( 'flux.kron' ) ) )

    # Kron Radius
    r_kron  = src.get( 'flux.kron.radius' )
    # FracDev
    fracdev = src.get( 'cmodel.fracDev' )

    # Quality Control
    extend  = src.get( 'classification.extendedness' )
    inner   = src.get( 'detect.is-patch-inner' )
    primary = src.get( 'detect.is-primary' )

    for i in range( num_obj ):

        id  = src[i].get( 'id' )
        ra  = src[i].getRa().asDegrees()
        dec = src[i].getDec().asDegrees()

        print id, comma, ra, comma, dec, comma, \
              args.patch, comma, args.tract, comma, \
              src[i].get( 'deblend.nchild' ),   \
              comma, mag_psf[i], comma, mag_mod[i], comma, mag_kro[i], \
              comma, err_psf[i], comma, err_mod[i], comma, err_kro[i], \
              comma, extend[i],  comma, r_kron[i],  comma, fracdev[i], \
              comma, inner[i],   comma, primary[i],    \
              comma, src[i].get( 'flags.negative' ),   \
              comma, src[i].get( 'flags.pixel.edge' ), \
              comma, src[i].get( 'flags.pixel.interpolated.any' ), \
              comma, src[i].get( 'flags.pixel.saturated.any' ),    \
              comma, src[i].get( 'flags.pixel.cr.any' ),           \
              comma, src[i].get( 'flags.pixel.bad' ),              \
              comma, src[i].get( 'flags.pixel.suspect.any' ),      \
              comma, src[i].get( 'flux.aperture.flags' ),          \
              comma, src[i].get( 'flux.kron.flags' ),              \
              comma, src[i].get( 'flux.psf.flags' ),               \
              comma, src[i].get( 'cmodel.flux.flags' ),            \
              comma, src[i].get( 'flags.badcentroid' ),            \
              comma, src[i].get( 'centroid.sdss.flags' ),          \
              comma, src[i].get( 'parent' ),                       \
              comma, src[i].get( 'correctfluxes.apcorr.flags' ),   \
              comma, src[i].get( 'shape.sdss.flags' ),             \
              comma, src[i].get( 'shape.sdss.flags.shift' ),       \
              comma, src[i].get( 'shape.sdss.flags.maxiter' )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("tract",  type=int, help="TRACT to get")
    parser.add_argument("patch",  type=str, help="PATCH to get")
    parser.add_argument("filter", type=str, help="FILTER to get")
    args = parser.parse_args()

    dataDir = "/lustre/HSC_DR/rerun/ssp_s14a0_udeep_20140523a/"
    #dataDir = "/lustre/HSC_DR/rerun/ssp_s14a0_wide_20140523a/"

    butler  = dafPersist.Butler( dataDir )

    dataId = {"tract" : args.tract, "patch": args.patch, 'filter':args.filter}

    printInfo(butler, dataId)
