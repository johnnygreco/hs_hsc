#!/bin/bash 

imfit redBCG_127_HSC-I_full_img.fits \
    --psf redBCG_127_HSC-I_full_psf.fits \
    --mask redBCG_127_HSC-I_full_mskfin.fits \
    --noise redBCG_127_HSC-I_full_sig.fits \
    -c redBCG_127_HSC-I_full_2ser.imfit \
    --save-params redBCG_127_HSC-I_full_2ser.out \
    --save-residual redBCG_127_HSC-I_full_2ser_imfit_res.fits \
    --nm
