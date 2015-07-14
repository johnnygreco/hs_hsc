prefix=$1

galSBP.py $prefix"_psf.fits" --step=0.08 --galQ=0.95 --galPA=0.0 --stage=3 --iniSma=5.0 --maxSma=28.0 --outThre=1e-6
