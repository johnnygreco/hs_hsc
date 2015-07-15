#!/bin/sh

if [ $1 ] 
then 

    prefix=$1

    coaddCutoutGalfitSimple.py $prefix --mag=18.0 --ser2Comp --ser3Comp \
        --skyGrad
    # --run --useF1 --useF4

else
    echo "cutout1Ser.sh prefix"
fi
