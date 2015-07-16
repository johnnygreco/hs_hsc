#!/bin/sh

if [ $1 ] 
then 

    for i in `cat $1`; do 
        coaddCutoutGalfitSimple.py $i --mag=18.0 --ser2Comp --ser3Comp \
            --skyGrad --run ;
    done

else
    echo "cutout1Ser.sh prefix"
fi
