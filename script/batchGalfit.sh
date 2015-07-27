#!/bin/sh

if [ $1 ] 
then 

    for i in `cat $1`; do 
        coaddCutoutGalfitSimple.py $i --mag=18.0 --ser2Comp --ser3Comp \
            --skyGrad --run1 --constrCen ;
    done

else
    echo "batchGalfit.sh prefix"
fi
