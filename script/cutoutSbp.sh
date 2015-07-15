#!/bin/sh

if [ $1 ] 
then 

    prefix=$1

    coaddCutoutSbp.py $prefix --step 0.12 --fracBad 0.70 \
        --lowClip 3.0 --uppClip 2.0 --nClip 3 --olthresh 0.5 \
        --intMode median --minIt 10 --maxIt 100 

else
    echo "cutoutSbp.sh prefix"
fi
