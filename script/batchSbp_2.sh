#!/bin/sh

if [ $1 ] 
then 

    lis=$1

    round=0

    for i in `cat $1`; do 
        round=`echo " $round + 1" | bc`
        echo "### ELLIPSE RUN NUMBER $round"

        coaddCutoutSbp.py $i --step 0.15 --fracBad 0.70 \
            --lowClip 3.0 --uppClip 3.0 --nClip 3 --olthresh 0.3 \
            --intMode median --minIt 10 --maxIt 200 --maxTry 4 \
            --outRatio 1.2 ; 
    done

else
    echo "batchSbp.sh list"
fi
