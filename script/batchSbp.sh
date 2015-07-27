#!/bin/sh

if [ $1 ] 
then 

    lis=$1

    round=0

    for i in `cat $1`; do 
        round=`echo " $round + 1" | bc`
        echo "    "
        echo "### ELLIPSE RUN NUMBER $round"
        echo "    "

        coaddCutoutSbp.py $i --step 0.10 --fracBad 0.50 \
            --lowClip 3.0 --uppClip 3.0 --nClip 2 --olthresh 0.3 \
            --intMode median --minIt 10 --maxIt 160 --maxTry 6 \
            --outRatio 1.25 ; 
    done

else
    echo "batchSbp.sh list"
fi
