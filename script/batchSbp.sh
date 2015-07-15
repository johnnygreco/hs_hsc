#!/bin/sh

if [ $1 ] 
then 

    lis=$1

    round=0

    for i in `cat $1`; do 
        round=`echo " $round + 1" | bc`
        echo "### ELLIPSE RUN NUMBER $round"

        coaddCutoutSbp.py $i --step 0.12 --fracBad 0.70 \
            --lowClip 3.0 --uppClip 2.0 --nClip 2 --olthresh 0.5 \
            --intMode median --minIt 10 --maxIt 100 ; 
    done

else
    echo "batchSbp.sh list"
fi
