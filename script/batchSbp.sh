#!/bin/sh

if [ $1 ] 
then 

    lis=$1

    round=0

    for i in `cat $1`; do 
        round=`echo " $round + 1" | bc`
        echo "### ELLIPSE RUN NUMBER $round"
        coaddCutoutSbp.py $i --step=0.10 ;
    done

else
    echo "cutoutSbp.sh prefix"
fi
