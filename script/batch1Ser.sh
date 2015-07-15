#!/bin/sh

if [ $1 ] 
then 

    lis=$1

    round=0

    for i in `cat $1`; do 
        round=`echo " $round + 1" | bc`
        echo "### GALFIT RUN : $round"
        coaddCutout1Ser.py $i --mag=18.0 ;
    done

else
    echo "cutout1Ser.sh prefix"
fi
