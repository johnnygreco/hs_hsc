#!/bin/sh

if [ $1 ] 
then 

    lis=$1

    for i in `cat $1`; do 
        coaddCutoutSbp.py $i --step=0.10 ;
    done

else
    echo "cutoutSbp.sh prefix"
fi
