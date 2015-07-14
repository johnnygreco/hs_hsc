#!/bin/sh

if [ $1 ] 
then 

    prefix=$1

    coaddCutoutSbp.py $prefix --step=0.10

else
    echo "cutoutSbp.sh prefix"
fi
