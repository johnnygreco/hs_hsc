#!/bin/sh

if [ $1 ] 
then 

    prefix=$1

    coaddCutout1Ser.py $prefix --mag=18.0

else
    echo "cutout1Ser.sh prefix"
fi
