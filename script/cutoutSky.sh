#!/bin/sh

if [ $1 ]
then
    if [ $2 ]
    then 
        root=$2
    else 
        root='./'
    fi 

    coaddCutoutSky.py $1 -r $root --verbose --visual 

else
    echo "cutoutSky.sh prefix root"
fi
