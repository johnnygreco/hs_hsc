#!/bin/sh

if [ $1 ]
then
    if [ $2 ]
    then 
        root=$2
    else 
        root='./'
    fi 

    for i in `cat $1`; do 
        
        coaddCutoutSky.py $i -r $root --verbose --visual ;

    done

else
    echo "batchSky.sh list_prefix root"
fi
