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
        coaddCutoutPrepare.py $i -r $root \
            -c 1 --bkgH 8 --bkgC 80 \
            --thrH 2.5 --thrC 1.2 \
            --growC 6.8 --growW 4.0 --growH 1.8 \
            --debConC 0.015 --debConH 0.004 \
            --debThrC 32.0 --debThrH 16.0 ;
    done

else
    echo "batchPrep.sh list_prefix root"
fi
