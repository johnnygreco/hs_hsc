#!/bin/sh

if [ $1 ]
then
    if [ $2 ]
    then 
        root=$2
    else 
        root='./'
    fi 
    coaddCutoutPrepare.py $1 -r $root \
        -c 1 --bkgH 8 --bkgC 80 \
        --thrH 2.5 --thrC 1.2 \
        --growC 7.0 --growW 5.0 --growH 1.9 \
        --debConC 0.015 --debConH 0.004 \
        --debThrC 32.0 --debThrH 16.0
else
    echo "cutoutPrep.sh prefix root"
fi
