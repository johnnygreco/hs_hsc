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
            -c 1 --growH 2.5 --growW 3.5 ;

    done

else
    echo "batchPrep.sh list_prefix root"
fi
