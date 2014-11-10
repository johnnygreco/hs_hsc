#!/bin/bash 

prefix=$1
touch $2

for i in `ls $prefix*.txt`; do 

    cat $i >> $2; 

done
