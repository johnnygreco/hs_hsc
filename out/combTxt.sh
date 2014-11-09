#!/bin/bash 

touch $2
for i in `cat $1`; do 

    cat $i >> $2; 

done
