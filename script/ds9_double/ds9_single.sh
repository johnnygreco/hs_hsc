#!/bin/bash 

img=$1
reg=$2

png_name="${1%.*}.png"

ds9 -nan black -view info no -view panner no -view magnifier no -view buttons no -view colorbar no -geometry 1000x1000 -scale asinh -scale mode 96.0 -zoom to fit $img -cmap heat -regions load $2 -zoom to fit -export png $png_name -exit  
