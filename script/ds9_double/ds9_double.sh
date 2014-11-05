#!/bin/bash 

if [ $1 ] 
then 
    img1=$1
else 
    echo "ds9_double.sh img1 img2" 
fi 

if [ $2 ] 
then 
    img2=$2
else 
    echo "ds9_double.sh img1 img2" 
fi 

if [ $3 ] 
then 
    jpeg_name=$2
else 
    jpeg_name='double.jpeg' 
fi 

ds9 -nan black -view info no -view panner no -view magnifier no -view buttons no -view colorbar no -geometry 1800x1000 -scale asinh -zoom to fit $img1 -cmap HSV -regions load img1.reg -zoom to fit $img2 -cmap HSV -regions load img2.reg -match scale -match frame wcs -saveimage jpeg $jpeg_name 99 -tile column



