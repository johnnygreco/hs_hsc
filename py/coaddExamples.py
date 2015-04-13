# Examples:

# Example galaxies in the WIDE field
# From redmapper ID=16474; z=0.38
wideRa  = 133.02153
wideDec = 1.3656108

# Exampe galaxies in the UDEEP field
# From COSMOS campact galaxies: ID=6026402
udeepRa  = 149.75626
udeepDec = 2.794748

# Get the corner (Ra,Dec) for a list of images

from coaddImgCornerRaDec import *
multi = coaddAllBandsCorners('/lustre/Subaru/SSP/rerun/yasuda/SSP3.4.1_20141224_cosmos/',prefix='ssp341_cosmos')
