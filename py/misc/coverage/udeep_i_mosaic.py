#!/usr/bin/env python

import atpy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
from matplotlib.patches import Polygon


tbl = atpy.Table( 'udeep_i_mosaic.fits' )
num = len( tbl.llcra )

fig = plt.figure( figsize=(12,10), dpi=70 )
ax1 = fig.add_subplot( 1, 1, 1 )

ax1.set_xlim( 148.9, 151.18 )
ax1.set_ylim(  1.31,   3.09 )
#ax1.grid( True )

loc1 = plticker.MultipleLocator( base=0.2 )
ax1.xaxis.set_major_locator( loc1 )
loc2 = plticker.MultipleLocator( base=0.1 )
ax1.yaxis.set_major_locator( loc2 )
ax1.tick_params( 'both', length=6, width=1.2, which='major' )

plt.xlabel(r'RA (J2000)',  fontsize=25, labelpad=20 )
plt.ylabel(r'Dec (J2000)', fontsize=25, labelpad=20 )
plt.title(r'HSC-UDEEP Mosaic Patches', fontsize=32, y=1.02 )


for ii in range( num ):

    ax1.add_patch( Polygon( [ [tbl.llcra[ii], tbl.llcdecl[ii]],
                            [tbl.lrcra[ii], tbl.lrcdecl[ii]],
                            [tbl.urcra[ii], tbl.urcdecl[ii]],
                            [tbl.ulcra[ii], tbl.ulcdecl[ii]],
                            [tbl.llcra[ii], tbl.llcdecl[ii]] ],
                            closed=True, fill=True, alpha=0.3, color='Blue' ) )

ax1.add_patch( Polygon( [ [tbl.llcra[0], tbl.llcdecl[0]],
                        [tbl.lrcra[0], tbl.lrcdecl[0]],
                        [tbl.urcra[0], tbl.urcdecl[0]],
                        [tbl.ulcra[0], tbl.ulcdecl[0]],
                        [tbl.llcra[0], tbl.llcdecl[0]] ],
                        closed=True, fill=False, color='Red', linewidth=3 ) )

plt.subplots_adjust( left=0.10, bottom=0.10, top=0.92, right=0.99 )

fig.savefig( "udeep_i_mosaic.png", dpi=120 )
plt.show()
