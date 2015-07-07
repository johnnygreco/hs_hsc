pro hsc_filter_diff

    ;--------------------------------------------------------------------------;
    ; Read in the magnitude file 
    mag1k = read_mags_hs( 'ssp_miles_z3_kro.out.mags' )
    ;mwrfits, mag1k, 'ssp_miles_z3_kro_mags.fits', /create
    mag1s = read_mags_hs( 'ssp_miles_z3_sal.out.mags' )
    ;mwrfits, mag1k, 'ssp_miles_z3_sal_mags.fits', /create
    mag2k = read_mags_hs( 'ssp_miles_z4_kro.out.mags' )
    ;mwrfits, mag1k, 'ssp_miles_z4_kro_mags.fits', /create
    mag2s = read_mags_hs( 'ssp_miles_z4_sal.out.mags' )
    ;mwrfits, mag1k, 'ssp_miles_z4_sal_mags.fits', /create
    mag3k = read_mags_hs( 'ssp_miles_z5_kro.out.mags' )
    ;mwrfits, mag1k, 'ssp_miles_z5_kro_mags.fits', /create
    mag3s = read_mags_hs( 'ssp_miles_z5_sal.out.mags' )
    ;mwrfits, mag1k, 'ssp_miles_z4_sal_mags.fits', /create
    ;--------------------------------------------------------------------------;
    psxsize  = 28 
    psysize  = 20
    position = [ 0.15, 0.15, 0.992, 0.992 ]
    ;--------------------------------------------------------------------------;

    ;--------------------------------------------------------------------------;
    ; HSC-I band i_hsc 
    ; Compared with i_sdss, i_suprime, i_megacam, i_old_megacam, i_des, 
    ;               f814w_acs, f814w_wfc3
    ; Using g_hsc - i_hsc; i_hsc - z_hsc as x-axis 
    ;--------------------------------------------------------------------------;

    ;--------------------------------------------------------------------------;
    iplot_1 = 'hsc_i_filter_compare_1.eps'

    hsc_ri_a = ( mag2s.r_hsc - mag2s.i_hsc )
    hsc_iz_a = ( mag2s.i_hsc - mag2s.z_hsc )
    i_diff1a = ( mag2s.i_hsc - mag2s.i_sdss ) 
    i_diff2a = ( mag2s.i_hsc - mag2s.i_suprime ) 
    i_diff3a = ( mag2s.i_hsc - mag2s.i_megacam )
    i_diff4a = ( mag2s.i_hsc - mag2s.i_old_megacam )
    i_diff5a = ( mag2s.i_hsc - mag2s.i_des )
    i_diff6a = ( mag2s.i_hsc - mag2s.f814w_acs )
    i_diff7a = ( mag2s.i_hsc - mag2s.f814w_wfc3 )

    hsc_ri_b = ( mag1s.r_hsc - mag1s.i_hsc )
    hsc_iz_b = ( mag1s.i_hsc - mag1s.z_hsc )
    i_diff1b = ( mag1s.i_hsc - mag1s.i_sdss ) 
    i_diff2b = ( mag1s.i_hsc - mag1s.i_suprime ) 
    i_diff3b = ( mag1s.i_hsc - mag1s.i_megacam )
    i_diff4b = ( mag1s.i_hsc - mag1s.i_old_megacam )
    i_diff5b = ( mag1s.i_hsc - mag1s.i_des )
    i_diff6b = ( mag1s.i_hsc - mag1s.f814w_acs )
    i_diff7b = ( mag1s.i_hsc - mag1s.f814w_wfc3 )

    hsc_ri_c = ( mag3s.r_hsc - mag3s.i_hsc )
    hsc_iz_c = ( mag3s.i_hsc - mag3s.z_hsc )
    i_diff1c = ( mag3s.i_hsc - mag3s.i_sdss ) 
    i_diff2c = ( mag3s.i_hsc - mag3s.i_suprime ) 
    i_diff3c = ( mag3s.i_hsc - mag3s.i_megacam )
    i_diff4c = ( mag3s.i_hsc - mag3s.i_old_megacam )
    i_diff5c = ( mag3s.i_hsc - mag3s.i_des )
    i_diff6c = ( mag3s.i_hsc - mag3s.f814w_acs )
    i_diff7c = ( mag3s.i_hsc - mag3s.f814w_wfc3 )

    ;--------------------------------------------------------------------------;
    mydevice = !d.name 
    !p.font=1
    set_plot, 'ps' 
    device, filename=iplot_1, font_size=9.0, /encapsulated, $
        /color, set_font='TIMES-ROMAN', /bold, xsize=psxsize, ysize=psysize

    xmin = min( [ min( hsc_iz_a ), min( hsc_iz_b ), min( hsc_iz_c ) ] )
    xmax = max( [ max( hsc_iz_a ), max( hsc_iz_b ), max( hsc_iz_c ) ] )
    xrange = [ xmin - 0.02, xmax + 0.02 ]
    ;xrange = [ min( alog10( mag2s.agegyr ) ) - 0.05, $
    ;           max( alog10( mag2s.agegyr ) ) + 0.10 ]

    ymin = min( [ min( i_diff1a ), min( i_diff2a ), min( i_diff3a ), $
                  min( i_diff4a ), min( i_diff5a ), min( i_diff6a ), $
                  min( i_diff7a ) ] )
    ymax = max( [ max( i_diff1a ), max( i_diff2a ), max( i_diff3a ), $
                  max( i_diff4a ), max( i_diff5a ), max( i_diff6a ), $
                  max( i_diff7a ) ] )
    yrange = [ ymin - 0.005, ymax + 0.005 ]

    index_a = sort( hsc_iz_a )
    index_b = sort( hsc_iz_b )
    index_c = sort( hsc_iz_c )

    cgPlot, hsc_iz_a[ index_a ], i_diff1a[ index_a ], $
        xstyle=1, ystyle=1, linestyle=0, color=cgColor( 'Dark Gray' ), $
        thick=3.0, charsize=3.5, charthick=8.0, xthick=10.0, ythick=10.0, $
        /noerase, position=position, yminor=5, $
        /nodata, xrange=xrange, yrange=yrange, $
        xticklen=0.025, yticklen=0.02, $
        xtitle=TeXToIDL( 'i_{HSC} - z_{HSC} ') + '(mag)', $
        ytitle=TeXToIDL( 'i_{HSC} - i_{Other} ') + '(mag)'

    cgOPlot, !X.Crange, [ 0.0, 0.0 ], linestyle=2, thick=5.0, $
        color=cgColor( 'Dark Gray' )

    ;cgOPlot, mag2s.agegyr, i_diff1a, $
    cgOPlot, hsc_iz_a[ index_a ], i_diff1a[ index_a ], $
        linestyle=0, color=cgColor( 'RED6' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_iz_a[ index_a ], i_diff3a[ index_a ], $
        linestyle=0, color=cgColor( 'GRN6' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_iz_a[ index_a ], i_diff4a[ index_a ], $
        linestyle=0, color=cgColor( 'BLK5' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_iz_a[ index_a ], i_diff5a[ index_a ], $
        linestyle=0, color=cgColor( 'TAN6' ), psym=16, thick=3.0, symsize=1.2
    cgOPlot, hsc_iz_a[ index_a ], i_diff6a[ index_a ], $
        linestyle=0, color=cgColor( 'Pink' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_iz_a[ index_a ], i_diff7a[ index_a ], $
        linestyle=0, color=cgColor( 'BLU3' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_iz_a[ index_a ], i_diff2a[ index_a ], $
        linestyle=0, color=cgColor( 'BLU7' ), psym=16, thick=3.0, symsize=1.2 

    cgOPlot, hsc_iz_b[ index_b ], i_diff1b[ index_b ], $
        linestyle=0, color=cgColor( 'RED6' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_iz_b[ index_b ], i_diff3b[ index_b ], $
        linestyle=0, color=cgColor( 'GRN6' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_iz_b[ index_b ], i_diff4b[ index_b ], $
        linestyle=0, color=cgColor( 'BLK5' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_iz_b[ index_b ], i_diff5b[ index_b ], $
        linestyle=0, color=cgColor( 'TAN6' ), psym=9, thick=3.0, symsize=1.5
    cgOPlot, hsc_iz_b[ index_b ], i_diff6b[ index_b ], $
        linestyle=0, color=cgColor( 'Pink' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_iz_b[ index_b ], i_diff7b[ index_b ], $
        linestyle=0, color=cgColor( 'BLU3' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_iz_b[ index_b ], i_diff2b[ index_b ], $
        linestyle=0, color=cgColor( 'BLU7' ), psym=9, thick=3.0, symsize=1.5 

    cgOPlot, hsc_iz_c[ index_c ], i_diff1c[ index_c ], $
        linestyle=0, color=cgColor( 'RED6' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_iz_c[ index_c ], i_diff3c[ index_c ], $
        linestyle=0, color=cgColor( 'GRN6' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_iz_c[ index_c ], i_diff4c[ index_c ], $
        linestyle=0, color=cgColor( 'BLK5' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_iz_c[ index_c ], i_diff5c[ index_c ], $
        linestyle=0, color=cgColor( 'TAN6' ), psym=4, thick=3.0, symsize=1.5
    cgOPlot, hsc_iz_c[ index_c ], i_diff6c[ index_c ], $
        linestyle=0, color=cgColor( 'Pink' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_iz_c[ index_c ], i_diff7c[ index_c ], $
        linestyle=0, color=cgColor( 'BLU3' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_iz_c[ index_c ], i_diff2c[ index_c ], $
        linestyle=0, color=cgColor( 'BLU7' ), psym=4, thick=3.0, symsize=1.5 

    cgPlots, 0.18, 0.95, psym=16, symsize=2.5, color=cgColor( 'BLU7' ), /normal
    cgPlots, 0.20, 0.95, psym=9,  symsize=2.5, color=cgColor( 'BLU7' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.95, psym=4,  symsize=2.5, color=cgColor( 'BLU7' ), $
        /normal, thick=5
    cgText,  0.25, 0.94, 'SuprimeCam i', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.18, 0.91, psym=16, symsize=2.5, color=cgColor( 'RED6' ), /normal
    cgPlots, 0.20, 0.91, psym=9,  symsize=2.5, color=cgColor( 'RED6' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.91, psym=4,  symsize=2.5, color=cgColor( 'RED6' ), $
        /normal, thick=5
    cgText,  0.25, 0.90, 'SDSS i', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.18, 0.87, psym=16, symsize=2.5, color=cgColor( 'GRN6' ), /normal
    cgPlots, 0.20, 0.87, psym=9,  symsize=2.5, color=cgColor( 'GRN6' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.87, psym=4,  symsize=2.5, color=cgColor( 'GRN6' ), $
        /normal, thick=5
    cgText,  0.25, 0.86, 'MegaCam i', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.18, 0.83, psym=16, symsize=2.5, color=cgColor( 'BLK5' ), /normal
    cgPlots, 0.20, 0.83, psym=9,  symsize=2.5, color=cgColor( 'BLK5' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.83, psym=4,  symsize=2.5, color=cgColor( 'BLK5' ), $
        /normal, thick=5
    cgText,  0.25, 0.82, 'MegaCam i (old)', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.47, 0.95, psym=16, symsize=2.5, color=cgColor( 'TAN6' ), /normal
    cgPlots, 0.49, 0.95, psym=9,  symsize=2.5, color=cgColor( 'TAN6' ), $
        /normal, thick=5
    cgPlots, 0.51, 0.95, psym=4,  symsize=2.5, color=cgColor( 'TAN6' ), $
        /normal, thick=5
    cgText,  0.54, 0.94, 'DES i', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.47, 0.91, psym=16, symsize=2.5, color=cgColor( 'PINK' ), /normal
    cgPlots, 0.49, 0.91, psym=9,  symsize=2.5, color=cgColor( 'PINK' ), $
        /normal, thick=5
    cgPlots, 0.51, 0.91, psym=4,  symsize=2.5, color=cgColor( 'PINK' ), $
        /normal, thick=5
    cgText,  0.54, 0.90, 'ACS/HST F814W', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.47, 0.87, psym=16, symsize=2.5, color=cgColor( 'BLU3' ), /normal
    cgPlots, 0.49, 0.87, psym=9,  symsize=2.5, color=cgColor( 'BLU3' ), $
        /normal, thick=5
    cgPlots, 0.51, 0.87, psym=4,  symsize=2.5, color=cgColor( 'BLU3' ), $
        /normal, thick=5
    cgText,  0.54, 0.86, 'WFC3/HST F814W', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlot, hsc_iz_a[ index_a ], i_diff1a[ index_a ], $
        xstyle=1, ystyle=1, linestyle=0, color=cgColor( 'Dark Gray' ), $
        thick=3.0, charsize=3.5, charthick=8.0, xthick=10.0, ythick=10.0, $
        /noerase, position=position, yminor=5, $
        /nodata, xrange=xrange, yrange=yrange, $
        xticklen=0.025, yticklen=0.02

    device, /close 
    set_plot, mydevice 
    ;--------------------------------------------------------------------------;

    ;--------------------------------------------------------------------------;
    ; HSC-Z band z_hsc 
    ; Compared with z_sdss, z_suprime, z_megacam, z_old_suprime, z_des, 
    ;               f850lp_acs, f850lp_wfc3
    ; Using i_hsc - z_hsc; z_hsc - y_hsc as x-axis 
    ;--------------------------------------------------------------------------;

    ;--------------------------------------------------------------------------;
    zplot_1 = 'hsc_z_filter_compare_1.eps'

    hsc_zy_a = ( mag2s.z_hsc - mag2s.y_hsc )
    hsc_iz_a = ( mag2s.i_hsc - mag2s.z_hsc )
    z_diff1a = ( mag2s.z_hsc - mag2s.z_sdss ) 
    z_diff2a = ( mag2s.z_hsc - mag2s.z_suprime ) 
    z_diff3a = ( mag2s.z_hsc - mag2s.z_megacam )
    z_diff4a = ( mag2s.z_hsc - mag2s.z_old_suprime )
    z_diff5a = ( mag2s.z_hsc - mag2s.z_des )
    z_diff6a = ( mag2s.z_hsc - mag2s.f850lp_acs )
    z_diff7a = ( mag2s.z_hsc - mag2s.f850lp_wfc3 )

    hsc_zy_b = ( mag1s.z_hsc - mag1s.y_hsc )
    hsc_iz_b = ( mag1s.i_hsc - mag1s.z_hsc )
    z_diff1b = ( mag1s.z_hsc - mag1s.z_sdss ) 
    z_diff2b = ( mag1s.z_hsc - mag1s.z_suprime ) 
    z_diff3b = ( mag1s.z_hsc - mag1s.z_megacam )
    z_diff4b = ( mag1s.z_hsc - mag1s.z_old_suprime )
    z_diff5b = ( mag1s.z_hsc - mag1s.z_des )
    z_diff6b = ( mag1s.z_hsc - mag1s.f850lp_acs )
    z_diff7b = ( mag1s.z_hsc - mag1s.f850lp_wfc3 )

    hsc_zy_c = ( mag3s.z_hsc - mag3s.y_hsc )
    hsc_iz_c = ( mag3s.i_hsc - mag3s.z_hsc )
    z_diff1c = ( mag3s.z_hsc - mag3s.z_sdss ) 
    z_diff2c = ( mag3s.z_hsc - mag3s.z_suprime ) 
    z_diff3c = ( mag3s.z_hsc - mag3s.z_megacam )
    z_diff4c = ( mag3s.z_hsc - mag3s.z_old_suprime )
    z_diff5c = ( mag3s.z_hsc - mag3s.z_des )
    z_diff6c = ( mag3s.z_hsc - mag3s.f850lp_acs )
    z_diff7c = ( mag3s.z_hsc - mag3s.f850lp_wfc3 )

    ;--------------------------------------------------------------------------;
    mydevice = !d.name 
    !p.font=1
    set_plot, 'ps' 
    device, filename=zplot_1, font_size=9.0, /encapsulated, $
        /color, set_font='TIMES-ROMAN', /bold, xsize=psxsize, ysize=psysize

    xmin = min( [ min( hsc_zy_a ), min( hsc_zy_b ), min( hsc_zy_c ) ] )
    xmax = max( [ max( hsc_zy_a ), max( hsc_zy_b ), max( hsc_zy_c ) ] )
    xrange = [ xmin - 0.01, xmax + 0.04 ]
    ;xrange = [ min( alog10( mag2s.agegyr ) ) - 0.05, $
    ;           max( alog10( mag2s.agegyr ) ) + 0.10 ]

    ymin = min( [ min( z_diff1a ), min( z_diff2a ), min( z_diff3a ), $
                  min( z_diff4a ), min( z_diff5a ), min( z_diff6a ), $
                  min( z_diff7a ) ] )
    ymax = max( [ max( z_diff1a ), max( z_diff2a ), max( z_diff3a ), $
                  max( z_diff4a ), max( z_diff5a ), max( z_diff6a ), $
                  max( z_diff7a ) ] )
    yrange = [ ymin - 0.005, ymax + 0.01 ]

    index_a = sort( hsc_zy_a )
    index_b = sort( hsc_zy_b )
    index_c = sort( hsc_zy_c )

    cgPlot, hsc_zy_a[ index_a ], z_diff1a[ index_a ], $
        xstyle=1, ystyle=1, linestyle=0, color=cgColor( 'Dark Gray' ), $
        thick=3.0, charsize=3.5, charthick=8.0, xthick=10.0, ythick=10.0, $
        /noerase, position=position, yminor=5, $
        /nodata, xrange=xrange, yrange=yrange, $
        xticklen=0.025, yticklen=0.02, $
        xtitle=TeXToIDL( 'z_{HSC} - y_{HSC} ') + '(mag)', $
        ytitle=TeXToIDL( 'z_{HSC} - z_{Other} ') + '(mag)'

    cgOPlot, !X.Crange, [ 0.0, 0.0 ], linestyle=2, thick=5.0, $
        color=cgColor( 'Dark Gray' )

    ;cgOPlot, mag2s.agegyr, z_diff1a, $
    cgOPlot, hsc_zy_a[ index_a ], z_diff1a[ index_a ], $
        linestyle=0, color=cgColor( 'RED6' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_zy_a[ index_a ], z_diff3a[ index_a ], $
        linestyle=0, color=cgColor( 'GRN6' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_zy_a[ index_a ], z_diff4a[ index_a ], $
        linestyle=0, color=cgColor( 'BLK5' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_zy_a[ index_a ], z_diff5a[ index_a ], $
        linestyle=0, color=cgColor( 'TAN6' ), psym=16, thick=3.0, symsize=1.2
    cgOPlot, hsc_zy_a[ index_a ], z_diff6a[ index_a ], $
        linestyle=0, color=cgColor( 'Pink' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_zy_a[ index_a ], z_diff7a[ index_a ], $
        linestyle=0, color=cgColor( 'BLU3' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_zy_a[ index_a ], z_diff2a[ index_a ], $
        linestyle=0, color=cgColor( 'BLU7' ), psym=16, thick=3.0, symsize=1.2 

    cgOPlot, hsc_zy_b[ index_b ], z_diff1b[ index_b ], $
        linestyle=0, color=cgColor( 'RED6' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_zy_b[ index_b ], z_diff3b[ index_b ], $
        linestyle=0, color=cgColor( 'GRN6' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_zy_b[ index_b ], z_diff4b[ index_b ], $
        linestyle=0, color=cgColor( 'BLK5' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_zy_b[ index_b ], z_diff5b[ index_b ], $
        linestyle=0, color=cgColor( 'TAN6' ), psym=9, thick=3.0, symsize=1.5
    cgOPlot, hsc_zy_b[ index_b ], z_diff6b[ index_b ], $
        linestyle=0, color=cgColor( 'Pink' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_zy_b[ index_b ], z_diff7b[ index_b ], $
        linestyle=0, color=cgColor( 'BLU3' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_zy_b[ index_b ], z_diff2b[ index_b ], $
        linestyle=0, color=cgColor( 'BLU7' ), psym=9, thick=3.0, symsize=1.5 

    cgOPlot, hsc_zy_c[ index_c ], z_diff1c[ index_c ], $
        linestyle=0, color=cgColor( 'RED6' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_zy_c[ index_c ], z_diff3c[ index_c ], $
        linestyle=0, color=cgColor( 'GRN6' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_zy_c[ index_c ], z_diff4c[ index_c ], $
        linestyle=0, color=cgColor( 'BLK5' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_zy_c[ index_c ], z_diff5c[ index_c ], $
        linestyle=0, color=cgColor( 'TAN6' ), psym=4, thick=3.0, symsize=1.5
    cgOPlot, hsc_zy_c[ index_c ], z_diff6c[ index_c ], $
        linestyle=0, color=cgColor( 'Pink' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_zy_c[ index_c ], z_diff7c[ index_c ], $
        linestyle=0, color=cgColor( 'BLU3' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_zy_c[ index_c ], z_diff2c[ index_c ], $
        linestyle=0, color=cgColor( 'BLU7' ), psym=4, thick=3.0, symsize=1.5 

    cgPlots, 0.18, 0.95, psym=16, symsize=2.5, color=cgColor( 'BLU7' ), /normal
    cgPlots, 0.20, 0.95, psym=9,  symsize=2.5, color=cgColor( 'BLU7' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.95, psym=4,  symsize=2.5, color=cgColor( 'BLU7' ), $
        /normal, thick=5
    cgText,  0.25, 0.94, 'SuprimeCam z', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.18, 0.91, psym=16, symsize=2.5, color=cgColor( 'RED6' ), /normal
    cgPlots, 0.20, 0.91, psym=9,  symsize=2.5, color=cgColor( 'RED6' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.91, psym=4,  symsize=2.5, color=cgColor( 'RED6' ), $
        /normal, thick=5
    cgText,  0.25, 0.90, 'SDSS z', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.18, 0.87, psym=16, symsize=2.5, color=cgColor( 'GRN6' ), /normal
    cgPlots, 0.20, 0.87, psym=9,  symsize=2.5, color=cgColor( 'GRN6' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.87, psym=4,  symsize=2.5, color=cgColor( 'GRN6' ), $
        /normal, thick=5
    cgText,  0.25, 0.86, 'MegaCam z', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.18, 0.83, psym=16, symsize=2.5, color=cgColor( 'BLK5' ), /normal
    cgPlots, 0.20, 0.83, psym=9,  symsize=2.5, color=cgColor( 'BLK5' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.83, psym=4,  symsize=2.5, color=cgColor( 'BLK5' ), $
        /normal, thick=5
    cgText,  0.25, 0.82, 'SuprimeCam z (old)', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.47, 0.95, psym=16, symsize=2.5, color=cgColor( 'TAN6' ), /normal
    cgPlots, 0.49, 0.95, psym=9,  symsize=2.5, color=cgColor( 'TAN6' ), $
        /normal, thick=5
    cgPlots, 0.51, 0.95, psym=4,  symsize=2.5, color=cgColor( 'TAN6' ), $
        /normal, thick=5
    cgText,  0.54, 0.94, 'DES z', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.47, 0.91, psym=16, symsize=2.5, color=cgColor( 'PINK' ), /normal
    cgPlots, 0.49, 0.91, psym=9,  symsize=2.5, color=cgColor( 'PINK' ), $
        /normal, thick=5
    cgPlots, 0.51, 0.91, psym=4,  symsize=2.5, color=cgColor( 'PINK' ), $
        /normal, thick=5
    cgText,  0.54, 0.90, 'ACS/HST F850LP', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.47, 0.87, psym=16, symsize=2.5, color=cgColor( 'BLU3' ), /normal
    cgPlots, 0.49, 0.87, psym=9,  symsize=2.5, color=cgColor( 'BLU3' ), $
        /normal, thick=5
    cgPlots, 0.51, 0.87, psym=4,  symsize=2.5, color=cgColor( 'BLU3' ), $
        /normal, thick=5
    cgText,  0.54, 0.86, 'WFC3/HST F850LP', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlot, hsc_zy_a[ index_a ], z_diff1a[ index_a ], $
        xstyle=1, ystyle=1, linestyle=0, color=cgColor( 'Dark Gray' ), $
        thick=3.0, charsize=3.5, charthick=8.0, xthick=10.0, ythick=10.0, $
        /noerase, position=position, yminor=5, $
        /nodata, xrange=xrange, yrange=yrange, $
        xticklen=0.025, yticklen=0.02

    device, /close 
    set_plot, mydevice 
    ;--------------------------------------------------------------------------;

    ;--------------------------------------------------------------------------;
    ; HSC-R band r_hsc 
    ; Compared with r_sdss, r_suprime, r_megacam, r_des 
    ;               f625w_acs
    ; Using r_hsc - i_hsc; i_hsc - z_hsc as x-axis 
    ;--------------------------------------------------------------------------;

    ;--------------------------------------------------------------------------;
    rplot_1 = 'hsc_r_filter_compare_1.eps'

    hsc_ri_a = ( mag2s.r_hsc - mag2s.i_hsc )
    hsc_iz_a = ( mag2s.i_hsc - mag2s.z_hsc )
    r_diff1a = ( mag2s.r_hsc - mag2s.r_sdss ) 
    r_diff2a = ( mag2s.r_hsc - mag2s.r_suprime ) 
    r_diff3a = ( mag2s.r_hsc - mag2s.r_megacam )
    r_diff5a = ( mag2s.r_hsc - mag2s.r_des )
    r_diff6a = ( mag2s.r_hsc - mag2s.f625w_acs )

    hsc_ri_b = ( mag1s.r_hsc - mag1s.i_hsc )
    hsc_iz_b = ( mag1s.i_hsc - mag1s.z_hsc )
    r_diff1b = ( mag1s.r_hsc - mag1s.r_sdss ) 
    r_diff2b = ( mag1s.r_hsc - mag1s.r_suprime ) 
    r_diff3b = ( mag1s.r_hsc - mag1s.r_megacam )
    r_diff5b = ( mag1s.r_hsc - mag1s.r_des )
    r_diff6b = ( mag1s.r_hsc - mag1s.f625w_acs )

    hsc_ri_c = ( mag3s.r_hsc - mag3s.i_hsc )
    hsc_iz_c = ( mag3s.i_hsc - mag3s.z_hsc )
    r_diff1c = ( mag3s.r_hsc - mag3s.r_sdss ) 
    r_diff2c = ( mag3s.r_hsc - mag3s.r_suprime ) 
    r_diff3c = ( mag3s.r_hsc - mag3s.r_megacam )
    r_diff5c = ( mag3s.r_hsc - mag3s.r_des )
    r_diff6c = ( mag3s.r_hsc - mag3s.f625w_acs )

    ;--------------------------------------------------------------------------;
    mydevice = !d.name 
    !p.font=1
    set_plot, 'ps' 
    device, filename=rplot_1, font_size=9.0, /encapsulated, $
        /color, set_font='TIMES-ROMAN', /bold, xsize=psxsize, ysize=psysize

    xmin = min( [ min( hsc_ri_a ), min( hsc_ri_b ), min( hsc_ri_c ) ] )
    xmax = max( [ max( hsc_ri_a ), max( hsc_ri_b ), max( hsc_ri_c ) ] )
    xrange = [ xmin - 0.01, xmax + 0.04 ]
    ;xrange = [ min( alog10( mag2s.agegyr ) ) - 0.05, $
    ;           max( alog10( mag2s.agegyr ) ) + 0.10 ]

    ymin = min( [ min( r_diff1a ), min( r_diff2a ), min( r_diff3a ), $
                  min( r_diff5a ), min( r_diff6a ) ] )
    ymax = max( [ max( r_diff1a ), max( r_diff2a ), max( r_diff3a ), $
                  max( r_diff5a ), max( r_diff6a ) ] )
    yrange = [ ymin - 0.005, ymax + 0.01 ]

    index_a = sort( hsc_ri_a )
    index_b = sort( hsc_ri_b )
    index_c = sort( hsc_ri_c )

    cgPlot, hsc_ri_a[ index_a ], r_diff1a[ index_a ], $
        xstyle=1, ystyle=1, linestyle=0, color=cgColor( 'Dark Gray' ), $
        thick=3.0, charsize=3.5, charthick=8.0, xthick=10.0, ythick=10.0, $
        /noerase, position=position, yminor=5, $
        /nodata, xrange=xrange, yrange=yrange, $
        xticklen=0.025, yticklen=0.02, $
        xtitle=TeXToIDL( 'r_{HSC} - i_{HSC} ') + '(mag)', $
        ytitle=TeXToIDL( 'r_{HSC} - r_{Other} ') + '(mag)'

    cgOPlot, !X.Crange, [ 0.0, 0.0 ], linestyle=2, thick=5.0, $
        color=cgColor( 'Dark Gray' )

    ;cgOPlot, mag2s.agegyr, r_diff1a, $
    cgOPlot, hsc_ri_a[ index_a ], r_diff1a[ index_a ], $
        linestyle=0, color=cgColor( 'RED6' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_ri_a[ index_a ], r_diff3a[ index_a ], $
        linestyle=0, color=cgColor( 'GRN6' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_ri_a[ index_a ], r_diff5a[ index_a ], $
        linestyle=0, color=cgColor( 'TAN6' ), psym=16, thick=3.0, symsize=1.2
    cgOPlot, hsc_ri_a[ index_a ], r_diff6a[ index_a ], $
        linestyle=0, color=cgColor( 'Pink' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_ri_a[ index_a ], r_diff2a[ index_a ], $
        linestyle=0, color=cgColor( 'BLU7' ), psym=16, thick=3.0, symsize=1.2 

    cgOPlot, hsc_ri_b[ index_b ], r_diff1b[ index_b ], $
        linestyle=0, color=cgColor( 'RED6' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_ri_b[ index_b ], r_diff3b[ index_b ], $
        linestyle=0, color=cgColor( 'GRN6' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_ri_b[ index_b ], r_diff5b[ index_b ], $
        linestyle=0, color=cgColor( 'TAN6' ), psym=9, thick=3.0, symsize=1.5
    cgOPlot, hsc_ri_b[ index_b ], r_diff6b[ index_b ], $
        linestyle=0, color=cgColor( 'Pink' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_ri_b[ index_b ], r_diff2b[ index_b ], $
        linestyle=0, color=cgColor( 'BLU7' ), psym=9, thick=3.0, symsize=1.5 

    cgOPlot, hsc_ri_c[ index_c ], r_diff1c[ index_c ], $
        linestyle=0, color=cgColor( 'RED6' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_ri_c[ index_c ], r_diff3c[ index_c ], $
        linestyle=0, color=cgColor( 'GRN6' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_ri_c[ index_c ], r_diff5c[ index_c ], $
        linestyle=0, color=cgColor( 'TAN6' ), psym=4, thick=3.0, symsize=1.5
    cgOPlot, hsc_ri_c[ index_c ], r_diff6c[ index_c ], $
        linestyle=0, color=cgColor( 'Pink' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_ri_c[ index_c ], r_diff2c[ index_c ], $
        linestyle=0, color=cgColor( 'BLU7' ), psym=4, thick=3.0, symsize=1.5 

    cgPlots, 0.18, 0.95, psym=16, symsize=2.5, color=cgColor( 'BLU7' ), /normal
    cgPlots, 0.20, 0.95, psym=9,  symsize=2.5, color=cgColor( 'BLU7' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.95, psym=4,  symsize=2.5, color=cgColor( 'BLU7' ), $
        /normal, thick=5
    cgText,  0.25, 0.94, 'SuprimeCam r', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.18, 0.91, psym=16, symsize=2.5, color=cgColor( 'RED6' ), /normal
    cgPlots, 0.20, 0.91, psym=9,  symsize=2.5, color=cgColor( 'RED6' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.91, psym=4,  symsize=2.5, color=cgColor( 'RED6' ), $
        /normal, thick=5
    cgText,  0.25, 0.90, 'SDSS r', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.18, 0.87, psym=16, symsize=2.5, color=cgColor( 'GRN6' ), /normal
    cgPlots, 0.20, 0.87, psym=9,  symsize=2.5, color=cgColor( 'GRN6' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.87, psym=4,  symsize=2.5, color=cgColor( 'GRN6' ), $
        /normal, thick=5
    cgText,  0.25, 0.86, 'MegaCam r', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.47, 0.95, psym=16, symsize=2.5, color=cgColor( 'TAN6' ), /normal
    cgPlots, 0.49, 0.95, psym=9,  symsize=2.5, color=cgColor( 'TAN6' ), $
        /normal, thick=5
    cgPlots, 0.51, 0.95, psym=4,  symsize=2.5, color=cgColor( 'TAN6' ), $
        /normal, thick=5
    cgText,  0.54, 0.94, 'DES r', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.47, 0.91, psym=16, symsize=2.5, color=cgColor( 'PINK' ), /normal
    cgPlots, 0.49, 0.91, psym=9,  symsize=2.5, color=cgColor( 'PINK' ), $
        /normal, thick=5
    cgPlots, 0.51, 0.91, psym=4,  symsize=2.5, color=cgColor( 'PINK' ), $
        /normal, thick=5
    cgText,  0.54, 0.90, 'ACS/HST F625W', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlot, hsc_ri_a[ index_a ], r_diff1a[ index_a ], $
        xstyle=1, ystyle=1, linestyle=0, color=cgColor( 'Dark Gray' ), $
        thick=3.0, charsize=3.5, charthick=8.0, xthick=10.0, ythick=10.0, $
        /noerase, position=position, yminor=5, $
        /nodata, xrange=xrange, yrange=yrange, $
        xticklen=0.025, yticklen=0.02

    device, /close 
    set_plot, mydevice 
    ;--------------------------------------------------------------------------;

    ;--------------------------------------------------------------------------;
    ; HSC-G band g_hsc 
    ; Compared with g_sdss, g_suprime, g_megacam, g_des 
    ;               f475w_acs, f475w_wfc3
    ; Using g_hsc - r_hsc; r_hsc - i_hsc as x-axis 
    ;--------------------------------------------------------------------------;

    ;--------------------------------------------------------------------------;
    gplot_1 = 'hsc_g_filter_compare_1.eps'

    hsc_gr_a = ( mag2s.g_hsc - mag2s.r_hsc )
    hsc_ig_a = ( mag2s.i_hsc - mag2s.g_hsc )
    g_diff1a = ( mag2s.g_hsc - mag2s.g_sdss ) 
    g_diff2a = ( mag2s.g_hsc - mag2s.g_suprime ) 
    g_diff3a = ( mag2s.g_hsc - mag2s.g_megacam )
    g_diff5a = ( mag2s.g_hsc - mag2s.g_des )
    g_diff6a = ( mag2s.g_hsc - mag2s.f475w_acs )
    g_diff7a = ( mag2s.g_hsc - mag2s.f475w_wfc3 )

    hsc_gr_b = ( mag1s.g_hsc - mag1s.r_hsc )
    hsc_ig_b = ( mag1s.i_hsc - mag1s.g_hsc )
    g_diff1b = ( mag1s.g_hsc - mag1s.g_sdss ) 
    g_diff2b = ( mag1s.g_hsc - mag1s.g_suprime ) 
    g_diff3b = ( mag1s.g_hsc - mag1s.g_megacam )
    g_diff5b = ( mag1s.g_hsc - mag1s.g_des )
    g_diff6b = ( mag1s.g_hsc - mag1s.f475w_acs )
    g_diff7b = ( mag1s.g_hsc - mag1s.f475w_wfc3 )

    hsc_gr_c = ( mag3s.g_hsc - mag3s.r_hsc )
    hsc_ig_c = ( mag3s.i_hsc - mag3s.g_hsc )
    g_diff1c = ( mag3s.g_hsc - mag3s.g_sdss ) 
    g_diff2c = ( mag3s.g_hsc - mag3s.g_suprime ) 
    g_diff3c = ( mag3s.g_hsc - mag3s.g_megacam )
    g_diff5c = ( mag3s.g_hsc - mag3s.g_des )
    g_diff6c = ( mag3s.g_hsc - mag3s.f475w_acs )
    g_diff7c = ( mag3s.g_hsc - mag3s.f475w_wfc3 )

    ;--------------------------------------------------------------------------;
    mydevice = !d.name 
    !p.font=1
    set_plot, 'ps' 
    device, filename=gplot_1, font_size=9.0, /encapsulated, $
        /color, set_font='TIMES-ROMAN', /bold, xsize=psxsize, ysize=psysize

    xmin = min( [ min( hsc_gr_a ), min( hsc_gr_b ), min( hsc_gr_c ) ] )
    xmax = max( [ max( hsc_gr_a ), max( hsc_gr_b ), max( hsc_gr_c ) ] )
    xrange = [ xmin - 0.01, xmax + 0.04 ]
    ;xrange = [ min( alog10( mag2s.agegyr ) ) - 0.05, $
    ;           max( alog10( mag2s.agegyr ) ) + 0.10 ]

    ymin = min( [ min( g_diff1a ), min( g_diff2a ), min( g_diff3a ), $
                  min( g_diff5a ), min( g_diff6a ), min( g_diff7a ) ] )
    ymax = max( [ max( g_diff1a ), max( g_diff2a ), max( g_diff3a ), $
                  max( g_diff7a ), max( g_diff5a ), max( g_diff6a ) ] )
    yrange = [ ymin - 0.005, ymax + 0.01 ]

    index_a = sort( hsc_gr_a )
    index_b = sort( hsc_gr_b )
    index_c = sort( hsc_gr_c )

    cgPlot, hsc_gr_a[ index_a ], g_diff1a[ index_a ], $
        xstyle=1, ystyle=1, linestyle=0, color=cgColor( 'Dark Gray' ), $
        thick=3.0, charsize=3.5, charthick=8.0, xthick=10.0, ythick=10.0, $
        /noerase, position=position, yminor=5, $
        /nodata, xrange=xrange, yrange=yrange, $
        xticklen=0.025, yticklen=0.02, $
        xtitle=TeXToIDL( 'g_{HSC} - r_{HSC} ') + '(mag)', $
        ytitle=TeXToIDL( 'g_{HSC} - g_{Other} ') + '(mag)'

    cgOPlot, !X.Crange, [ 0.0, 0.0 ], linestyle=2, thick=5.0, $
        color=cgColor( 'Dark Gray' )

    ;cgOPlot, mag2s.agegyr, g_diff1a, $
    cgOPlot, hsc_gr_a[ index_a ], g_diff1a[ index_a ], $
        linestyle=0, color=cgColor( 'RED6' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_gr_a[ index_a ], g_diff3a[ index_a ], $
        linestyle=0, color=cgColor( 'GRN6' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_gr_a[ index_a ], g_diff5a[ index_a ], $
        linestyle=0, color=cgColor( 'TAN6' ), psym=16, thick=3.0, symsize=1.2
    cgOPlot, hsc_gr_a[ index_a ], g_diff6a[ index_a ], $
        linestyle=0, color=cgColor( 'Pink' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_gr_a[ index_a ], g_diff7a[ index_a ], $
        linestyle=0, color=cgColor( 'BLU3' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_gr_a[ index_a ], g_diff2a[ index_a ], $
        linestyle=0, color=cgColor( 'BLU7' ), psym=16, thick=3.0, symsize=1.2 

    cgOPlot, hsc_gr_b[ index_b ], g_diff1b[ index_b ], $
        linestyle=0, color=cgColor( 'RED6' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_gr_b[ index_b ], g_diff3b[ index_b ], $
        linestyle=0, color=cgColor( 'GRN6' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_gr_b[ index_b ], g_diff5b[ index_b ], $
        linestyle=0, color=cgColor( 'TAN6' ), psym=9, thick=3.0, symsize=1.5
    cgOPlot, hsc_gr_b[ index_b ], g_diff6b[ index_b ], $
        linestyle=0, color=cgColor( 'Pink' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_gr_b[ index_b ], g_diff7b[ index_b ], $
        linestyle=0, color=cgColor( 'BLU3' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_gr_b[ index_b ], g_diff2b[ index_b ], $
        linestyle=0, color=cgColor( 'BLU7' ), psym=9, thick=3.0, symsize=1.5 

    cgOPlot, hsc_gr_c[ index_c ], g_diff1c[ index_c ], $
        linestyle=0, color=cgColor( 'RED6' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_gr_c[ index_c ], g_diff3c[ index_c ], $
        linestyle=0, color=cgColor( 'GRN6' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_gr_c[ index_c ], g_diff5c[ index_c ], $
        linestyle=0, color=cgColor( 'TAN6' ), psym=4, thick=3.0, symsize=1.5
    cgOPlot, hsc_gr_c[ index_c ], g_diff6c[ index_c ], $
        linestyle=0, color=cgColor( 'Pink' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_gr_c[ index_c ], g_diff7c[ index_c ], $
        linestyle=0, color=cgColor( 'BLU3' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_gr_c[ index_c ], g_diff2c[ index_c ], $
        linestyle=0, color=cgColor( 'BLU7' ), psym=4, thick=3.0, symsize=1.5 

    cgPlots, 0.18, 0.95, psym=16, symsize=2.5, color=cgColor( 'BLU7' ), /normal
    cgPlots, 0.20, 0.95, psym=9,  symsize=2.5, color=cgColor( 'BLU7' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.95, psym=4,  symsize=2.5, color=cgColor( 'BLU7' ), $
        /normal, thick=5
    cgText,  0.25, 0.94, 'SuprimeCam g', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.18, 0.91, psym=16, symsize=2.5, color=cgColor( 'RED6' ), /normal
    cgPlots, 0.20, 0.91, psym=9,  symsize=2.5, color=cgColor( 'RED6' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.91, psym=4,  symsize=2.5, color=cgColor( 'RED6' ), $
        /normal, thick=5
    cgText,  0.25, 0.90, 'SDSS g', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.18, 0.87, psym=16, symsize=2.5, color=cgColor( 'GRN6' ), /normal
    cgPlots, 0.20, 0.87, psym=9,  symsize=2.5, color=cgColor( 'GRN6' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.87, psym=4,  symsize=2.5, color=cgColor( 'GRN6' ), $
        /normal, thick=5
    cgText,  0.25, 0.86, 'MegaCam g', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.47, 0.95, psym=16, symsize=2.5, color=cgColor( 'TAN6' ), /normal
    cgPlots, 0.49, 0.95, psym=9,  symsize=2.5, color=cgColor( 'TAN6' ), $
        /normal, thick=5
    cgPlots, 0.51, 0.95, psym=4,  symsize=2.5, color=cgColor( 'TAN6' ), $
        /normal, thick=5
    cgText,  0.54, 0.94, 'DES g', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.47, 0.91, psym=16, symsize=2.5, color=cgColor( 'PINK' ), /normal
    cgPlots, 0.49, 0.91, psym=9,  symsize=2.5, color=cgColor( 'PINK' ), $
        /normal, thick=5
    cgPlots, 0.51, 0.91, psym=4,  symsize=2.5, color=cgColor( 'PINK' ), $
        /normal, thick=5
    cgText,  0.54, 0.90, 'ACS/HST F475W', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.47, 0.87, psym=16, symsize=2.5, color=cgColor( 'BLU3' ), /normal
    cgPlots, 0.49, 0.87, psym=9,  symsize=2.5, color=cgColor( 'BLU3' ), $
        /normal, thick=5
    cgPlots, 0.51, 0.87, psym=4,  symsize=2.5, color=cgColor( 'BLU3' ), $
        /normal, thick=5
    cgText,  0.54, 0.86, 'WFC3/HST F475W', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlot, hsc_gr_a[ index_a ], g_diff1a[ index_a ], $
        xstyle=1, ystyle=1, linestyle=0, color=cgColor( 'Dark Gray' ), $
        thick=3.0, charsize=3.5, charthick=8.0, xthick=10.0, ythick=10.0, $
        /noerase, position=position, yminor=5, $
        /nodata, xrange=xrange, yrange=yrange, $
        xticklen=0.025, yticklen=0.02

    device, /close 
    set_plot, mydevice 
    ;--------------------------------------------------------------------------;

    ;--------------------------------------------------------------------------;
    ; HSC-Y band y_hsc 
    ; Compared with y_des, y_vircam, y_wfcam, 
    ;               f105w_wfc3
    ; Using z_hsc - y_hsc as x-axis 
    ;--------------------------------------------------------------------------;

    ;--------------------------------------------------------------------------;
    yplot_1 = 'hsc_y_filter_compare_1.eps'

    hsc_zy_a = ( mag2s.z_hsc - mag2s.y_hsc )
    y_diff1a = ( mag2s.y_hsc - mag2s.y_vircam ) 
    y_diff3a = ( mag2s.y_hsc - mag2s.y_wfcam )
    y_diff5a = ( mag2s.y_hsc - mag2s.y_des )
    y_diff7a = ( mag2s.y_hsc - mag2s.f105w_wfc3 )

    hsc_zy_b = ( mag1s.z_hsc - mag1s.y_hsc )
    y_diff1b = ( mag1s.y_hsc - mag1s.y_vircam ) 
    y_diff3b = ( mag1s.y_hsc - mag1s.y_wfcam )
    y_diff5b = ( mag1s.y_hsc - mag1s.y_des )
    y_diff7b = ( mag1s.y_hsc - mag1s.f105w_wfc3 )

    hsc_zy_c = ( mag3s.z_hsc - mag3s.y_hsc )
    y_diff1c = ( mag3s.y_hsc - mag3s.y_vircam ) 
    y_diff3c = ( mag3s.y_hsc - mag3s.y_wfcam )
    y_diff5c = ( mag3s.y_hsc - mag3s.y_des )
    y_diff7c = ( mag3s.y_hsc - mag3s.f105w_wfc3 )

    ;--------------------------------------------------------------------------;
    mydevice = !d.name 
    !p.font=1
    set_plot, 'ps' 
    device, filename=yplot_1, font_size=9.0, /encapsulated, $
        /color, set_font='TIMES-ROMAN', /bold, xsize=psxsize, ysize=psysize

    xmin = min( [ min( hsc_zy_a ), min( hsc_zy_b ), min( hsc_zy_c ) ] )
    xmax = max( [ max( hsc_zy_a ), max( hsc_zy_b ), max( hsc_zy_c ) ] )
    xrange = [ xmin - 0.01, xmax + 0.04 ]
    ;xrange = [ min( alog10( mag2s.agegyr ) ) - 0.05, $
    ;           max( alog10( mag2s.agegyr ) ) + 0.10 ]

    ymin = min( [ min( y_diff1a ), min( y_diff3a ), min( y_diff5a ), $
                  min( y_diff7a ) ] )
    ymax = max( [ max( y_diff1a ), max( y_diff3a ), max( y_diff5a ), $
                  max( y_diff7a ) ] )
    yrange = [ ymin - 0.005, ymax + 0.025 ]

    index_a = sort( hsc_zy_a )
    index_b = sort( hsc_zy_b )
    index_c = sort( hsc_zy_c )

    cgPlot, hsc_zy_a[ index_a ], y_diff1a[ index_a ], $
        xstyle=1, ystyle=1, linestyle=0, color=cgColor( 'Dark Gray' ), $
        thick=3.0, charsize=3.5, charthick=8.0, xthick=10.0, ythick=10.0, $
        /noerase, position=position, yminor=5, $
        /nodata, xrange=xrange, yrange=yrange, $
        xticklen=0.025, yticklen=0.02, $
        xtitle=TeXToIDL( 'z_{HSC} - Y_{HSC} ') + '(mag)', $
        ytitle=TeXToIDL( 'Y_{HSC} - Y_{Other} ') + '(mag)'

    cgOPlot, !X.Crange, [ 0.0, 0.0 ], linestyle=2, thick=5.0, $
        color=cgColor( 'Dark Gray' )

    ;cgOPlot, mag2s.agegyr, y_diff1a, $
    cgOPlot, hsc_zy_a[ index_a ], y_diff1a[ index_a ], $
        linestyle=0, color=cgColor( 'RED6' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_zy_a[ index_a ], y_diff3a[ index_a ], $
        linestyle=0, color=cgColor( 'GRN6' ), psym=16, thick=3.0, symsize=1.2 
    cgOPlot, hsc_zy_a[ index_a ], y_diff5a[ index_a ], $
        linestyle=0, color=cgColor( 'TAN6' ), psym=16, thick=3.0, symsize=1.2
    cgOPlot, hsc_zy_a[ index_a ], y_diff7a[ index_a ], $
        linestyle=0, color=cgColor( 'BLU3' ), psym=16, thick=3.0, symsize=1.2 

    cgOPlot, hsc_zy_b[ index_b ], y_diff1b[ index_b ], $
        linestyle=0, color=cgColor( 'RED6' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_zy_b[ index_b ], y_diff3b[ index_b ], $
        linestyle=0, color=cgColor( 'GRN6' ), psym=9, thick=3.0, symsize=1.5 
    cgOPlot, hsc_zy_b[ index_b ], y_diff5b[ index_b ], $
        linestyle=0, color=cgColor( 'TAN6' ), psym=9, thick=3.0, symsize=1.5
    cgOPlot, hsc_zy_b[ index_b ], y_diff7b[ index_b ], $
        linestyle=0, color=cgColor( 'BLU3' ), psym=9, thick=3.0, symsize=1.5 

    cgOPlot, hsc_zy_c[ index_c ], y_diff1c[ index_c ], $
        linestyle=0, color=cgColor( 'RED6' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_zy_c[ index_c ], y_diff3c[ index_c ], $
        linestyle=0, color=cgColor( 'GRN6' ), psym=4, thick=3.0, symsize=1.5 
    cgOPlot, hsc_zy_c[ index_c ], y_diff5c[ index_c ], $
        linestyle=0, color=cgColor( 'TAN6' ), psym=4, thick=3.0, symsize=1.5
    cgOPlot, hsc_zy_c[ index_c ], y_diff7c[ index_c ], $
        linestyle=0, color=cgColor( 'BLU3' ), psym=4, thick=3.0, symsize=1.5 

    cgPlots, 0.18, 0.95, psym=16, symsize=2.5, color=cgColor( 'BLU7' ), /normal
    cgPlots, 0.20, 0.95, psym=9,  symsize=2.5, color=cgColor( 'BLU7' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.95, psym=4,  symsize=2.5, color=cgColor( 'BLU7' ), $
        /normal, thick=5
    cgText,  0.25, 0.94, 'VIRCAM Y', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.18, 0.91, psym=16, symsize=2.5, color=cgColor( 'RED6' ), /normal
    cgPlots, 0.20, 0.91, psym=9,  symsize=2.5, color=cgColor( 'RED6' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.91, psym=4,  symsize=2.5, color=cgColor( 'RED6' ), $
        /normal, thick=5
    cgText,  0.25, 0.90, 'WFCAM Y', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.18, 0.87, psym=16, symsize=2.5, color=cgColor( 'TAN6' ), /normal
    cgPlots, 0.20, 0.87, psym=9,  symsize=2.5, color=cgColor( 'TAN6' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.87, psym=4,  symsize=2.5, color=cgColor( 'TAN6' ), $
        /normal, thick=5
    cgText,  0.25, 0.86, 'DES Y', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlots, 0.18, 0.83, psym=16, symsize=2.5, color=cgColor( 'BLU3' ), /normal
    cgPlots, 0.20, 0.83, psym=9,  symsize=2.5, color=cgColor( 'BLU3' ), $
        /normal, thick=5
    cgPlots, 0.22, 0.83, psym=4,  symsize=2.5, color=cgColor( 'BLU3' ), $
        /normal, thick=5
    cgText,  0.25, 0.82, 'WFC3/HST F105W', charthick=3.0, charsize=3.0, $
        color=cgColor( 'Black' ), alignment=0, /normal

    cgPlot, hsc_zy_a[ index_a ], y_diff1a[ index_a ], $
        xstyle=1, ystyle=1, linestyle=0, color=cgColor( 'Dark Gray' ), $
        thick=3.0, charsize=3.5, charthick=8.0, xthick=10.0, ythick=10.0, $
        /noerase, position=position, yminor=5, $
        /nodata, xrange=xrange, yrange=yrange, $
        xticklen=0.025, yticklen=0.02

    device, /close 
    set_plot, mydevice 
    ;--------------------------------------------------------------------------;

end
