pro udeep_i_mosaic 

    cat = mrdfits( 'udeep_i_mosaic.fits', 1, head ) 
    
    num = n_elements( cat.tract ) 

    ;; Size of the EPS output 
    psxsize=50 
    psysize=50
    ;; Location of each sub-figures 
    position = [ 0.12, 0.12, 0.99, 0.99 ] 

    psname = 'udeep_i_mosaic.eps'

    xrange = [ 148.9, 151.2 ] 
    yrange = [   1.3,   3.1 ] 

    transparent = 40 
    alpha = 1 - ( transparent / 100.0 )

    mydevice = !d.name 
    !p.font=1
    set_plot, 'ps' 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    device, filename=psname, font_size=12.0, /encapsulated, $
        /color, set_font='TIMES-ROMAN', /bold, $
        xsize=psxsize, ysize=psysize

    cgPlot, cat.llcra, cat.llcdecl, xstyle=1, ystyle=1, $
        xrange=xrange, yrange=yrange, $ 
        xthick=12.0, ythick=12.0, charsize=5.0, charthick=12.0, $
        xtitle='RA (J2000)', ytitle='Dec (J2000)', $
        color=cgColor( 'Black' ), /noerase, /nodata, $
        position=position, xticklen=0.01, yticklen=0.01 

    for ii = 0, ( num - 1 ), 1 do begin 

        patch_x = [ cat[ii].llcra, cat[ii].lrcra, cat[ii].urcra, $
            cat[ii].ulcra, cat[ii].llcra ]
        patch_y = [ cat[ii].llcdecl, cat[ii].lrcdecl, cat[ii].urcdecl, $
            cat[ii].ulcdecl, cat[ii].llcdecl ]

        if ( ( ii MOD 2 ) EQ 0 ) then begin 
            patch_color = 'Red' 
        endif else begin 
            patch_color = 'Blue' 
        endelse

        cgPlots, patch_x, patch_y, color=cgColor( patch_color ), thick=3.0 

    endfor 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    device, /close 
    set_plot, mydevice

end
