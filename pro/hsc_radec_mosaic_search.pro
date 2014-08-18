function hsc_radec_mosaic_search, ra, dec, radius, mosaic_catalog 

    if NOT file_test( mosaic_catalog ) then begin 
        print, 'Can not find the catalog for mosaic information, Check !'
        return, -1 
    endif else begin 
        data = mrdfits( mosaic_catalog, 1, head, /silent ) 
        n_mosaic = n_elements( data.skytile_id )
    endelse

    if ( ( n_elements( ra ) NE 1 ) OR ( n_elements( dec ) NE 1 ) ) then begin 
        print, 'One RA-DEC at a time please !' 
        return, -1
    endif 

    if ( radius LE 0 ) then begin 
        print, 'The search radius should larger than 0!'
        return, -1 
    endif 

    overlap   = intarr( n_mosaic ) 
    overcount = 0

    center = [ ra, dec ]

    for ii = 0, ( n_mosaic - 1 ), 1 do begin 

        vec1 = [ data[ii].llcra, data[ii].llcdecl ]
        vec2 = [ data[ii].lrcra, data[ii].lrcdecl ]
        vec3 = [ data[ii].urcra, data[ii].urcdecl ]
        vec4 = [ data[ii].ulcra, data[ii].ulcdecl ]

        if ( hs_circle_rectangle_overlap( center, radius, vec1, vec2, vec3, $
            vec4 ) ) then begin 
            overcount += 1
            overlap[ii] = 1
            print, '# Overlap ' + string( overcount ) + ' : ' + $
                string( data[ii].tract ) + '  ' + string( data[ii].patch )
        endif else begin 
            overlap[ii] = 0 
        endelse

    endfor 
    
    if ( overcount EQ 0 ) then begin 
        print, 'No overlapped mosaic is found !!'
    endif else begin 
        if ( overcount EQ 1 ) then begin 
            print, '# 1 overlapped mosaic is found !'
        endif else begin 
            print, '# ' + string( overcount ) + ' overlapped mosaics are found !'
        endelse
    endelse

    if ( overcount GT 0 ) then begin 
        mosaic_overlap = data[ where( overlap EQ 1 ) ]
    endif else begin 
        mosaic_overlap = -1 
    endelse

    return, mosaic_overlap 

end


pro test 

    catalog = '/home/hs/work/hsc_test/catalog/mosaic/hsc_udeep_i_20140523a_mosaic.fits' 

    ra  = 150.23983 
    dec =   2.56283
    radius = 0.06

    overlap = hsc_radec_mosaic_search( ra, dec, radius, catalog ) 

end
