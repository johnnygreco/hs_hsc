function hsc_radec_frame_search, ra, dec, radius, frame_catalog 

    if NOT file_test( frame_catalog ) then begin 
        print, 'Can not find the catalog for frame information, Check !'
        return, -1 
    endif else begin 
        data = mrdfits( frame_catalog, 1, head, /silent ) 
        n_frame = n_elements( data.frame_id )
    endelse

    if ( ( n_elements( ra ) NE 1 ) OR ( n_elements( dec ) NE 1 ) ) then begin 
        print, 'One RA-DEC at a time please !' 
        return, -1
    endif 

    if ( radius LE 0 ) then begin 
        print, 'The search radius should larger than 0!'
        return, -1 
    endif 

    overlap   = intarr( n_frame ) 
    overcount = 0

    center = [ ra, dec ]
    purpos = strcompress( strlowcase( data.purpos ), /remove_all )

    for ii = 0, ( n_frame - 1 ), 1 do begin 

        vec1 = [ data[ii].llcra, data[ii].llcdecl ]
        vec2 = [ data[ii].lrcra, data[ii].lrcdecl ]
        vec3 = [ data[ii].urcra, data[ii].urcdecl ]
        vec4 = [ data[ii].ulcra, data[ii].ulcdecl ]

        if ( hs_circle_rectangle_overlap( center, radius, vec1, vec2, vec3, $
            vec4 ) AND ( purpos[ii] EQ 'object' ) ) then begin 
            overcount += 1
            overlap[ii] = 1
            print, '# Overlap ' + string( overcount ) + ' : ' + $
                string( data[ii].visit ) + '  ' + string( data[ii].ccd )
        endif else begin 
            overlap[ii] = 0 
        endelse

    endfor 
    
    if ( overcount EQ 0 ) then begin 
        print, 'No overlapped frame is found !!'
    endif else begin 
        if ( overcount EQ 1 ) then begin 
            print, '# 1 overlapped frame is found !'
        endif else begin 
            print, '# ' + string( overcount ) + ' overlapped frames are found !'
        endelse
    endelse

    if ( overcount GT 0 ) then begin 
        frame_overlap = data[ where( overlap EQ 1 ) ]
    endif else begin 
        frame_overlap = -1 
    endelse

    return, frame_overlap 

end


pro test 

    catalog = '/home/hs/work/hsc_test/catalog/frame/hsc_udeep_i_20140523a_frameinfo.fits' 

    ra  = 150.23983 
    dec =   2.56283
    radius = 0.06

    overlap = hsc_radec_frame_search( ra, dec, radius, catalog ) 

end
