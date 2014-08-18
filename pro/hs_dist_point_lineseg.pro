function dist2, vec1, vec2
    return, ( vec1[0] - vec2[0] )^2.0 + ( vec1[1] - vec2[1] )^2.0 
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function hs_dist_point_lineseg, point, vec1, vec2 

    length = dist2( vec1, vec2 ) 

    if ( length EQ 0 ) then begin 
        pl_dist = dist2( point, vec1 ) 
    endif else begin 
        tt = ( ( point[0] - vec1[0] ) * ( vec2[0] - vec1[0] ) + $ 
               ( point[1] - vec1[1] ) * ( vec2[1] - vec1[1] ) ) / length 
        if ( tt LE 0 ) then begin 
            pl_dist = dist2( point, vec1 ) 
        endif else begin 
            if ( tt GT 1.0 ) then begin 
                pl_dist = dist2( point, vec2 ) 
            endif else begin 
                xx = ( vec1[0] + tt * ( vec2[0] - vec1[0] ) ) 
                yy = ( vec1[1] + tt * ( vec2[1] - vec1[1] ) ) 
                pp = [ xx, yy ]
                pl_dist = dist2( point, pp )
            endelse 
        endelse
    endelse

    return, sqrt( float( pl_dist ) ) 

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro test 

    point = [ 142.5, -6.0 ]
    vec1  = [ 60.0,   2.5 ] 
    vec2  = [ 150.0, -2.5 ] 

    print, hs_dist_point_lineseg( point, vec1, vec2 ) 

end 
