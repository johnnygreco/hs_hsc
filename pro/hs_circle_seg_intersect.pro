function hs_circle_seg_intersect, center, radius, vec1, vec2 

    shortest = hs_dist_point_lineseg( center, vec1, vec2 ) 

    if ( shortest LE radius ) then begin 
        return, 1 
    endif else begin 
        return, 0 
    endelse 

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro test 

    center = [ 5.0, 5.0 ] 
    radius = 4.01 
    vec1   = [ -6.0, 1.0 ] 
    vec2   = [  0.0, 1.0 ] 
    if ( hs_circle_seg_intersect( center, radius, vec1, vec2 ) ) then begin 
        print, 'Intersect !' 
    endif else begin 
        print, 'Not Intersect !' 
    endelse

end
