function dist2, vec1, vec2
    return, ( vec1[0] - vec2[0] )^2.0 + ( vec1[1] - vec2[1] )^2.0 
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function hs_circle_rectangle_overlap, center, radius, vec1, vec2, vec3, vec4 

    if ( n_elements( center ) NE 2 ) then begin 
        print, 'The center of the Circle should be described as [x,y]'
        return, -1 
    endif 

    if ( ( n_elements( vec1 ) NE 2 ) OR ( n_elements( vec2 ) NE 2 ) OR $ 
         ( n_elements( vec3 ) NE 2 ) OR ( n_elements( vec4 ) NE 2 ) ) then begin 
        print, 'The vectors for the four corners of the rectangle should be '
        print, '    described as [x,y] !!'
        return, -1 
    endif 

    if ( radius LE 0.0 ) then begin 
        print, 'The radius should be larger than 0!'
        return, -1 
    endif

    ;; 1. If the center of the circle is inside the rectangle 
    rect_x = [ vec1[0], vec2[0], vec3[0], vec4[0] ] 
    rect_y = [ vec1[1], vec2[1], vec3[1], vec4[1] ] 
    center_circle_inside = inside( center[0], center[1], rect_x, rect_y )

    ;; 2. If the rectangle is inside the circle 
    rect_center = [ ( vec1[0] + vec3[0] ) / 2.0, ( vec1[1] + vec3[1] ) / 2.0 ]
    dist_center = dist2( rect_center, center ) 
    if ( dist_center LE radius^2.0 ) then begin 
        center_rect_inside = 1 
    endif else begin 
        center_rect_inside = 0 
    endelse

    ;; 3. If any of the four edges of the rectangle is overlaped with the circle
    edge1_overlap = hs_circle_seg_intersect( center, radius, vec1, vec2 )
    edge2_overlap = hs_circle_seg_intersect( center, radius, vec2, vec3 )
    edge3_overlap = hs_circle_seg_intersect( center, radius, vec3, vec4 )
    edge4_overlap = hs_circle_seg_intersect( center, radius, vec4, vec1 )

    if ( center_circle_inside OR center_rect_inside OR $
         edge1_overlap OR edge2_overlap OR edge3_overlap OR $
         edge4_overlap ) then begin 
        circle_rect_overlap = 1 
    endif else begin 
        circle_rect_overlap = 0 
    endelse

    return, circle_rect_overlap

end

pro test 

    vec1 = [ -3, -2 ] 
    vec2 = [  3, -2 ]
    vec3 = [  3,  2 ]
    vec4 = [ -3,  2 ]
    
    center = [ 4.0, 0.0 ]
    radius = 1.01 
    if ( hs_circle_rectangle_overlap( center, radius, vec1, vec2, vec3, $
        vec4 ) ) then begin 
        print, 'Overlaped !'
    endif else begin 
        print, 'Not overlaped !'
    endelse

end
