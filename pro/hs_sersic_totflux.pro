;+
; NAME:
;             HS_SERSIC_TOTFLUX
;
; PURPOSE:
;             Get the total flux of a surface brightness distribution that is 
;             modeled by the Sersic function
;
; USAGE:
; 
; OUTPUT: 
;
; AUTHOR:
;             Song Huang
;
; HISTORY:
;             Song Huang, 2014/08/14 - First version 
;-
; CATEGORY:    HS_GALAXY
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; The first two functions are based on the Get_SersicB() from John Moustakas
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sersicb_func, bb
common sersicb, nn
return, gamma(2.0*nn)-2D*igamma(2*nn,bb)*gamma(2*nn)
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function get_sersicb, sersicn
common sersicb, nn
nn = sersicn
return, zbrent(0D,20D,func_name='sersicb_func',max_iterations=50)
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function hs_sersic_totflux, i_e, r_e, n_ser, q 

    compile_opt idl2
    on_error, 2
 
    if n_Params() LT 4 then begin 
        print, 'Syntax, flux = HS_SERSIC_TOTFLUX( I_e, R_e, n_ser, q )' 
        return, -1 
    endif 

    ; Based on equations from Ciotti & Bertin 1999
    flux_tot = 2.0D * !Pi * n_ser * GAMMA( 2.0 * n_ser ) * $
        EXP( get_sersicb( n_ser ) ) * q * ( r_e^2.0 ) * I_e / $
        ( ( get_sersicb( n_ser ) )^( 2.0 * n_ser ) )
    return, flux_tot 

end
