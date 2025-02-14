FUNCTION charis_pdi_mueller_rotation, theta
;+
; INPUTS:
;	theta = float, degrees, component angle for rotation matrix 
; OUTPUTS:
;	T = array, 4x4 floating point rotation matrix
;-
;--------------------------------------------------------------------
;
    trad = theta * (!DPI / 180.d) ; theta in radians
    T = make_array(4,4,/double)
    T[3,3] = 1
    T[0,0] = 1
    T[1,1] = cos(2.*trad)
    T[2,2] = cos(2.*trad)
    T[2,1] = sin(2.*trad)
    T[1,2] = -1.*T[2,1]
    return, T
END