FUNCTION charis_pdi_mueller_component, epsilon, delta
;+
; INPUTS:
;	epsilon = float, linear diattenuation (ideally 0.0); [-1.0, 1.0]
;   delta = float, linear retardance in degrees (ideally 180.0); [0.0, 360.0]
; OUTPUTS:
;	M = array, 4x4 floating point component Mueller matrix (Mtel, Mhwp, Mder)
;-
;--------------------------------------------------------------------
;
    delta_rad = delta * (!DPI / 180.d)
    M = make_array(4,4,/double)
    M[0,0] = 1
    M[1,1] = 1
    M[0,1] = epsilon
    M[1,0] = epsilon
    M[2,2] = sqrt(1-epsilon^2) * cos(delta_rad)
    M[3,3] = M[2,2]
    M[3,2] = sqrt(1-epsilon^2) * sin(delta_rad)
    M[2,3] = -1*M[3,2]
    return, M
END