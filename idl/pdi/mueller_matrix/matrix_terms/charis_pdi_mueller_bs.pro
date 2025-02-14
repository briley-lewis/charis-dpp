FUNCTION charis_pdi_mueller_bs, epsilon_bs, half
;+
; INPUTS:
;	epsilon_bs = float, diattenuation of the Wollastron prism
;   half = int, +1 for right half or -1 for left half
; OUTPUTS:
;	M = array, 4x4 floating point matrix for Wollaston prism
;-
;--------------------------------------------------------------------
;
    M = make_array(4,4,/double)
    M[0,0] = 0.5
    M[1,1] = M[0,0]
    M[0,1] = 0.5*half*epsilon_bs
    M[1,0] = M[0,1]
    M[2,2] = 0.5*sqrt(1.-epsilon_bs^2)
    M[3,3] = M[2,2]
    return, M
END