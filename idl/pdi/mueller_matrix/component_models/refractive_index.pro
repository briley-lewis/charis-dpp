FUNCTION refractive_index, lambda, A, B, C, D, F
;+
; INPUTS:
;	lambda = float, wavelength ***IN MICRONS***
;   A,B,C,D,F = coefficients, e.g. from Ghosh 1999 for SiO2
; ...
; OUTPUTS:
;	n = float, the refractive index corresponding to the inputs
;-
;--------------------------------------------------------------------
;
L2 = lambda^2d
n_squared = A + ((B*L2)/(L2-C)) + ((D*L2)/(L2-F))
return, sqrt(n_squared)
END