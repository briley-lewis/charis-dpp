FUNCTION complex_refractive_index, eps_r1, eps_r2
;+
; INPUTS:
;	eps_r1 = float, dielectric constant 1 for the material
;   eps_r2 = float, dielectric constant 2 for the material
; ...
; OUTPUTS:
;	n_complex = float, the complex refractive index of the material
;-
;--------------------------------------------------------------------
;
eps_rmag = sqrt(eps_r1^2d + eps_r2^2d)
scalar = (1d/sqrt(2d))
n = scalar*sqrt(eps_rmag + eps_r1)
k =  scalar*sqrt(eps_rmag - eps_r1)
n_complex = complex(n, k)
return, n_complex
end