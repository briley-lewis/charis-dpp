FUNCTION s_reflection_coefficient, cos_i, cos_t, n1, n2
;+
; INPUTS:
;	cos_i = float, the cosine of the angle of incidence
;   cos_t = float, the cosine of the angle of transmission
;   n1 = float, refractive index of material 1
;   n2 = float, refractive index of material 2
; ...
; OUTPUTS:
;	r_s = float, the Fresnel reflection coefficient for s polarization
;-
;--------------------------------------------------------------------
;
r_s = (n1*cos_i - n2*cos_t) / (n1*cos_i + n2*cos_t)
return, r_s
end