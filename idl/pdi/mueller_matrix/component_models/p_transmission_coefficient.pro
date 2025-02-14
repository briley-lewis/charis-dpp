FUNCTION p_transmission_coefficient, cos_i, cos_t, n1, n2
;+
; INPUTS:
;	cos_i = float, the cosine of the angle of incidence
;   cos_t = float, the cosine of the angle of transmission
;   n1 = float, refractive index of material 1
;   n2 = float, refractive index of material 2
; ...
; OUTPUTS:
;	t_p = float, the Fresnel transmission coefficient for p polarization
;-
;--------------------------------------------------------------------
;
t_p = 2.*n1*cos_i / (n2*cos_i + n1*cos_t)
return, t_p
end