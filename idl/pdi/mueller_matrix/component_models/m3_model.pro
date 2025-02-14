FUNCTION m3_model, lambda, epsilon_m3 = epsilon_m3
;+
; INPUTS:
;	lambda = float, observing wavelength(s) in nm
; ...
; OUTPUTS:
;	delta_m3 = float, the approximate retardance of Subaru's M3 mirror
;                based on findings in Hart+2021
;	epsilon_m3 = float, the diattenuation of Subaru's M3 mirror
;                based on findings in Hart+2021
;-
;--------------------------------------------------------------------
;

lambda_meters = lambda * 1e-09

m1 = 2.104d
b1 = 14.20d
m2 = 2.100d
b2 = 13.20d

eps_r1 = -(10.d ^b1)*(lambda_meters^m1)
eps_r2 = (10.d ^b2)*(lambda_meters^m2)

n1 = 1.0d ; index of refraction for air assumed by Hart+2021
n2 = complex_refractive_index(eps_r1, eps_r2)

th_i = !DPI / 4.d

cos_i = cos(th_i)
cos_t = cosine_transmitted(th_i, n1, n2)

r_s = s_reflection_coefficient(cos_i, cos_t, n1, n2)
r_p = p_reflection_coefficient(cos_i, cos_t, n1, n2)

epsilon_m3 = (abs(r_s)^2 - abs(r_p)^2) / (abs(r_s)^2 + abs(r_p)^2)

arg_rs = atan(r_s, /phase)
arg_rp = atan(r_p, /phase)
delta_m3 = ((arg_rs - arg_rp)+2*!DPI) mod (2*!DPI)
delta_m3 = delta_m3*!RADEG
return, delta_m3
END