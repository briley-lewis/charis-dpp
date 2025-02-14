FUNCTION single_mirror_rtot, lambda, d, n0, n1, n2, th_i01

cos_i01 = cos(th_i01)
cos_t01 = cosine_transmitted(th_i01, n0, n1)

beta = 2.*(2.*!PI/lambda)*d*n1*cos_t01
e_ibeta = exp(complex(0,1)*beta) ; This is e^(i*beta)

rs01 = s_reflection_coefficient(cos_i01, cos_t01, n0, n1)
ts01 = s_transmission_coefficient(cos_i01, cos_t01, n0, n1)
rp01 = p_reflection_coefficient(cos_i01, cos_t01, n0, n1)
tp01 = p_transmission_coefficient(cos_i01, cos_t01, n0, n1)

cos_i10 = cos_t01
th_i10 = acos(real_part(cos_i10))
cos_t10 = cosine_transmitted(th_i10, n1, n0)

ts10 = s_transmission_coefficient(cos_i10, cos_t10, n1, n0)
rs10 = s_reflection_coefficient(cos_i10, cos_t10, n1, n0)
tp10 = p_transmission_coefficient(cos_i10, cos_t10, n1, n0)
rp10 = p_reflection_coefficient(cos_i10, cos_t10, n1, n0)

cos_i12 = cos_t01
th_i12 = acos(real_part(cos_i12))
cos_t12 = cosine_transmitted(th_i12, n1, n2)

rs12 = s_reflection_coefficient(cos_i12, cos_t12, n1, n2)
rp12 = p_reflection_coefficient(cos_i12, cos_t12, n1, n2)

rs_tot = rs01+(ts01*ts10*rs12*e_ibeta)/(1.-rs10*rs12*e_ibeta)
rp_tot = rp01+(tp01*tp10*rp12*e_ibeta)/(1.-rp10*rp12*e_ibeta)
return, [[rs_tot], [rp_tot]]
end