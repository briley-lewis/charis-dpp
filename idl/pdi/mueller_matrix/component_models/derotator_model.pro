FUNCTION derotator_model, lambda, epsilon_der=epsilon_der

d = 262.56d ;width of Quartz layer in nm, from Hart+2021

n0 = 1.0d

; Quartz index of refraction. Hard coding these for now until I can access the functional forms...
n1_real = double([1.466673, 1.4665226, 1.4663808, 1.466249, 1.4661282, 1.4660134, 1.4659056000000001, 1.465806, 1.4657124, 1.465625, 1.4655444, 1.4654690000000001, 1.4653972, 1.465331, 1.465269, 1.4652112000000002, 1.4651576, 1.4651068, 1.4650596, 1.4650158, 1.464974, 1.4649364])
k1 = 0.d
n1 = complex(n1_real, k1)

; Silver
n2_real = double([0.27433, 0.28631, 0.29935, 0.31319, 0.32788, 0.34373, 0.3609, 0.37905, 0.3986, 0.41961, 0.44179, 0.46559, 0.49152, 0.5188, 0.5485, 0.57968, 0.61342, 0.64938, 0.68823, 0.72949, 0.77389, 0.82102])
k2 = double([7.4359, 7.7047, 7.9867, 8.2744, 8.569, 8.8761, 9.1961, 9.5224, 9.8619, 10.214, 10.572, 10.943, 11.333, 11.73, 12.146, 12.569, 13.01, 13.464, 13.938, 14.424, 14.93, 15.448])
n2 = complex(n2_real, k2)

r1 = single_mirror_rtot(lambda, d, n0, n1, n2, !DPI / 3.d)
rs1 = r1[*,0]
rp1 = r1[*,1]

r2 = single_mirror_rtot(lambda, d, n0, n1, n2, !DPI / 6.d)
rs2 = r2[*,0]
rp2 = r2[*,1]

rs3 = rs1
rp3 = rp1

rs_tot = rs1*rs2*rs3
rp_tot = rp1*rp2*rp3

arg_rs = atan(rs_tot, /phase)
arg_rp = atan(rp_tot, /phase)
delta_der = ((arg_rs - arg_rp) + 2.d * !DPI) mod (2.d * !DPI)
delta_der = delta_der*!RADEG

epsilon_der = (abs(rs_tot)^2 - abs(rp_tot)^2) / (abs(rs_tot)^2 + abs(rp_tot)^2)

return, delta_der
end