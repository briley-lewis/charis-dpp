function charis_pdi_mueller_elliptical_retarder, epsilon_der, delta_h, delta_45, delta_r
  compile_opt idl2
  ;+
  ; Updated derotator model from McIntosh+2026
  ; INPUTS:
  ;   epsilon_der = float, diattenuation of the derotator
  ;   delta_h  = float, horizontal linear retardance in radians
  ;   delta_45 = float, 45-degree linear retardance in radians
  ;   delta_r  = float, right-handed circular retardance in radians
  ; OUTPUTS:
  ;   mm = array, 4x4 floating point component Mueller matrix
  ;-

  dh_rad = delta_h
  dr_rad = delta_r
  d45_rad = delta_45

  phi_sq = dh_rad^2 + dr_rad^2 + d45_rad^2
  phi = sqrt(phi_sq)

  mm = dblarr(4, 4)
  mm[0, 0] = 1.d

  if (phi eq 0.d) then begin
    mm[1, 1] = 1.d
    mm[2, 2] = 1.d
    mm[3, 3] = 1.d
    return, mm
  endif

  c = cos(phi)
  s = sin(phi)
  t = 1.d - c

  two_two   = (dh_rad^2 + (d45_rad^2 + dr_rad^2) * c) / phi_sq
  two_three = (d45_rad * dh_rad * t) / phi_sq + (dr_rad * s) / phi
  two_four  = (dr_rad * dh_rad * t) / phi_sq - (d45_rad * s) / phi

  three_two   = (dh_rad * d45_rad * t) / phi_sq - (dr_rad * s) / phi
  three_three = (d45_rad^2 + (dh_rad^2 + dr_rad^2) * c) / phi_sq
  three_four  = (dr_rad * d45_rad * t) / phi_sq + (dh_rad * s) / phi

  four_two   = (dh_rad * dr_rad * t) / phi_sq + (d45_rad * s) / phi
  four_three = (d45_rad * dr_rad * t) / phi_sq - (dh_rad * s) / phi
  four_four  = (dr_rad^2 + (dh_rad^2 + d45_rad^2) * c) / phi_sq

  
  ; Row 2 (Index 1)
  mm[1, 1] = two_two
  mm[2, 1] = two_three   ; M23 -> Column index 2, Row index 1
  mm[3, 1] = two_four    ; M24 -> Column index 3, Row index 1

  ; Row 3 (Index 2)
  mm[1, 2] = three_two   ; M32 -> Column index 1, Row index 2
  mm[2, 2] = three_three
  mm[3, 2] = three_four  ; M34 -> Column index 3, Row index 2

  ; Row 4 (Index 3)
  mm[1, 3] = four_two    ; M42 -> Column index 1, Row index 3
  mm[2, 3] = four_three  ; M43 -> Column index 2, Row index 3
  mm[3, 3] = four_four

  return, mm
end
