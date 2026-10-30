function charis_pdi_mueller_matrix, $
  epsilon_bs, $
  epsilon_der, delta_h_der, delta_45_der, delta_r_der,$
  epsilon_hwp, delta_hwp, $
  epsilon_tel, delta_tel, $
  theta_der, theta_hwp, $
  altitude, parang
  compile_opt idl2
  ;+
  ; INPUTS:
  ;
  ; ...
  ; OUTPUTS:
  ;	M = array, 4x4 floating point Mueller matrix
  ;-
  ; --------------------------------------------------------------------
  Mbs_left = charis_pdi_mueller_bs(epsilon_bs, -1d)
  Mbs_right = charis_pdi_mueller_bs(epsilon_bs, 1d)

  Mder = charis_pdi_mueller_elliptical_retarder(epsilon_der, delta_h_der, delta_45_der, delta_r_der)
  Mhwp = charis_pdi_mueller_component(epsilon_hwp, delta_hwp)
  Mtel = charis_pdi_mueller_component(epsilon_tel, delta_tel)

  Tmder = charis_pdi_mueller_rotation(-theta_der)
  Tder = charis_pdi_mueller_rotation(theta_der)
  Tmhwp = charis_pdi_mueller_rotation(-theta_hwp)
  Thwp = charis_pdi_mueller_rotation(theta_hwp)
  Talt = charis_pdi_mueller_rotation(-altitude)
  Tpa = charis_pdi_mueller_rotation(parang)

  rMder = Tmder ## (Mder ## Tder)
  rMhwp = Tmhwp ## (Mhwp ## Thwp)
  rMtel = Talt ## (Mtel ## Tpa)
  M0 = rMder ## (rMhwp ## rMtel) ; note we don't need the diag matrix in McIntosh+2026 because of the already applied parity flip
  M_right = Mbs_right ## M0
  M_left = Mbs_left ## M0

  mueller_matrix = M_right - M_left
  return, mueller_matrix
end