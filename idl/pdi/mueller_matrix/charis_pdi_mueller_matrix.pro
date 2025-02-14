FUNCTION charis_pdi_mueller_matrix,$
         epsilon_bs,$
         epsilon_der, delta_der,$
         epsilon_hwp, delta_hwp,$
         epsilon_tel, delta_tel,$
         theta_der, theta_hwp,$
         altitude, parang
;+
; INPUTS:
;	
; ...
; OUTPUTS:
;	M = array, 4x4 floating point Mueller matrix
;-
;--------------------------------------------------------------------
    Mbs_left = charis_pdi_mueller_bs(epsilon_bs, -1d)
    Mbs_right = charis_pdi_mueller_bs(epsilon_bs, 1d)   

    Mder = charis_pdi_mueller_component(epsilon_der, delta_der)
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
    M0 = rMder ## (rMhwp ## rMtel)
    M_right = Mbs_right ## M0
    M_left = Mbs_left ## M0

    mueller_matrix = M_right - M_left
    return, mueller_matrix
END