FUNCTION get_charis_mueller_matrix,$
         h1, hwpang, epsilon_bs,$
         epsilon_der, delta_der,$
         epsilon_hwp, delta_hwp,$
         epsilon_tel, delta_tel,$
         dtheta_der, dtheta_hwp,$
         no_adisync = no_adisync,$
         theta_hwp_offset=theta_hwp_offset

    nL = sxpar(h1, 'naxis3')
    h3 = headfits(('./data/raw/CRSA'+(sxpar(h1, 'origname'))+'_cube.fits'), ext=3)

    altitude = sxpar(h1,'altitude')
    parang = sxpar(h1,'parang')
    theta_der = sxpar(h3, 'D_IMRANG')

    if keyword_set(no_adisync) then begin
        theta_hwp = hwpang
    endif else begin 
        deg2rad = !DPI / 180.d
        rad2deg = 180.d / !DPI
        
        a_tel = 19.823806d * deg2rad ; Subaru's latitude in radians
        azimuth = sxpar(h1,'azimuth')

        az_rad = (azimuth - 180.d) * deg2rad
        alt_rad = altitude * deg2rad

        arg1 = sin(az_rad)
        arg2 = sin(alt_rad)*cos(az_rad) + cos(alt_rad)*sin(a_tel)/cos(a_tel)
        fr = atan(arg1, arg2) * rad2deg
        theta_hwp_offset  = 0.5d*fr + altitude
        theta_hwp  = theta_hwp_offset + hwpang ; in degrees
    endelse

    theta_der_adj = theta_der + dtheta_der
    theta_hwp_adj = theta_hwp + dtheta_hwp

    M = make_array(nL, 4, /double)
    for Li=0,nL-1 do begin
        Mi = charis_pdi_mueller_matrix(epsilon_bs, epsilon_der[Li], delta_der[Li],$
                                        epsilon_hwp, delta_hwp[Li],$
                                        epsilon_tel[Li], delta_tel[Li],$
                                        theta_der_adj, theta_hwp_adj,$
                                        altitude, parang)
        M[Li, *] = Mi[*,0]
    endfor
    return, M
END