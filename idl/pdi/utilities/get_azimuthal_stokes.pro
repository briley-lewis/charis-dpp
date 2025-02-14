FUNCTION get_azimuthal_stokes, qarr, uarr, phi, uphi=uphi
sz = size(qarr, /dim)
if n_elements(sz) eq 3 then begin
    phiarr = fltarr(sz[0],sz[1],sz[2])
    for Li=0,sz[2]-1 do phiarr[*,*,Li] = phi
endif else begin
    phiarr = phi
endelse
qphi = (-qarr * cos(2.*phiarr)) - (uarr * sin(2.*phiarr))
uphi = (qarr * sin(2.*phiarr)) - (uarr * cos(2.*phiarr))
return, qphi
end