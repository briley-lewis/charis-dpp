FUNCTION adjust_pol_axis,qcube,ucube,da,phi,ucube_out=ucube_out
    sz = size(qcube, /dim)
    qphicube = get_azimuthal_stokes(qcube, ucube, phi, uphi=uphicube)
    adj_phi = (da*!DTOR)+phi
    phicube = fltarr(sz[0],sz[1],sz[2])
    for Li=0,sz[2]-1 do phicube[*,*,Li] = adj_phi
    qcube_out = uphicube*sin(2*phicube) - qphicube*cos(2*phicube)
    ucube_out = -1*uphicube*cos(2*phicube) - qphicube*sin(2*phicube)
    return, qcube_out
END
         