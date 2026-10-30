function get_xycind, dimx, dimy, xc, yc, rad
  compile_opt idl2
  ; [BUGS.md #79] 1-D indices of the pixels within `rad` of (xc, yc) in a
  ; dimx x dimy image (a GPI-pipeline routine that was never shipped with the
  ; CHARIS DPP).  Returns -1 if none.
  x = findgen(dimx) # replicate(1., dimy)
  y = replicate(1., dimx) # findgen(dimy)
  r = sqrt((x - xc) ^ 2 + (y - yc) ^ 2)
  return, where(r le rad)
end
