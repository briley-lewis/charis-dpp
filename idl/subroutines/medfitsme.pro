function medfitsme, flist, mean = mean, dir = directory, cube = cube, nwght = nwght, varcube = varcube
  compile_opt idl2
  ; [BUGS.md #154] alias of charis_combinefits (the routine called by
  ; charis_ssdi.pro was never shipped).
  return, charis_combinefits(flist, mean = mean, dir = directory, cube = cube, nwght = nwght, varcube = varcube)
end
