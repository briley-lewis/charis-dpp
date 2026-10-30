function calc_centroid_spots4, image, xguess, yguess, fwhm, boxsize = boxsize
  compile_opt idl2
  ; [BUGS.md #51] 2-D Gaussian centroid of a satellite spot (called by
  ; charis_register_cube with /fitgauss; never shipped).  Fits an elliptical
  ; Gaussian with mpfit2dpeak inside a box of `boxsize` (default 4*fwhm)
  ; pixels around the guess and returns [x, y]; [-1, -1] on failure.
  if ~keyword_set(boxsize) then boxsize = long(4 * fwhm) > 7
  half = boxsize / 2
  sz = size(image, /dim)
  x0 = (long(round(xguess)) - half) > 0
  y0 = (long(round(yguess)) - half) > 0
  x1 = (long(round(xguess)) + half) < (sz[0] - 1)
  y1 = (long(round(yguess)) + half) < (sz[1] - 1)
  sub = image[x0 : x1, y0 : y1]
  bad = where(~finite(sub), nbad)
  if nbad gt 0 then sub[bad] = median(sub, /even)
  catch, err
  if err ne 0 then begin
    catch, /cancel
    return, [-1., -1.]
  endif
  yfit = mpfit2dpeak(sub, a, /gaussian, /tilt, /quiet)
  catch, /cancel
  if n_elements(a) lt 6 then return, [-1., -1.]
  xc = a[4] + x0
  yc = a[5] + y0
  if abs(xc - xguess) gt boxsize or abs(yc - yguess) gt boxsize then return, [-1., -1.]
  return, [xc, yc]
end
