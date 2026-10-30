function compute_average_angle, $
  angles, $
  dimension = dimension, $
  radians = radians$
  compile_opt idl2

  deg2rad = !dpi / 180.d
  rad2deg = 180.d / !dpi
  if ~keyword_set(radians) then angs = angles * deg2rad else angs = angles

  x = mean((cos(angs)), dimension = dimension)
  y = mean((sin(angs)), dimension = dimension)
  ang_avg = (atan(y, x) + (2d * !dpi)) mod (2d * !dpi)
  if ~keyword_set(radians) then ang_avg = ang_avg * rad2deg
  RETURN, ang_avg
end