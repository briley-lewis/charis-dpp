pro charis_set_epoch, pfname = pfname, mjd = mjd, header = header, date = date, quiet = quiet
  compile_opt idl2
  ; Set the session default epoch used by charis_get_constant when a caller
  ; does not pass mjd=.  [BUGS.md #322]
  ;
  ; Exactly one of
  ;   pfname=  dataset .info file; uses its OBSDATE ('YYYYMMDD')
  ;   mjd=     a Modified Julian Date
  ;   header=  a FITS header carrying MJD (or MJD-STR)
  ;   date=    'YYYYMMDD'
  ; The top-level reduction routines call charis_set_epoch, pfname=pfname so
  ; that every constant (pixel scale, north offset, satellite-spot contrast,
  ; M3 coefficients, parity, ...) matches the night being reduced rather than
  ; the latest epoch.  Call it yourself with mjd= to override in a session.
  common charis_epoch, charis_epoch_mjd

  mjd_new = -1d
  if n_elements(mjd) gt 0 then begin
    mjd_new = double(mjd[0])
  endif else if n_elements(header) gt 0 then begin
    mjd_new = double(sxpar(header, 'MJD', count = c))
    if c eq 0 then mjd_new = double(sxpar(header, 'MJD-STR', count = c))
    if c eq 0 then message, 'charis_set_epoch: header has no MJD or MJD-STR'
  endif else begin
    if n_elements(date) eq 0 then begin
      if n_elements(pfname) eq 0 then message, 'charis_set_epoch: give pfname=, mjd=, header= or date='
      param, 'obsdate', date, /get, pfname = pfname
    endif
    d = strtrim(string(date), 2)
    if strlen(d) ne 8 then message, 'charis_set_epoch: OBSDATE must be YYYYMMDD, got ' + d
    yyyy = fix(strmid(d, 0, 4)) & mm = fix(strmid(d, 4, 2)) & dd = fix(strmid(d, 6, 2))
    mjd_new = julday(mm, dd, yyyy, 0, 0, 0) - 2400000.5d
  endelse

  if mjd_new lt 50000d or mjd_new gt 70000d then message, 'charis_set_epoch: implausible MJD ' + string(mjd_new)
  charis_epoch_mjd = mjd_new
  if ~keyword_set(quiet) then $
    print, 'charis_set_epoch: MJD ', strtrim(string(mjd_new, format = '(f10.2)'), 2), $
    ' -> epoch "', charis_get_constant(name = 'epoch', mjd = mjd_new), '"'
end
