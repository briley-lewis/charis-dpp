;+
; NAME:
;   GTI2MASK
;
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;
; PURPOSE:
;   Convert Good Time Interval (GTI) to evenly sampled mask array
;
; CALLING SEQUENCE:
;   MASK = GTI2MASK(GTI, [TIME, TLIMITS=, TIMEDEL=, NTBINS=,
;                   GOOD=, BAD=, /FILL, /INVERT])
;
; DESCRIPTION:
;
;   The function GTI2MASK converts an existing valid Good Time
;   Interval (GTI) array to a mask array.  By definition a GTI
;   indicates an array of intervals which are not on an evenly sampled
;   array.  This routine either accepts a time grid from the user, or
;   the grid must be described by the TLIMITS and TIMEDEL keywords.
;
;   The output mask array describes whether each grid point lies
;   within a good interval or not.  There is full control over the
;   values of the good and bad values.
;
;   This routine is the inverse of MASK2GTI.
;
;   It should be noted that this function is not constrained to
;   operation only on time arrays.  It should work on any
;   one-dimensional quantity with intervals.
;
; INPUTS:
;
;   GTI - a 2xNINTERVAL array where NINTERVAL is the number of
;         intervals.  GTI(*,i) represent the start and stop times of
;         interval number i.  The intervals must be non-overlapping
;         and time-ordered (use GTITRIM to achieve this).
;
;         A scalar value of zero indicates that the GTI is empty, ie,
;         there are no good intervals.
;
;   TIME - optional time array that specifies the time grid for the
;          mask array.  If TIME is not specified then the user must
;          give the TLIMITS and TIMEDEL keywords to fully describe the
;          grid spacing.  The TIME array is overwritten if the FILL
;          keyword is used.
;
; KEYWORDS:
;
;   TLIMITS - a 2-element array giving the start and stop limits over
;             which the mask array is to be generated.  The TLIMITS
;             and TIMEDEL keywords are required if the TIME parameter
;             is not given.
;
;   TIMEDEL - a scalar specifying the interval between grid points.
;             The TLIMITS and TIMEDEL keywords are required if the
;             TIME parameter is not given.
;
;   NTBINS - upon return, this keyword contains the number of time
;            samples created.
;
;   GOOD - the value of "good" in the output mask array.
;          Default: 1b
;
;   BAD - the value of "bad" in the output mask array.
;         Default: 0b
;
;   INVERT - if set, the array GTI is treated as a "bad" time
;            interval, ie, the GOOD and BAD values are swapped.
;
;   FILL - if set, the array TIME is filled with values determined
;          from the TLIMITS and TIMEDEL keyword.
;
; RETURNS:
;
;   A mask array, either sampled at the points specified by TIME, or
;   by the grid specified by TLIMITS and TIMEDEL.  The "good" value
;   indicates that the point lies within the good interval, while a
;   "bad" value indicates the point was outside.
;
; SEE ALSO:
;
;   MASK2GTI, GTITRIM, GTIMERGE
;
; MODIFICATION HISTORY:
;   Written, CM, 1997-2001
;   Documented, CM, Apr 2001
;   Add internal OVERLAP and MINFRACEXP keywords, CM, 03 Feb 2007
;   Refine and simplify the OVERLAP processing, CM, 14 Feb 2007
;   Handle case of /OVERLAP when there is no intersection, CM, 22 Aug 2007
;   Use VALUE_LOCATE for performance (in non-OVERLAP case), CM, 04 May 2008
;
;  $Id: gti2mask.pro,v 1.7 2008/05/04 21:35:27 craigm Exp $
;
;-
; Copyright (C) 1997-2001, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
; -

function gti2mask, gti, time, fill = fill, query = query, $
  tlimits = tlimits, good = good, bad = bad, invert = invert, $
  timepixr = timepixr, timedel = timedel, ntbins = nbins, $
  overlap = overlap, minfracexp = minfracexp, exposure = expo
  compile_opt idl2

  if keyword_set(query) then return, 1
  if keyword_set(overlap) and n_elements(timedel) eq 0 then begin
    message, 'ERROR: must specify OVERLAP and TIMEDEL together'
    return, 0
  endif

  if n_elements(good) eq 0 then good = 1b
  if n_elements(bad) eq 0 then bad = 0b
  if n_elements(invert) eq 0 then invert = 0
  if n_elements(tlimits) lt 2 and n_elements(gti) lt 2 then begin
    message, 'ERROR: could not determine TLIMITS'
    return, 0l
  endif
  if n_elements(minfracexp) eq 0 then minfracexp = 0d
  if n_elements(tlimits) lt 2 then tlimits = [min(gti), max(gti)]
  if n_elements(timepixr) eq 0 then timepixr = 0.d
  if n_elements(timedel) eq 0 then timedel = 1.d
  nbins = long((tlimits[1] - tlimits[0]) / timedel)
  if n_elements(time) eq 0 or keyword_set(fill) then $
    time = dindgen(nbins) * timedel + tlimits[0]
  nbins = n_elements(time)

  if keyword_set(invert) then rep = good else rep = bad
  if keyword_set(invert) then new = bad else new = good
  mask = make_array(nbins, value = rep)
  if n_elements(gti) lt 2 then begin
    if keyword_set(overlap) then expo = dblarr(nbins)
    return, mask
  endif

  ngti = n_elements(gti) / 2
  nt = n_elements(time)

  if keyword_set(overlap) then begin
    expo = dblarr(nbins)
    i0 = value_locate([time, tlimits[1]], gti[0, *])
    i1 = value_locate([time, tlimits[1]], gti[1, *])
    for i = 0l, ngti - 1 do begin
      tmin = gti[0, i]
      tmax = gti[1, i]
      i0i = i0[i]
      i1i = i1[i]
      if i0i ge nt or i1i le -1 then goto, next_gti_1
      if i0i eq -1 then begin
        tmin = tlimits[0]
        i0i = 0
      endif
      if i1i eq nt then begin
        tmax = tlimits[1]
        i1i = nt - 1
      endif
      expo[i0i] = expo[i0i : i1i] + timedel
      expo[i0i] = expo[i0i] + (time[i0i] - tmin)
      expo[i1i] = expo[i1i] - (time[i1i] + timedel - tmax)
      next_gti_1:
    endfor
  endif else begin
    ; ; Arghhh!  One-element arrays confuse VALUE_LOCATE
    if n_elements(gti) eq 2 then begin
      ii = (time ge gti[0]) - 1
    endif else begin
      ii = value_locate(gti[0, *], time)
    endelse
    wh = where(ii ge 0 and time lt gti[1, ii], ct)
    if ct gt 0 then mask[wh] = new
  endelse

  ; ; Finally determine where overlap is larger than required
  if keyword_set(overlap) then begin
    mask1 = expo gt timedel * minfracexp
    wh1 = where(mask1, ct1)
    if ct1 gt 0 then mask[wh1] = good
    wh1 = where(mask1 eq 0, ct1)
    if ct1 gt 0 then mask[wh1] = bad
  endif

  return, mask
end