;+
; NAME:
;   GTIWHERE
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;
; PURPOSE:
;   Locate a list of times within a Good Time Interval (GTI) array
;
; CALLING SEQUENCE:
;   WH = GTIWHERE(TIME, GTI, COUNT=, INTERVALS=, /INVERT, /INCLUDE)
;
; DESCRIPTION:
;
;   The function GTIWHERE is an efficient method to determine which
;   good time interval (GTI) a particular time sample falls into, if
;   any.  The user passes an array of one or more times in the TIME
;   parameter, and GTIWHERE determines which of these times are in a
;   good interval.  The INTERVALS keyword returns an identification of
;   which interval the time fell into.
;
;   It should be noted that this function is not constrained to
;   operation only on time arrays.  It should work on any
;   one-dimensional quantity with intervals.
;
;   The definition of "inside" a good-time interval is
;      TIME GE GTI(0,i) AND TIME LT GTI(1,i)     (INCLUDE=0)
;      TIME GE GTI(0,i) AND TIME LE GTI(1,i)     (INCLUDE=1)
;   where GTI(*,i) is the i'th interval.
;
;   When using /INVERT, the definition of "outside" a bad-time
;   interval is
;      TIME LT GTI(0,i) OR TIME GE GTI(1,i)    (INCLUDE=0)
;      TIME LT GTI(0,i) OR TIME GT GTI(1,i)    (INCLUDE=0)
;
; INPUTS:
;
;   TIME - an array of times, in no particular order.
;
;   GTI - a 2xNINTERVAL array where NINTERVAL is the number of
;         intervals.  GTI(*,i) represents the start and stop times of
;         interval number i.  The intervals must be non-overlapping
;         and time-ordered (use GTITRIM to achieve this).
;
;         A scalar value of zero indicates that the GTI is empty, ie,
;         there are no good intervals.
;
; KEYWORDS:
;
;   COUNT - upon return, the number of resulting intervals.  A value
;           of zero indicates no good time intervals.
;
;   INTERVALS - upon return, an array of integers specifying the
;               interval number each time falls into.  The number of
;               elements of INTERVALS is COUNT (ie, there is one
;               INTERVAL for each returned index).
;
;   INCLUDE - if set then the endpoints are considered closed at both
;             ends; by convention the endpoints are normally treated
;             as closed at the left and open at the right.
;
;   INVERT - if set, the array GTI is treated as a "bad" time
;            interval.  Only times *outside* of the GTI are selected.
;            In this case there are NGTI+1 possible intervals, where
;            NGTI is the number of time intervals passed in GTI.  The
;            indices returned in the INTERVALS keyword start from the
;            "left."
;
; RETURNS:
;
;   An array of integer indices indicating which elements of TIME are
;   within the "good" intervals.  The number of selected elements is
;   COUNT.  If COUNT is zero, then a scalar value of -1L is returned,
;   indicating no selected values.
;
;
; SEE ALSO:
;
;   GTISEG, GTIMERGE, GTIENLARGE
;
; MODIFICATION HISTORY:
;   Written, CM, 1997-2001
;   Documented, CM, Apr 2001
;   Added usage message, CM, 2006 Aug 18
;   Handle the case if /INVERT and /INCLUDE, which changes the
;     boundary conditions slightly, CM, 2007 Dec 15
;   Handle case of empty input GTI and /INVERT; also handle
;     INVERT in case where MIN(TIME) touches MIN(GTI), CM, 2008 Jul 08
;
;  $Id: gtiwhere.pro,v 1.6 2008/07/08 20:24:18 craigm Exp $
;
;-
; Copyright (C) 1997-2001, 2006, 2007, 2008 Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
; -
function gtiwhere, time, gti, count = count, intervals = vv, $
  invert = invert, include = include, query = query
  compile_opt idl2

  if keyword_set(query) then return, 1

  count = 0l
  vv = -1l

  if n_params() eq 0 then begin
    message, 'USAGE: GTIWHERE(TIME, GTI, COUNT=, INTERVALS=, /INVERT, /INCLUDE)', /info
    goto, finish
  endif

  if n_elements(gti) lt 2 then begin
    if keyword_set(invert) then begin
      ; ; All times are accepted when /INVERT is set and GTI is null
      count = n_elements(time)
      vv = lonarr(count)
      wh = lindgen(count)
    endif
    goto, finish
  endif

  ngti = n_elements(gti) / 2

  if keyword_set(invert) then begin
    vgti = [min([gti[*], time - 1]), gti[*], max([gti[*], time + 1])]
    ngti = ngti + 1
  endif else begin
    vgti = gti
  endelse
  vgti = reform(vgti, 2, ngti, /overwrite)

  if ngti lt 5 then begin
    ; ; Optimization for a few GTIs, which is faster by using WHERE
    ; ; rather than VALUE_LOCATE
    mint = min(time, max = maxt)
    vv = lonarr(n_elements(time))

    for i = 0l, ngti - 1 do begin
      if mint gt vgti[1, i] then goto, next_interval
      if maxt lt vgti[0, i] then goto, next_interval
      if keyword_set(include) and keyword_set(invert) then begin
        wh1 = where(time gt vgti[0, i] and time lt vgti[1, i], ct1)
      endif else if keyword_set(include) then begin
        wh1 = where(time ge vgti[0, i] and time le vgti[1, i], ct1)
      endif else begin
        wh1 = where(time ge vgti[0, i] and time lt vgti[1, i], ct1)
      endelse

      if ct1 gt 0 then begin
        if count eq 0 then wh = wh1 else wh = [wh, wh1]
        count = count + ct1
        vv[wh1] = i
      endif
      next_interval:
    endfor
  endif else begin
    ; ; Standard uses VALUE_LOCATE which is faster for all but the
    ; ; cases NGTI between 1 and 5.
    vv = value_locate(vgti[0, *], time)
    if keyword_set(include) and keyword_set(invert) then begin
      wh = where((vv gt 0) and (time lt vgti[1, vv > 0]), count)
    endif else if keyword_set(include) then begin
      wh = where((vv ge 0) and (time le vgti[1, vv > 0]), count)
    endif else begin
      wh = where((vv ge 0) and (time lt vgti[1, vv > 0]), count)
    endelse
  endelse

  finish:
  if count le 0 then return, 0l
  vv = vv[wh]
  return, wh
end