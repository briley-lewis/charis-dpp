;+
; NAME:
;   GTITRIM
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;
; PURPOSE:
;   Normalize a Good Time Interval (GTI) - no overlapping and adjoining
;
; CALLING SEQUENCE:
;   NEWGTI = GTITRIM(GTI, COUNT=, MAXGAP=, MINGTI=)
;
; DESCRIPTION:
;
;   A good time interval is by definition a list of intervals which
;   represent "good" or acceptable portions of the real number line.
;   In this library a GTI is a 2xNINTERVAL array where NINTERVAL is
;   the number of intervals.
;
;   The numbers in the array represent the start and stop times of
;   each interval.  Thus, the array [[0,10],[20,30],[40,50]] represent
;   intervals ranging from 0-10, 20-30 and 40-50.  Formally, this
;   example GTI represents times which would satisfy the following
;   expression, for each time T and interval number i:
;       T GE GTI(0,i) AND T LT GTI(1,i)
;   Note that the endpoint is closed on the left but open on the
;   right.
;
;   However, not every 2xNINTERVAL array is a valid or "normalized"
;   GTI as used by this library.  The array must satisfy several
;   conditions:
;      * time ordered (ascending)
;      * no overlapping intervals
;      * no adjoining intervals (intervals that start and stop at the
;        same point; e.g. the point 10 in this array [[0,10],[10,20]])
;
;   A user who desires to create his or her own GTI array can proceed
;   as follows.
;
;   First, the array is placed in time order.  This can be
;   accomplished simply using the built-in function SORT.  This
;   statement sorts the array by start times.
;
;       GTI = GTI(*, SORT(GTI(0,*)))
;
;   Second, the GTITRIM function is used to fully normalize the set of
;   intervals:
;
;       GTI = GTITRIM(GTI)
;
;   After these two procedures the GTI is considered valid and can be
;   passed to the other routines of the library.  Of course if the
;   user can guarantee the above requirements without using GTITRIM
;   then this is acceptable as well.
;
;   It should be noted that this function is not constrained to
;   operation only on time arrays.  It should work on any
;   one-dimensional quantity with intervals.
;
; INPUTS:
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
;   MAXGAP - Maximum allowable gap for merging existing good time
;            intervals.  Intervals with gaps smaller than MAXGAP will
;            be combined into a single interval.
;            Default: 0   (any gap keeps intervals separate)
;
;   MINGTI - Minimum size interval.  If any interval is smaller than
;            MINGTI then it is discarded.
;            Default: 0   (all intervals are preserved)
;
;   COUNT - upon return, the number of resulting intervals.  A value
;           of zero indicates no good time intervals.
;
;
; RETURNS:
;
;   A new GTI array containing the normalized intervals.  The array is
;   2xCOUNT where COUNT is the number of resulting intervals.
;   GTI(*,i) represents the start and stop times of interval number i.
;   The intervals are non-overlapping and time-ordered.
;
;   If COUNT is zero then the returned array is a scalar value of
;   zero, indicating no good intervals were found.
;
; SEE ALSO:
;
;   GTIMERGE
;
; MODIFICATION HISTORY:
;   Written, CM, 1997-2001
;   Documented, CM, Apr 2001
;   Corrected bug which bypassed MIN/MAXGTI, CM, 20 Jul 2003
;
;  $Id: gtitrim.pro,v 1.7 2006/10/22 09:50:09 craigm Exp $
;
;-
; Copyright (C) 1997-2001, 2003, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
; -

function gtitrim, gti, maxgap = maxgap, mingti = mingti, maxgti = maxgti, $
  count = count, query = query
  compile_opt idl2

  if keyword_set(query) then return, 1

  count = 0l
  ngti = n_elements(gti)
  if n_elements(maxgap) eq 0 then maxgap = 0d
  if n_elements(mingti) eq 0 then mingti = 0d

  ; ; Case of empty GTI
  if ngti lt 2 then return, 0l

  ; ; Case of single GTI
  if ngti eq 2 then begin
    if gti[1] - gti[0] lt mingti then return, 0l
    count = 1l
    newgti = reform(gti, 2, 1)
    goto, gti_checks
  endif

  newgti = gti
  ngti = ngti / 2
  ; ; Check for gaps, spaces between GTIs
  gaps = newgti[0, 1 : *] - newgti[1, *]
  wh = where(gaps gt maxgap, ct)
  if ct eq 0 then begin
    ; ; All of the gaps are too small
    newgti = reform([min(newgti), max(newgti)], 2, 1)
  endif else if ct gt 0 and ct lt ngti - 1 then begin
    ; ; Remake the GTI, removing the gaps
    vgti = reform(make_array(2, ct + 1, value = newgti[0]), 2, ct + 1, /overwrite)
    vgti[0, 0] = min(newgti)
    vgti[1, ct] = max(newgti)
    vgti[1, 0 : ct - 1] = newgti[1, wh]
    vgti[0, 1 : ct] = newgti[0, wh + 1]
    newgti = temporary(vgti)
  endif

  gti_checks:
  ; ; Remove too-small GTIs
  dur = newgti[1, *] - newgti[0, *]
  wh = where(dur ge mingti, count)
  if count gt 0 then newgti = newgti[*, wh] else newgti = 0l
  if count gt 0 then newgti = reform(newgti, 2, count, /overwrite)

  ; ;
  if count gt 0 and n_elements(maxgti) gt 0 then begin
    maxgti1 = maxgti[0]
    nper = (newgti[1, *] - newgti[0, *]) / maxgti1
    iper = ceil(nper)
    totgti = total(iper)
    if maxgti1 gt 0 and totgti gt count then begin
      vzero = newgti[0]
      vzero[0] = 0
      newgti2 = make_array(2, totgti, value = vzero)
      j = 0l
      for i = 0l, count - 1 do begin
        if iper[i] eq 1 then begin
          newgti2[*, j] = newgti[*, i]
        endif else begin
          for k = 0, iper[i] - 1 do $
            newgti2[*, j + k] = [newgti[0, i] + maxgti1 * k, $
              newgti[1, i] < (newgti[0, i] + maxgti1 * (k + 1))]
        endelse
        j = j + iper[i]
      endfor

      count = totgti
      newgti = reform(newgti2, 2, totgti)
    endif
  endif

  return, newgti
end