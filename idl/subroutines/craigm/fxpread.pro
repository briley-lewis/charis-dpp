;+
; NAME:
;   FXPREAD
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE:
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Internal routine to read from a Unix pipe.
;
; DESCRIPTION:
;
;   See the following procedures for appropriate documentation.
;
;     FXGOPEN  - open resource
;     FXGCLOSE - close resource
;     FXGREAD  - read from resource
;     FXGWRITE - write to resource
;     FXGSEEK  - seek on resource (i.e., perform POINT_LUN)
;
;     FXGFILTERED - determine if resource is a normal file.
;
; PARAMETERS
;
;   unit - LUN of the pipe command, *not* the cache file.
;
;   buffer - the buffer to accept the data.  Data is read in
;            *unformatted*.
;
; Side Effects
;
;   The pipe is read as needed and the cache is populated.
;   The file pointer advances.
;
; MODIFICATIONS
;   Corrected error message, 21 Sep 2000, CM
;   Changed copyright notice, 21 Sep 2000, CM
;
;  $Id: fxpread.pro,v 1.2 2001/03/25 18:10:46 craigm Exp $
;
;-
; Copyright (C) 1999-2000, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
; -

pro FXPREAD, UNIT, BUFFER, transfer_count = TC
  compile_opt idl2

  @fxpcommn

  SZ = size(BUFFER)
  TYPECODE = SZ[SZ[0] + 1] < 12
  NBYTES = SZ[SZ[0] + 2] * BYTELENS(TYPECODE)

  if NBYTES lt 0 then begin
    TYPENAME = ['UNDEFINED', 'BYTE', 'INTEGER', 'LONG', $
      'FLOAT', 'DOUBLE', 'COMPLEX', 'STRING', $
      'STRUCTURE', 'DCOMPLEX', 'POINTER', $
      'OBJECT', 'UNKNOWN']
    message, 'ERROR: Cannot read unformatted ' + TYPENAME[TYPECODE] + ' data'
    RETURN
  endif

  if POINTER(UNIT) + NBYTES gt CACHE_LEN(UNIT) then begin
    FXPBUFFR, UNIT, POINTER(UNIT) + NBYTES
  endif

  readu, CACHE_UNIT(UNIT), BUFFER, transfer_count = TC
  point_lun, -CACHE_UNIT(UNIT), P
  POINTER(UNIT) = P

  RETURN
end