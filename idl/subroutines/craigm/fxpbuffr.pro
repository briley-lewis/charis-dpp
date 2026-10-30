;+
; NAME:
;   FXPBUFFR
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE:
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Internal routine to read data from the pipe and store it in a
;            cache file.
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
;   unit - LUN of the pipe, *not* the cache file.
;
;   newlen - the new desired length of the file.  After a successful
;            call to FXPBUFFR, the length of the cache file will be at
;            least this long.
;
; Side Effects
;
;  The pipe is read and the cache increases in size.
;
; MODIFICATION HISTORY:
;   Changed copyright notice, 21 Sep 2000, CM
;
;  $Id: fxpbuffr.pro,v 1.3 2007/09/01 23:05:31 craigm Exp $
;
;-
; Copyright (C) 1999-2000, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
; -

pro FXPBUFFR, UNIT, NEWLEN
  compile_opt idl2

  @fxfilter
  @fxpcommn

  CUNIT = CACHE_UNIT(UNIT)
  CLEN = CACHE_LEN(UNIT)

  ; ; Compute the size of the buffer, rounded to the nearest granularity
  ; ; size.
  NTOT = ceil(double(NEWLEN - CLEN) / BUFFER_GRAN) * BUFFER_GRAN

  ; ; Advance the cache file pointer to the end of the file
  if not EOF_REACHED(UNIT) then $
    point_lun, CUNIT, CLEN

  ; ; This strange ON_IOERROR layout is because IDL has problems
  ; ; resetting the error handler inside the loop.  Lots of data is
  ; ; lost, presumably because it is buffered.  Therefore, I moved the
  ; ; ON_IOERROR command outside of the loop.
  READING = 0
  on_ioerror, io_failed

  ; ; Read from pipe until we are done.  To avoid the "drinking from a
  ; ; fire-hose effect, we read in 32K chunks.

  while not EOF_REACHED(UNIT) and NTOT gt 0 do begin
    ; ; Don't read more that BUFFER_MAX at a time, in case we fill
    ; ; the memory while seeking to the end of the pipe file.
    NREAD = NTOT < BUFFER_MAX
    BUFFER = bytarr(NREAD)

    ; ; You may ask, why this funny treatment of I/O errors here.
    ; ; Well, pipes don't have a well defined endpoint.  A read on
    ; ; the pipe with a fixed buffer size might read past the end of
    ; ; the pipe data stream.  We have to catch that error here since
    ; ; an EOF is not fatal.  See above for the initializing
    ; ; ON_IOERROR.  When such an error does occur, processing
    ; ; proceeds to the end of the loop, and if READING EQ 1, then
    ; ; control returns to READ_FAILED.  Ugghh.  But that's the way
    ; ; IDL (mis-)handles I/O.
    CC = 0
    READING = 1
    readu, UNIT, BUFFER, transfer_count = CC

    ; ; IO Operation could have failed.  If the read operation
    ; ; failed, that may actually be a signal that the pipe has
    ; ; finished.  There may still be good data to extract from the
    ; ; most recent READU command, so we may restart.
    if 0 then begin
      io_failed:
      if not READING then goto, read_done
    endif

    read_failed:
    READING = 0
    ; ; If the transfer count is zero, that is our clue that the pipe
    ; ; has ended.
    if CC eq 0 then begin
      ; ; If there really was a failure here, and not enough bytes
      ; ; were read, then that will be picked up when we read from
      ; ; the cache file.
      CC = fstat(UNIT)
      CC = CC.transfer_count
      if CC eq 0 then EOF_REACHED(UNIT) = 1
    endif

    ; ; Write the buffer out to the cache file ...
    if CC gt 0 then writeu, CUNIT, BUFFER[0 : CC - 1]

    ; ; ... and note the new cache length.
    CACHE_LEN(UNIT) = CACHE_LEN(UNIT) + CC
    CACHE_MAX(UNIT) = CACHE_MAX(UNIT) + CC

    NTOT = NTOT - CC
  endwhile

  read_done:
  on_ioerror, null

  ; ; Finally, seek back to the point we started at, since that is what
  ; ; the caller is expecting.
  point_lun, CUNIT, POINTER(UNIT)

  RETURN
end