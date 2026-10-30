;+
; NAME:
;   CMSVWRITE
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;
; PURPOSE:
;   Write a single variable to an open SAVE file
;
; CALLING SEQUENCE:
;
;   CMSVWRITE, UNIT, DATA [ , NAME=NAME, COMPATIBILITY=COMPAT ]
;
; DESCRIPTION:
;
;   CMSVWRITE writes a single IDL variable to an open IDL SAVE file.
;   The file should already have been opened for writing as a normal
;   file using OPENW or OPENU.
;
;   CMSVWRITE is a simplified version of the CMSVLIB package, and as
;   such is not capable of writing heap data (pointers) or object
;   data, or structures that contain them.  Strings, structures, and
;   all array types are supported.
;
;   This procedure is part of the CMSVLIB SAVE library for IDL by
;   Craig Markwardt.  You must have the full CMSVLIB core package
;   installed in order for this procedure to function properly.
;
; ==================================================================
;   Research Systems, Inc. has issued a separate license intended
;   to resolve any potential conflict between this software and the
;   IDL End User License Agreement. The text of that license
;   can be found in the file LICENSE.RSI, included with this
;   software library.
; ==================================================================
;
; INPUTS:
;
;   UNIT - the open file unit.
;
;   DATA - the data to be written.
;
; KEYWORDS:
;
;   NAME - the optional name of the variable to be written (must be a
;          valid variable name).
;          Default: CMSVWRITE automatically creates a valid name.
;
;   COMPATIBILITY - a string, which describes the format to be used in
;          the output file.  Possible values are:
;
;                  'IDL4' - format of IDL version 4;
;                  'IDL5' - format of IDL versions 5.0-5.3;
;                  'IDL6' - not supported yet, for versions 5.4-above;
;                  'RIVAL1' - same as 'IDL5'
;           Note that files written in IDL5 format may still be
;           readable by IDL v.4.
;           Default: 'IDL5'
;
;   NO_END - a save file must terminate with an "end" record.  By
;            default, CMSVWRITE will append such a record after the
;            variable is written, and then rewind the file pointer.
;            The end record must be written after the last variable,
;            but is optional otherwise.  Set this keyword to disable
;            writing the end record (for performance reasons).
;
;   QUIET - if set, error messages are not printed.
;           Default: an error causes errors to be printed with MESSAGE
;
;   STATUS - upon return, this keyword will contain 1 for success and
;            0 for failure.
;
;   ERRMSG - upon return with a failure, this keyword will contain the
;            error condition as a string.
;
; EXAMPLE:
;
;   Write variables A, B, C and D to a file.
;
;   openw, 50, 'test.sav'       ;; Add /STREAM under VMS !
;   cmsvwrite, 50, a, name='a'
;   cmsvwrite, 50, b, name='b'
;   cmsvwrite, 50, c, name='c'
;   close, 50
;
; SEE ALSO:
;
;   CMSVREAD, CMRESTORE, CMSAVE, SAVE, CMSVLIB
;
; MODIFICATION HISTORY:
;   Written and documented, 11 Jan 2001, CM
;   Make version checks with correct precision, 19 Jul 2001, CM
;   Added notification about RSI License, 13 May 2002, CM
;
; $Id: cmsvwrite.pro,v 1.11 2002/05/13 06:41:10 craigm Exp $
;
;-
; Copyright (C) 2001, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
; -
pro cmsvwrite, unit0, data, name = name0, compat = compat0, no_end = noend, $
  quiet = quiet, status = status, errmsg = errmsg
  compile_opt idl2

  status = 0
  catch, catcherr
  if catcherr eq 0 then lib = cmsvlib(/query) else lib = 0
  catch, /cancel
  if lib eq 0 then begin
    errmsg = 'ERROR: The CMSVLIB library must be in your IDL path.'
    if keyword_set(quiet) then return else message, errmsg
  endif

  if n_elements(compat0) eq 0 then compat = 'IDL5' $
  else compat = strtrim(compat0, 2)

  if compat ne 'IDL4' and compat ne 'IDL5' and compat ne 'IDL6' and $
    compat ne 'RIVAL1' then begin
    errmsg = 'ERROR: unrecognized COMPAT value'
  endif

  ; ; Do type checking on the data
  if n_elements(data) eq 0 then begin
    errmsg = 'ERROR: DATA must be defined'
    if keyword_set(quiet) then return else message, errmsg
  endif

  if double(!version.release) gt 4d then begin
    has_objects = 0
    cmsv_ptrsum, null, /null
    cmsv_ptrsum, data, pheap, has_objects = has_objects

    if n_elements(pheap) gt 1 then goto, ptr_error
    if pheap[0] ne null or has_objects then begin
      ptr_error:
      errmsg = 'ERROR: CMSVWRITE cannot save pointers or objects'
      if keyword_set(quiet) then return else message, errmsg
    endif
  endif

  if n_elements(unit0) eq 0 then begin
    errmsg = 'ERROR: UNIT is not defined'
    if keyword_set(quiet) then return else message, errmsg
  endif
  unit = floor(unit0[0])

  stat = fstat(unit)
  if stat.write eq 0 or stat.open eq 0 then begin
    errmsg = 'ERROR: UNIT is not open for writing'
    if keyword_set(quiet) then return else message, errmsg
  endif

  ; ; We are at the beginning of the file, make sure this is a proper
  ; ; IDL save file.
  if stat.cur_ptr eq 0 then begin
    ; ; Open the file
    cmsv_open, unit, 'FILENAME', off0, access = 'W', /reopen, $
      status = status, errmsg = errmsg
    if status eq 0 then begin
      if keyword_set(quiet) then return else message, errmsg
    endif

    ; ; Write the opening records which establish the version number
    ; ; and time stamp.
    pp = 0l
    cmsv_wrec, block, pp, block_name = 'TIMESTAMP', offset = off0, $
      status = status, errmsg = errmsg
    if (status ne 0) and (compat ne 'IDL4') then $
      cmsv_wrec, block, pp, block_name = 'VERSION', offset = off0, $
      compat = compat, status = status, errmsg = errmsg
    if status eq 0 then begin
      if keyword_set(quiet) then return else message, errmsg
    endif

    ; ; Write the block
    writeu, unit, block[0 : pp - 1]
  endif
  point_lun, -unit, off0

  ; ; Construct a name for this variable, if one wasn't provided
  if n_elements(name0) eq 0 then begin
    default_name:
    name = string(off0, format = '("CMSV_",Z8.8)')
  endif else begin
    name = strupcase(strtrim(name0[0], 2))
    if name eq '' then goto, default_name
  endelse

  pp = 0l
  block = 0
  dummy = temporary(block)

  ; ; Make a block and write out the variable type and data.  The
  ; ; following procedure is required.  First, we create a new block,
  ; ; which establishes the header and record type (VARIABLE).  Second
  ; ; we write out the variable type.  Finally we write the actual
  ; ; data, and call the "finblock" procedure which rewrites the header
  ; ; record to have the correct "next" pointer.
  ; ;
  ; ; Note that all these procedures write to a block in memory first,
  ; ; and only at the end is the block written to disk.
  sz = size(data)
  off1 = pp
  cmsv_wrec, block, pp, data, name, block_name = 'VARIABLE', $
    /init, offset = off0, $
    status = status, errmsg = errmsg
  if status eq 0 then begin
    if keyword_set(quiet) then return else message, errmsg
  endif

  ; ; Now write out an end marker block, which will serve to terminate
  ; ; the file in case the user closes it.  We will rewind to this
  ; ; point in case new writes occur, so they will overwrite the end
  ; ; block.

  eoff = off0 + pp
  if not keyword_set(noend) then begin
    cmsv_wrec, block, pp, block_name = 'END_MARKER', offset = off0, $
      status = status, errmsg = errmsg
    if status eq 0 then begin
      if keyword_set(quiet) then return else message, errmsg
    endif
  endif

  writeu, unit, block[0 : pp - 1]
  block = 0

  ; ; Now rewind the file pointer so that it points at the end marker.
  ; ; Any new writes will overwrite the marker.
  if not keyword_set(noend) then point_lun, unit, eoff

  status = 1
  return
end