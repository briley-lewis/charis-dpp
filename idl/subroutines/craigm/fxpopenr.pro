;+
; NAME:
;   FXPOPENR
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE:
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Internal routine to open a Unix pipe command for read access.
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
; Usage: FXPOPENR, UNIT, COMMAND, ERRMSG=ERRMSG
;
; PARAMETERS
;
;   unit - FXPOPENR returns the pipe LUN, created by GET_LUN in this
;          parameter.  The LUN should not be "pre-opened".
;          Unformatted reads on this LUN should be performed with
;          FXPREAD.
;
;   command - a scalar string, the pipe command to execute.  The
;             standard output of the command is redirected into UNIT.
;             Standard error is not redirected.
;
;             A failure of the command can only be discovered upon
;            trying to read from the LUN with FXPREAD.
;
; Keywords
;
;   errmsg - If set to defined value upon input, an error message is
;            returned upon output.  If no error occurs then ERRMSG is
;            not changed.  If an error occurs and ERRMSG is not
;            defined, then FXPOPENR issues a MESSAGE.
;
; Side Effects
;
;   The pipe command is opened with SPAWN, and an additional cache file
;   is opened with read/write access.
;
;   The FXFILTER family of commons is updated.
;
; MODIFICATION HISTORY:
;   Changed copyright notice, 21 Sep 2000, CM
;   Added the OPEN,/DELETE keyword, so that tmp-files are
;     automatically deleted when closed, 18 Feb 2006, CM
;   Added quotation marks to the list of characters which are
;     protected, while making a tmpfile name, 22 Oct 2006, CM
;
;  $Id: fxpopenr.pro,v 1.7 2009/02/12 02:32:50 craigm Exp $
;
;-
; Copyright (C) 1999-2000,2006 Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
; -

; Utility program to protect the pipe command
pro FXPOPENR_WRAP_CMD, CMD, SHELL
  compile_opt idl2

  ; SHELL = '/bin/sh -c '

  while strmid(CMD, 0, 1) eq '|' or strmid(CMD, 0, 1) eq ' ' do $
    CMD = strmid(CMD, 1, strlen(CMD) - 1)
  ; CMD = SHELL + '"' + CMD + ' 2>/dev/null"'
  CMD = CMD + ' 2>/dev/null'

  RETURN
end

; Utility program to generate a name for the cache file.  It is
; uniquely generated based on the time, the command string, and the
; current call number.
;
function FXPOPENR_TMPNAME, CMD
  compile_opt idl2
  @fxfilter
  common FXPOPEN_TMPNAME, RANDOM_SEED, SEQ_COUNTER
  if n_elements(RANDOM_SEED) eq 0 then begin
    RANDOM_VAL = long(systime(1))
    RANDOM_SEED = long(randomu(RANDOM_VAL) * double(ishft(1l, 31)))
    SEQ_COUNTER = 0l
  endif

  ; ; Take the first fifteen and  characters of the command
  TMPNAME = strcompress(CMD, /remove_all)

  ; ; Build a unique hash name based on the command, the current time,
  ; ; and a session-specific seed.  Possible problem here: if several
  ; ; sessions are started at the same time with the same command, and
  ; ; the commands are executed at the same second, then the temporary
  ; ; name will be the same.  I judge the likelihood of all of these
  ; ; events to be small.

  B = byte(CMD)
  N = n_elements(B)
  ; ; Construct a semi-unique hash value for the command string
  HASH = 0l
  for I = 0l, N - 1 do HASH = ishft(HASH, 2) xor B[I]
  HASH = HASH xor long(systime(1)) xor RANDOM_SEED xor ishft(SEQ_COUNTER, 16)
  SEQ_COUNTER = SEQ_COUNTER + 1

  if strlen(TMPNAME) gt 20 then begin
    TMPNAME = strmid(TMPNAME, 0, 15) + strmid(TMPNAME, N - 6, 5)
    N = 20l
  endif
  NEWNAME = ''
  ; ; Strip away any non-alpha characters
  for I = 0l, N - 1 do begin
    CC = strmid(TMPNAME, I, 1)
    if not (CC eq ' ' or CC eq '>' or CC eq '&' or CC eq '|' or $
      CC eq '/' or CC eq '*' or CC eq '?' or CC eq '<' or $
      CC eq '\' or $
      CC eq '"' or CC eq '''') then $
      NEWNAME = NEWNAME + CC
  endfor
  if NEWNAME eq '' then NEWNAME = 'fxp'

  RETURN, SCRATCH_DIR + NEWNAME + string(abs(HASH), format = '(Z8.8)')
end

; ; Main entry
pro FXPOPENR, UNIT, CMD, errmsg = ERRMSG, error = error, compress = compress
  compile_opt idl2

  ; ; Access the general FXFILTER family of commons, and the
  ; ; FXPIPE_COMMON, which has pipe-specific info.
  error = -1
  @fxfilter
  @fxpcommn

  if n_params() lt 2 then begin
    MESSAGE = 'Syntax:  FXPOPEN, UNIT, COMMAND'
    goto, err_return
  endif

  ; ; Initialize filter flags
  FFLAGS = 1l

  if not keyword_set(compress) then begin
    ; ; Sorry, useful pipes are only available under Unix.
    if strupcase(!version.os_family) ne 'UNIX' then begin
      MESSAGE = 'ERROR: FXPOPENR ONLY FUNCTIONS ON UNIX SYSTEMS.'
      goto, err_return
    endif

    ; ; --------- Begin pipe section
    ; ; Wrap the command to make sure it is safe
    NEWCMD = CMD
    FXPOPENR_WRAP_CMD, NEWCMD, SHELL

    ; ; Run the program
    OLDSHELL = getenv('SHELL')
    on_ioerror, spawn_failed
    if OLDSHELL ne '/bin/sh' then setenv, 'SHELL=/bin/sh'
    spawn, NEWCMD, unit = UNIT, pid = PID
    on_ioerror, null
    setenv, 'SHELL=' + OLDSHELL[0]

    ; ; Check for error conditions
    if UNIT lt 1l or UNIT gt 128l then begin
      spawn_failed:
      setenv, 'SHELL=' + OLDSHELL[0]
      MESSAGE = 'ERROR: SPAWN of "' + NEWCMD + '" FAILED'
      goto, err_return
    endif

    FFLAGS = FFLAGS or 2 ; ; This is a pipe
    ; ; ---- End pipe section
  endif else begin
    ; ; Compressed data - no PID
    PID = 0l

    openr, UNIT, CMD, /get_lun, /compress, error = error
    if error ne 0 then begin
      MESSAGE = 'ERROR: OPEN of compressed file "' + CMD + '" FAILED'
      goto, err_return
    endif

    ; ; FFLAGS (unchanged since it is not a pipe)
  endelse

  ; ; Prepare the FXFILTER dispatch table for function calls
  FILTERFLAG(UNIT) = 1 ; ; Flags: XXX will be updated below!
  SEEK_CMD(UNIT) = 'FXPSEEK'
  READ_CMD(UNIT) = 'FXPREAD'
  WRITE_CMD(UNIT) = '-' ; ; This pipe is not writable
  CLOSE_CMD(UNIT) = 'FXPCLOSE'

  ; ; Start filling in the FXPIPE_COMMON
  POINTER(UNIT) = 0l ; ; Start of pipe
  PROCESS_ID(UNIT) = PID ; ; Save process ID of pipe

  ; ; Build a unique cache name
  CACHE_FILENAME = FXPOPENR_TMPNAME(CMD)

  ; ; Open the output cache file, retrieving a LUN
  on_ioerror, open_error
  openw, CACHE, CACHE_FILENAME, /get_lun, /delete
  on_ioerror, null

  FFLAGS = FFLAGS or 4 ; ; On-disk backing store

  ; ; Error condition on the cache file
  if CACHE lt 1 or CACHE gt 128 then begin
    open_error:

    ; ; Reset to default behavior
    FILTERFLAG(UNIT) = 0
    SEEK_CMD(UNIT) = ''
    READ_CMD(UNIT) = ''
    WRITE_CMD(UNIT) = ''
    CLOSE_CMD(UNIT) = ''
    free_lun, UNIT

    MESSAGE = 'ERROR: Unable to open cache file ' + strtrim(CACHE_FILENAME, 2)
    goto, err_return
  endif

  ; ; Finish filling the pipe information
  CACHE_UNIT(UNIT) = CACHE
  CACHE_FILE(UNIT) = CACHE_FILENAME
  POINTER(UNIT) = 0 ; ; At beginning of pipe
  CACHE_LEN(UNIT) = 0 ; ; Currently no data cached
  CACHE_MAX(UNIT) = 0 ; ; Currently no storage allocated for cache
  EOF_REACHED(UNIT) = 0 ; ; Not known to be at end-of-file

  ; ; Update filter flags
  FILTERFLAG(UNIT) = FFLAGS

  good_return:
  error = 0
  RETURN

  err_return:
  if n_elements(ERRMSG) ne 0 then begin
    ERRMSG = MESSAGE
    RETURN
  endif else message, MESSAGE
  RETURN
end