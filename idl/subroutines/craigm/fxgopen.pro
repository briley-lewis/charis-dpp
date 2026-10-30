;+
; NAME:
;   FXGOPEN
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE:
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Open generic resource as a seekable file.
;
; MAJOR TOPICS:
;   File I/O, Pipes, URLs, FITS
;
; CALLING SEQUENCE:
;   FXGOPEN, UNIT, RESOURCE, ACCESS=ACCESS, ERRMSG=ERRMSG
;
; DESCRIPTION:
;
;   FXGOPEN opens a generic "resource" for reading or writing.  A
;   "resource" can be a file or a Unix pipe, or a standard network
;   URL for the http, https, or ftp protocols.  Networked URLs are
;   handled using the Unix command-line program called 'curl'.
;
;   Readable resources are fully random access.  You are permitted to
;   perform seek operations on both files and streams such as Unix
;   pipes.  In the case of a stream, the stream is read upon demand
;   and saved to an on-disk cache.
;
;   FXGOPEN also automatically recognizes some standard Unix file
;   extensions and operates on them.  For example, files ending with
;   '.gz' are recognized as being compressed with gzip, and are passed
;   through gzcat to uncompress them.  You can display existing
;   filename extension mappings and add new ones using the FXMAKEMAP
;   procedure.  This feature also worked with files retrieved over the
;   network, as long as the processing command declared with FXMAKEMAP
;   is able to accept '-' to indicate the data is supplied via
;   standard input.
;
;   The UNIT number is allocated using GET_LUN; however, the internal
;   implementation may allocate more LUNs.  Therefore you must use
;   FXGCLOSE to close the LUN and be sure that all resources are
;   deallocated.
;
;   You must use the specialized 'FXG' style functions to read, write
;   and seek on the resulting unit number:
;
;     FXGOPEN  - open resource
;     FXGCLOSE - close resource
;     FXGREAD  - read from resource
;     FXGWRITE - write to resource
;     FXGSEEK  - seek on resource (i.e., perform POINT_LUN)
;
;     FXGFILTERED - determine if resource is a normal file.
;
; INPUTS:
;
;   UNIT - FXGOPEN will return a LUN in this variable.  It should be
;          subsequently read and written with FXGREAD, FXGWRITE, and
;          closed with FXGCLOSE.
;
;   RESOURCE - a string, describing the resource to be opened.
;              FXGOPEN will automatically determine how to open it
;              according to:
;
;              * If a filename the suffix may be mapped according to
;              FXMAKEMAP.  In that case the appropriate pipe command
;              is opened as a Unix pipe with FXPOPENR.
;
;              * If a string beginning with "|" then the remaining
;              part of the string is interpretted as a Unix pipe
;              command, to be opened with FXPOPENR.
;
;              * If a URL (uniform resource locator), then it is
;              accessed.  Currently supported protocols are:
;
;                file - a local file
;                http - a file served by a web (HTTP) server
;                ftp - a file served an FTP server
;
;              I would like to add some sort of in-memory files,
;              probably with a "mem" protocol.
;
;
; KEYWORD PARAMETERS:
;
;   ACCESS - a string, set to the access privileges of the resource.
;            Possible values are:
;
;              'R'  - read-only
;              'W'  - write/create
;              'RW' - write/update
;
;            Not all protocols support writing (for example, none of
;            the "pipe" or network protocols supports writing).
;            DEFAULT: 'R'
;
;   ERRMSG - If a named variable is passed with this keyword, an error
;            message is returned: the empty string indicates success;
;            a non-empty string indicates failure.  If a named
;            variable is not passed, and the ERROR keyword is not
;            used, then execution is stopped upon an error.
;
;   ERROR - If a named variable is passed with this keyword, the error
;           status is returned: a zero indicates success; non-zero
;           indicates failure.  If a named variable is not passed, and
;           the ERRMSG keyword is not used, then execution is stopped
;           upon an error.
;
;   SUFFIX - Force a particular file type by specifying the suffix.
;            Default is to extract the suffix from the file name
;            itself.
;
; EXAMPLE:
;
;   fxgopen, unit, 'myfile.gz', errmsg=errmsg
;   if errmsg NE '' then do_error_message
;   bb = bytarr(1000)  ;; Read 1000 bytes
;   fxgread, unit, bb
;   fxgclose, unit
;
;   This example opens the file myfile.gz using FXGOPEN.  It is
;   automatically gunzip'ed on demand as the request for a 1000-byte
;   read is made.
;
; MODIFICATION HISTORY:
;   Written, 1999, CM
;   Documented, 02 Oct 1999, CM
;   Added correct ERROR keyword behavior, 04 Oct 1999, CM
;   Changed copyright notice, 21 Sep 2000, CM
;   Modified to use ARG_PRESENT for ERRMSG and ERROR, 21 Sep 2000, CM
;   Added SUFFIX keyword, 31 Oct 2000, CM
;   Added the HTTP and FTP protocols using curl, 22 Oct 2006, CM
;
; TODO:
;   * Make more windows friendly
;
;  $Id: fxgopen.pro,v 1.5 2009/02/12 02:32:50 craigm Exp $
;
;-
; Copyright (C) 1999-2000,2006 Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
; -

; ; Utility routine: open a network resource using the 'curl' command
pro fxgopen_curl, unit, resource, suffix, errmsg = errmsg, error = error, _extra = extra
  compile_opt idl2
  ; ; The curl command automatically redirects to stdout
  cmd = string(resource[0], format = '("curl -s ''",A0,"''")')
  @fxfilter
  wh = where(suffix eq filters(0, *), ct)
  ; ; Handle the case where the remote file is gzipped, compressed, etc
  ; ; XXX: this assumes that the command can take '-' to mean 'stdin'
  if ct gt 0 then cmd = cmd + ' | ' + string('-', format = filters(1, wh[0]))
  FXPOPENR, unit, cmd, errmsg = errmsg, error = error, _extra = extra
  return
end

pro FXGOPEN, UNIT, RESOURCE, access = ACCESS0, errmsg = errmsg, $
  error = error, suffix = suffix0, _extra = extra
  compile_opt idl2

  on_error, 2
  error = -1
  errmsg = ''

  ; ; Default the parameters
  if n_elements(ACCESS0) eq 0 then ACCESS0 = 'R'
  ACCESS = strupcase(ACCESS0)
  if ACCESS ne 'R' and ACCESS ne 'W' and ACCESS ne 'RW' then begin
    MESSAGE = 'ERROR: ACCESS must be R, W, or RW.'
    goto, err_return
  endif

  ; ; Check that the resource is at least a string.
  sz = size(RESOURCE)
  if sz[sz[0] + 1] ne 7 then begin
    MESSAGE = 'ERROR: RESOURCE must be a string.'
    goto, err_return
  endif

  ; ; Separate the protocol component of a URL
  len = strlen(RESOURCE)
  i = 0l
  while i lt len and strmid(RESOURCE, i, 1) ne ':' $
    and strmid(RESOURCE, i, 1) ne '/' do i = i + 1
  if i eq len or (i lt len and strmid(RESOURCE, i, 1) eq '/') then begin
    protocol = 'file'
    location = RESOURCE
  endif else begin
    if i eq 0 or i eq len - 1 then begin
      MESSAGE = 'ERROR: incorrect resource name format'
      goto, err_return
    endif
    protocol = strmid(RESOURCE, 0, i)
    location = strmid(RESOURCE, i + 1, strlen(RESOURCE) - i - 1)
  endelse

  ; ; An ode to DOS: single-letter protocols are probably disk drives
  if strlen(protocol) eq 1 then begin
    protocol = 'file'
    location = RESOURCE
  endif

  ; ; Separate the server component
  len = strlen(location)
  i = 0l
  while i lt len and strmid(location, i, 1) eq '/' do i = i + 1
  if i eq 0 or i eq 1 then begin ; ; No slash, or a single slash -- a local file

    if i eq len then begin
      MESSAGE = 'ERROR: incorrect resource name format'
      goto, err_return
    endif
    server = ''
    path = location
  endif else if i eq 3 then begin ; ; Three slashes -- a local file
    server = ''
    path = strmid(location, 2, len - 2)
  endif else if i gt 3 then begin ; ; Too many slashes
    MESSAGE = 'ERROR: incorrect resource name format'
    goto, err_return
  endif else begin ; ; Format proto://server[/path]
    path = strmid(location, 2, len - 2)
    slash = strpos(path, '/')
    if slash eq -1 then begin ; ; No path
      server = path
      path = ''
    endif else begin ; ; Server and path
      server = strmid(path, 0, slash)
      path = strmid(path, slash, strlen(path) - slash)
    endelse
  endelse

  ; ; Determine the suffix of the path
  components = STR_SEP(path, '.')
  len = n_elements(components)
  if len gt 1 then suffix = components[len - 1] else suffix = ''
  if n_elements(suffix0) gt 0 then $
    suffix = strtrim(suffix0[0], 2)

  ; ; Find out if this is a pipe
  if strmid(path, 0, 1) eq '|' then begin
    if ACCESS ne 'R' then begin
      MESSAGE = 'ERROR: pipes may only be opened with READ access.'
      goto, err_return
    endif
    FXPOPENR, UNIT, path, errmsg = errmsg, error = error, _extra = extra
    return
  endif

  @fxfilter
  case strlowcase(protocol) of
    ; ; FILE access is the only supported protocol currently.
    'file': begin
      wh = where(suffix eq filters(0, *), ct)
      if ct gt 0 then begin ; ; A filtered file must spawn a pipe

        ; ; This file suffix is associated with a PIPE
        if ACCESS ne 'R' then begin
          MESSAGE = 'ERROR: pipes may only be opened with READ access.'
          goto, err_return
        endif

        ; ; Check that the file itself is read-openable.
        openr, UNIT, path, /get_lun, error = error
        if error ne 0 then goto, open_error
        free_lun, UNIT

        ; ; If it is, then open a pipe on it.
        fmt = filters(1, wh[0])
        flags = filters(2, wh[0])
        flags = strtrim(strcompress(strupcase(flags)), 2)
        compress = 0
        if flags eq '' then begin
          cmd = string(path, format = fmt)
        endif else begin
          case 1 of
            (strpos(flags, 'COMPRESS') ge 0): compress = 1
          endcase

          cmd = path
        endelse

        FXPOPENR, UNIT, cmd, compress = compress, $
          errmsg = errmsg, error = error, _extra = extra
        return
      endif else begin
        ; ; General file access is achieved through trusty
        ; ; OPEN[RWU]

        case ACCESS of
          'R': openr, UNIT, path, /block, /get_lun, error = error
          'W': openw, UNIT, path, /block, /get_lun, error = error
          'RW': openu, UNIT, path, /block, /get_lun, error = error
        end
        if error ne 0 then begin
          open_error:
          ; ; Deal with the error condition
          MESSAGE = 'ERROR: could not open file "' + path + '"'
          goto, err_return
        endif

        ; ; Make sure the FXFILTER entry is zeroed.  We don't
        ; ; want trouble!
        filterflag(UNIT) = 0
        seek_cmd(UNIT) = ''
        read_cmd(UNIT) = ''
        write_cmd(UNIT) = ''
        close_cmd(UNIT) = ''
        return
      endelse
    end
    'http': fxgopen_curl, UNIT, RESOURCE, suffix, errmsg = errmsg, error = error, _extra = extra
    'https': fxgopen_curl, UNIT, RESOURCE, suffix, errmsg = errmsg, error = error, _extra = extra
    'ftp': fxgopen_curl, UNIT, RESOURCE, suffix, errmsg = errmsg, error = error, _extra = extra
    else: begin
      ; ; Sorry... we need more protocols here, but probably with
      ; ; an external program such as CURL
      MESSAGE = 'ERROR: protocol "' + protocol + '" is not supported'
      goto, err_return
    end
  endcase

  return

  err_return:
  forward_function arg_present ; ; For IDL versions before 5

  if arg_present(errmsg) or arg_present(error) then begin
    errmsg = MESSAGE
    return
  endif

  if double(!version.release) lt 5 then begin
    if n_elements(errmsg) ne 0 then begin
      errmsg = MESSAGE
      return
    endif
  endif
  message, MESSAGE
end