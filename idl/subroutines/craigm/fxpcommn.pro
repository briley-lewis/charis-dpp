;
; FXPIPE_COMMON - information related to pipes in the FXFILTER package.
;
; POINTER     - for each LUN, the current file pointer (a long integer)
; CACHE_UNIT  - for each pipe LUN, the LUN of the cache file.
; CACHE_LEN   - for each cache LUN, the number of cached bytes.
; CACHE_MAX   - for each cache LUN, the maximum size of the cache.
; CACHE_FILE  - for each pipe LUN, the name of the cache file. (a string)
; EOF_REACHED - for each pipe LUN, a flag indicating whether the end
; of the pipe has been reached.
; PROCESS_ID  - for each pipe LUN, the process ID of the pipe
; command.
; BYTELENS    - the length in bytes of each IDL data type.
;
; MODIFICATION HISTORY:
; Changed copyright notice, 21 Sep 2000, CM
;
; $Id: fxpcommn.pro,v 1.3 2007/09/01 23:01:17 craigm Exp $
;
; -
; Copyright (C) 1999-2000, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
; -
compile_opt idl2
common FXPIPE_COMMON, POINTER, CACHE_UNIT, CACHE_LEN, CACHE_MAX, $
  CACHE_FILE, BYTELENS, EOF_REACHED, PROCESS_ID

if n_elements(POINTER) eq 0 then begin
  POINTER = lonarr(256)
  CACHE_UNIT = lonarr(256)
  CACHE_LEN = lonarr(256)
  CACHE_MAX = lonarr(256)
  CACHE_FILE = strarr(256)
  EOF_REACHED = lonarr(256)
  PROCESS_ID = lonarr(256)

  ; ; This is the length in bytes of each IDL type
  ; 0   1   2   3   4   5   6   7   8   9  10  11  12
  BYTELENS = [-1l, 1, 2, 4, 4, 8, 16, -1, -1, 16, $
    - 1, -1, 2, 4, 8, 8, -1]
endif