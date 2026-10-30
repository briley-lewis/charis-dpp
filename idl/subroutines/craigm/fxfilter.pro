; $Id: fxfilter.pro,v 1.4 2009/02/12 02:32:50 craigm Exp $
; -
; Copyright (C) 1999-2000, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
; -

; COMMON FXFILTER
; * Defines mapping between normal IDL I/O and FILTER I/O
; FILTERFLAG - for each LUN, = 1*FILTERED + 2*PIPE + 4*DISK_STORE
; FILTERED - use FXG library? 1=yes 0=no
; PIPE     - is a pipe?       1=yes 0=no
; DISK_STORE - backing store on disk=1 mem=0
; SEEK_CMD   - for each LUN, procedure to execute when performing POINT_LUN
; READ_CMD   - for each LUN, procedure to execute when performing READU
; WRITE_CMD  - for each LUN, procedure to execute when performing WRITEU
; CLOSE_CMD  - for each LUN, procedure to execute when performing CLOSE
;
compile_opt idl2

common FXFILTER, FILTERFLAG, SEEK_CMD, READ_CMD, WRITE_CMD, CLOSE_CMD
if n_elements(FILTERFLAG) eq 0 then begin
  FILTERFLAG = intarr(256)
  SEEK_CMD = strarr(256)
  READ_CMD = strarr(256)
  WRITE_CMD = strarr(256)
  CLOSE_CMD = strarr(256)
endif
FXFILTER_MAX_LUN = 256

; COMMON FXFILTER_CONFIG
; * Defines general configuration of the filter package
; * Can be manipulated by FXMAKEMAP
; SCRATCH_DIR - directory where cache files are stored   (string)
; BUFFER_MAX  - maximum buffer size (in bytes) of a pipe read (long)
; BUFFER_GRAN - buffer size granularity (in bytes), should be a large
; power of 2, probably >= 512 (long)
; RM_COMMAND  - unix command to use to remove a file (string)
; CACHE_MAX   - maximum in-memory cache of a filtered file
;
common FXFILTER_CONFIG, SCRATCH_DIR, BUFFER_MAX, BUFFER_GRAN, $
  RM_COMMAND, MEMCACHE_MAX
if n_elements(SCRATCH_DIR) eq 0 then begin
  SCRATCH_DIR = getenv('IDL_TMPDIR')
  BUFFER_GRAN = 4096l
  BUFFER_MAX = 8l * BUFFER_GRAN
  RM_COMMAND = '/bin/rm'
  MEMCACHE_MAX = 10l * 1024l * 1024l
endif

; COMMON FXFILTER_FILTERS
; * Defines mapping between suffixes and commands used to read them.
; * Can be manipulated by FXMAKEMAP
; FILTERS - an array of pairs.  The first of the pair gives the
; filename suffix to be mapped (without leading '.'), and
; the second of the pair gives command to be executed when
; the suffix is encountered.  The command should be in the
; form of an IDL format statement which transfers the
; filename into the command.
;
common FXFILTER_FILTERS, FILTERS
if n_elements(FILTERS) eq 0 then begin
  ; SUFF     FORMAT_COMMAND     FLAGS
  FILTERS = [ $
    ['gz', '("gzip -dc ",A0)', 'COMPRESS'], $
    ['Z', '("zcat ",A0)', ''] $
    ]
endif