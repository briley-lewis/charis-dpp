;+
; NAME:
;   HELPFORM
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;
; PURPOSE:
;   Generate a descriptive string in IDL HELP format
;
; CALLING SEQUENCE:
;   STRINGS = HELPFORM(NAME, VALUE, [/SHORTFORM,] [/SINGLE,] [WIDTH=width])
;
; DESCRIPTION:
;
;   The HELPFORM function converts an IDL data value into a
;   representation very similar to the format produced by the built-in
;   command HELP.  Programmers can thus present data types and values
;   to users in a format they are familiar with.
;
;   For example, if the variable A is defined in the following manner,
;   and HELP is called, then the following transcript will result:
;
;     IDL> a = [1,2]
;     IDL> help, a
;     A               INT       = Array[2]
;
;   The same result can be achieved with the HELPFORM function:
;
;     IDL> print, helpform('A', a)
;     A               INT       = Array[2]
;
;   The benefit is that the output of HELPFORM is a string that can be
;   outputted or reformatted.  This capability is not available in all
;   versions of IDL.
;
;   HELPFORM actually produces *two* forms of output.  The above
;   output is considered the "long" form, as it appears in the IDL
;   HELP command, and is the default.  A "short" form can also be
;   produced, and is very similar to the information that appears in
;   certain terse IDL error messages.  It is activated by setting the
;   SHORTFORM keyword.
;
;   If the variable name is too long, the HELPFORM may be forced to be
;   two lines long to have consistent formatting.  In that case a
;   two-element string is returned.  If a single line is desired, use
;   the SINGLE keyword, but this comes at the expense of consistent
;   output formatting.
;
; INPUTS:
;
;   NAME - A scalar string containing the name of the IDL variable.
;          An empty string is permitted.  The name is ignored if the
;          SHORTFORM keyword is set.
;
;   VALUE - Any IDL value to be examined.  VALUE is optional if the
;           SIZE keyword is passed and uniquely describes the data.
;           VALUE should be passed for scalars and structures, since
;           the help form for these values requires additional
;           information beyond the SIZE.
;
; KEYWORDS:
;
;   SIZE - the IDL SIZE descriptor for the value to be printed.
;          Default: information is taken from VALUE.
;
;   SINGLE - if set, then output which would normally
;            appear on two lines for consistent formatting, appears on
;            one single line instead.
;
;   FULL_STRUCT - if set, then a detailed output is printed for
;                 structures, similar to HELP, VALUE, /STRUCTURE.
;

;   RECURSIVE_STRUCT - if both this keyword and FULL_STRUCT are set,
;                      and if VALUE itself has sub-structures, then
;                      print the full contents of those sub-structures
;                      as well.  The contents will be slightly indented.
;
;   SHORTFORM - set this keyword for a shorter output format that can
;               be used in error messages.
;
;   WIDTH - the width of the terminal in characters (used for
;           formatting).
;           Default: 80
;
; RETURNS:
;
;   An array of strings containing the HELPFORM output, which may have
;   more than one element depending on the length of NAME, SHORTFORM
;   and SINGLE.  The helpforms of pointer- and object-typed values
;   does not include the sequence number, but are otherwise correct.
;
; EXAMPLE:
;
;     IDL> print, helpform('A', size=[1,2,1,2])
;     A               BYTE      = Array[2]
;     ;; Do not pass VALUE and instead use SIZE to specify the type
;
;     IDL> print, helpform('A', size=[1,2,1,2], /shortform)
;     BYTE     (Array[2])
;     ;; Compare to the short form, which is meant to be placed in
;     ;; error messages
;
;     IDL> print, helpform('fjsldkfjsldfkjslkdfjslkdfjslkdfjsldkfjk',a)
;     fjsldkfjsldfkjslkdfjslkdfjslkdfjsldkfjk
;                      INT       = Array[2]
;     IDL> print, helpform('fjsldkfjsldfkjslkdfjslkdfjslkdfjsldkfjk',a,/single)
;     fjsldkfjsldfkjslkdfjslkdfjslkdfjsldkfjk INT       = Array[2]
;     ;; Compare the long and short forms
;
;
; SEE ALSO:
;
;   INPUTFORM, HELP
;
; MODIFICATION HISTORY:
;   Written, CM, 13 May 2000
;   Documented, 04 Jul 2000
;   Improved output for objects, CM, 11 Jan 2001
;   Added support for full structure output, CM 08 Feb 2001
;   Added forward_function declaration for safety, CM 08 Apr 2001
;   Print more info about POINTER type, CM 08 Apr 2001
;   Add the RECURSIVE_STRUCT keyword, CM 04 Jan 2009
;
; $Id: helpform.pro,v 1.6 2009/01/04 09:18:18 craigm Exp $
;
;-
; Copyright (C) 2000-2001, 2009, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
; -
forward_function helpform
function helpform, name0, value, size = sz, single = single, shortform = short, $
  width = width0, structure_name = stname, tagform = tagform, $
  full_struct = struct, recursive_struct = recstruct
  compile_opt idl2

  ; ; Names of all the known IDL types, as of IDL 5.2
  typenames = ['UNDEFINED', 'BYTE', 'INT', 'LONG', 'FLOAT', 'DOUBLE', $
    'COMPLEX', 'STRING', 'STRUCT', 'DCOMPLEX', 'POINTER', $
    'OBJREF', 'UINT', 'ULONG', $
    'LONG64', 'ULONG64', 'UNKNOWN']
  blanks = string(replicate(32b, 80))
  if n_elements(sz) lt 3 then sz = size(value)
  tp = sz[sz[0] + 1] < 16

  if n_elements(name0) eq 0 then name0 = ''
  name = strtrim(name0[0], 2)

  nlen = 15 ; ; Length of name
  tlen = 9 ; ; Length of type name

  if n_elements(width0) eq 0 then width0 = 80
  width = floor(width0[0])

  ; ; ================================ STRUCTURES
  if tp eq 8 and keyword_set(struct) then begin
    sz1 = size(value)
    if sz1[sz1[0] + 1] ne 8 then goto, not_struct
    nt = n_tags(value)
    len = n_tags(value, /length)
    tn = tag_names(value)
    sn = tag_names(value, /structure_name)

    if sn eq '' then sn = '<Anonymous>'

    a = string(sn, nt, len, $
      format = '("** Structure ",A0,", ",I0," tags, length=",I0,":")')

    for i = 0, nt - 1 do begin
      stri = helpform(tn[i], value[0].(i), /tagform)
      szi = size(value[0].(i))
      tpi = szi[szi[0] + 1]

      a = [a, '   ' + stri]
      ; ; Recursive structures
      if keyword_set(recstruct) and tpi eq 8 then begin
        stri = helpform(tn[i], value[0].(i), $
          /full_struct, /recursive_struct)
        a = [a, '     ' + stri]
      endif
    endfor

    return, a
  endif
  not_struct:

  if not keyword_set(short) then begin
    ; ; Pad the name out, or else put the name on a line by itself
    if strlen(name) gt nlen then begin
      if keyword_set(single) then begin
        a1 = name + ' '
      endif else begin
        a0 = name
        a1 = strmid(blanks, 0, nlen) + ' '
      endelse
    endif else begin
      a1 = strmid(name + blanks, 0, nlen) + ' '
    endelse

    a1 = a1 + strmid(typenames[tp] + blanks, 0, tlen)
    if not keyword_set(tagform) then $
      a1 = a1 + ' = '
  endif else begin
    a1 = strmid(typenames[tp] + blanks, 0, tlen)
  endelse

  ndims = sz[0]
  if ndims gt 0 then begin
    ; ; It is an array, compose the dimensions
    dims = sz[1 : ndims]
    v = 'Array['
    for i = 0l, ndims - 1 do begin
      v = v + strtrim(dims[i], 2)
      if i lt ndims - 1 then v = v + ', '
    endfor
    v = v + ']'

    ; ; If it is a structure, add the structure name (structures are
    ; ; never scalars)
    if not keyword_set(short) and tp eq 8 then begin
      ; ; Protect against empty value
      if n_elements(stname) eq 0 then begin
        if n_elements(value) gt 0 then v0 = value[0] else v0 = {dummy: 0}
        sn = tag_names(v0, /structure_name)
        sn = sn[0]
      endif else begin
        sn = strtrim(stname[0], 2)
      endelse
      if sn eq '' then sn = '<Anonymous>'
      v = '-> ' + sn + ' ' + v
    endif
  endif else begin
    ; ; It is a scalar

    ; ; Protect against empty or vector value
    if n_elements(value) gt 0 then begin
      v0 = value[0]
    endif else begin
      if tp ne 10 and tp ne 11 then tp = 0
    endelse

    case tp < 16 of
      0: v = '<Undefined>'
      1: v = string(v0, format = '(I4)')
      7: begin
        w = (width - 35) > 5
        if strlen(v0) gt w then $
          v = '''' + strmid(v0, 0, w) + '''...' $
        else $
          v = '''' + v0 + ''''
      end
      10: begin
        sz = size(v0)
        if sz[sz[0] + 1] eq 10 then v = string(v0[0], /print) $
        else v = '<PtrHeapVar>'
      end
      11: begin
        if n_elements(stname) eq 0 then begin
          forward_function obj_class
          sz = size(v0)
          if sz[sz[0] + 1] eq 11 then sn = '(' + obj_class(v0) + ')' $
          else sn = ''
        endif else begin
          sn = '(' + strupcase(strtrim(stname[0], 2)) + ')'
        endelse
        v = '<ObjHeapVar' + sn + '>'
      end
      16: v = ''
      else: v = string(v0)
    endcase
  endelse

  if keyword_set(short) then return, a1 + '(' + v + ')'

  a1 = a1 + v
  if n_elements(a0) gt 0 then a1 = [a0, a1]

  return, a1
end