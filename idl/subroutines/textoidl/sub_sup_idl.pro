;
;+
; NAME:
;       SUB_SUP_IDL
; PURPOSE:
;       Return the proper IDL font positioning command for TeX
;       sub/superscripts.
; CATEGORY:
; CALLING SEQUENCE:
;       fnt = sub_sup_idl( strn )
; INPUTS:
;       strn -- Either '^' or '_', the TeX super/subscript       in
;               characters
; KEYWORD PARAMETERS:
;       /FORCE_UD -- Set this to use !U/!D instead of !E/!I for
;                    sub/superscripts .
; OUTPUTS:
;       fnt -- Either '!U' or !E' for superscripts,              out
;              or '!D' or '!I' for subscripts.
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
; EXAMPLE:
; LIBRARY FUNCTIONS CALLED:
;
; MODIFICATION HISTORY:
;       $Id: sub_sup_idl.pro,v 1.1 1996/01/31 18:47:37 mcraig Exp $
;       $Log: sub_sup_idl.pro,v $
;       Revision 1.1  1996/01/31 18:47:37  mcraig
;       Initial revision
;
; RELEASE:
;       $Name: Rel_2_1_2 $
;-
function Sub_sup_idl, token, force_ud = force_ud
  compile_opt idl2

  if keyword_set(force_ud) then begin
    if (token eq '^') then return, '!U'
    if (token eq '_') then return, '!D'
    return, ''
  endif else begin
    if (token eq '^') then return, '!E'
    if (token eq '_') then return, '!I'
    return, ''
  endelse
end