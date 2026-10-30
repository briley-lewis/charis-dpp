;
;+
; NAME:
;       SHOWTEX
; PURPOSE:
;       Display TeX sequence translation table on current graphics device.
; CATEGORY:
;       text/strings
; CALLING SEQUENCE:
;       showtex
; INPUTS:
; KEYWORD PARAMETERS:
;       /HELP -- print out info on use of the function
;                and exit.
;       FONT  -- Set to 0 to use hardware font, -1 to use vector.
;                Note that the only hardware font supported is
;                Postscript
; OUTPUTS:
; COMMON BLOCKS:
; SIDE EFFECTS:
;       Plot is created.
; NOTES:
;       Hardware fonts are supported only for device PS (PostScript)
; EXAMPLE:
; MODIFICATION HISTORY:
;       $Id: showtex.pro,v 1.4 2004/06/15 17:25:54 mcraig Exp $
;       $Log: showtex.pro,v $
;       Revision 1.4  2004/06/15 17:25:54  mcraig
;       Fixed bug in regular expression, changed array notation to square brackets
;
;       Revision 1.3  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.2  1996/05/09 00:22:17  mcraig
;       Added error handling and updated built in help.
;
;       Revision 1.1  1996/02/08 18:55:12  mcraig
;       Initial revision
;
; RELEASE:
;       $Name: Rel_2_1_2 $
;
; COPYRIGHT:
;  Copyright (C) 1996 The Regents of the University of California, All
;  Rights Reserved.  Written by Matthew W. Craig.
;  See the file COPYRIGHT for restrictions on distrubting this code.
;  This code comes with absolutely NO warranty; see DISCLAIMER for details.
;-
;
pro Showtex, font = fnt, help = help
  compile_opt idl2

  ; Return to caller on error.
  on_error, 2

  ; Print help if needed.
  if keyword_set(help) then begin
    print, '    Display TeX sequence translation table on current graphics device.'
    print, '    showtex'
    print, '    Keywords:'
    print, '       /HELP       print this message and return'
    print, '       FONT        set to 0 to use hardware fonts for current device,'
    print, '                   -1 to use vector fonts (DEFAULT)'
    print, '    NOTES:  - The only hardware font supported is PostScript.'
    print, '            - The FONT keyword overrides the font selected in !p.font'
    return
  endif

  ; We begin by deciding on the font.  PostScript = 0 means use vector.
  PostScript = 0
  PlotTitle = 'Vector Fonts'
  if n_elements(fnt) eq 0 then begin ; get font from !p.font
    if !p.font ne -1 then begin ; User wants hardware font.
      PostScript = 1
      PlotTitle = 'PostScript Fonts'
    endif
  endif else begin ; get font from FONT keyword
    if fnt ne -1 then begin
      PostScript = 1
      PlotTitle = 'PostScript Fonts'
    endif
  endelse

  ; Bomb out if user wants hardware font for non-PostScript device.
  if (PostScript eq 1) and (strupcase(!d.name) ne 'PS') then begin
    ; Device isn't postscript
    ; and user wants hardware
    ; font.  Not good.
    print, 'Warning: No translation for device: ', !d.name
    return
  endif

  ; Set !P.font to value indicated by FONT keyword, saving surrent
  ; setting to reset at end.
  OldPFont = !p.font
  !p.font = PostScript - 1

  erase
  seq = Textoidl(/tex)
  DisplayString = seq + '  ' + Textoidl(seq)

  nseq = n_elements(seq)
  nrows = nseq / 5 + 1 ; Five sequences per row.
  dx = .9 / 5.
  dy = .9 / nrows
  y = .95
  xyouts, .5, y, PlotTitle, align = .5, /norm, size = 2.5
  count = 0
  for i = 1l, nrows do begin
    y = y - dy
    x = .1
    for j = 1, 5 do begin
      if (count lt nseq) then xyouts, x, y, DisplayString[count], align = .5, /norm
      count = count + 1
      x = x + dx
    endfor
  endfor

  ; Restore old !P.font.
  !p.font = OldPFont
end