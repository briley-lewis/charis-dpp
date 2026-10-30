pro sxaddpar_dl, Header, Name, Value, Comment, Location, before = before, $
  savecomment = savecom, after = after, format = format, pdu = pdu
  compile_opt idl2
  ;+
  ; NAME:
  ;       SXADDPAR
  ; PURPOSE:
  ;       Add or modify a parameter in a FITS header array.
  ;
  ; CALLING SEQUENCE:
  ;       SXADDPAR, Header, Name, Value, [ Comment,  Location, /SaveComment,
  ;                               BEFORE =, AFTER = , FORMAT= , /PDU]
  ;
  ; INPUTS:
  ;       Header = String array containing FITS or STSDAS header.    The
  ;               length of each element must be 80 characters.    If not
  ;               defined, then SXADDPAR will create an empty FITS header array.
  ;
  ;       Name = Name of parameter. If Name is already in the header the value
  ;               and possibly comment fields are modified.  Otherwise a new
  ;               record is added to the header.  If name is equal to 'COMMENT'
  ;               or 'HISTORY' or a blank string then the value will be added to
  ;               the record without replacement.  For these cases, the comment
  ;               parameter is ignored.
  ;
  ;       Value = Value for parameter.  The value expression must be of the
  ;               correct type, e.g. integer, floating or string.  String values
  ;                of 'T' or 'F' are considered logical values.
  ;
  ; OPTIONAL INPUT PARAMETERS:
  ;       Comment = String field.  The '/' is added by this routine.  Added
  ;               starting in position 31.    If not supplied, or set equal to
  ;               '', or /SAVECOMMENT is set, then the previous comment field is
  ;               retained (when found)
  ;
  ;       Location = Keyword string name.  The parameter will be placed before the
  ;               location of this keyword.    This parameter is identical to
  ;               the BEFORE keyword and is kept only for consistency with
  ;               earlier versions of SXADDPAR.
  ;
  ; OPTIONAL INPUT KEYWORD PARAMETERS:
  ;       BEFORE  = Keyword string name.  The parameter will be placed before the
  ;               location of this keyword.  For example, if BEFORE='HISTORY'
  ;               then the parameter will be placed before the first history
  ;               location.  This applies only when adding a new keyword;
  ;               keywords already in the header are kept in the same position.
  ;
  ;       AFTER   = Same as BEFORE, but the parameter will be placed after the
  ;               location of this keyword.  This keyword takes precedence over
  ;               BEFORE.
  ;
  ;       FORMAT  = Specifies FORTRAN-like format for parameter, e.g. "F7.3".  A
  ;               scalar string should be used.  For complex numbers the format
  ;               should be defined so that it can be applied separately to the
  ;               real and imaginary parts.  If not supplied then the default is
  ;               'G19.12' for double precision, and 'G14.7' for floating point.
  ;
  ;       /PDU    = specifies keyword is to be added to the primary data unit
  ;               header. If it already exists, it's current value is updated in
  ;               the current position and it is not moved.
  ;       /SAVECOMMENT = if set, then any existing comment is retained, i.e. the
  ;               COMMENT parameter only has effect if the keyword did not
  ;               previously exist in the header.
  ; OUTPUTS:
  ;       Header = updated FITS header array.
  ;
  ; EXAMPLE:
  ;       Add a keyword 'TELESCOP' with the value 'KPNO-4m' and comment 'Name
  ;       of Telescope' to an existing FITS header h.
  ;
  ;       IDL> sxaddpar, h, 'TELESCOPE','KPNO-4m','Name of Telescope'
  ; NOTES:
  ;       The functions SXADDPAR() and FXADDPAR() are nearly identical, with the
  ;       major difference being that FXADDPAR forces required FITS keywords
  ;       BITPIX, NAXISi, EXTEND, PCOUNT, GCOUNT to appear in the required order
  ;       in the header, and FXADDPAR supports the OGIP LongString convention.
  ;       There is no particular reason for having two nearly identical
  ;       procedures, but both are too widely used to drop either one.
  ;
  ;       All HISTORY records are inserted in order at the end of the header.
  ;
  ;       All COMMENT records are also inserted in order at the end of the header
  ;       header, but before the HISTORY records.  The BEFORE and AFTER keywords
  ;       can override this.
  ;
  ;       All records with no keyword (blank) are inserted in order at the end of
  ;       the header, but before the COMMENT and HISTORY records.  The BEFORE and
  ;       AFTER keywords can override this.

  ; RESTRICTIONS:
  ;       Warning -- Parameters and names are not checked
  ;               against valid FITS parameter names, values and types.
  ;
  ; MODIFICATION HISTORY:
  ;       DMS, RSI, July, 1983.
  ;       D. Lindler Oct. 86  Added longer string value capability
  ;       Converted to NEWIDL  D. Lindler April 90
  ;       Added Format keyword, J. Isensee, July, 1990
  ;       Added keywords BEFORE and AFTER. K. Venkatakrishna, May '92
  ;       Pad string values to at least 8 characters   W. Landsman  April 94
  ;       Aug 95: added /PDU option and changed routine to update last occurence
  ;               of an existing keyword (the one SXPAR reads) instead of the
  ;               first occurence.
  ;       Comment for string data can start after column 32 W. Landsman June 97
  ;       Make sure closing quote supplied with string value  W. Landsman  June 98
  ;       Converted to IDL V5.0    W. Landsman   June 98
  ;       Increase precision of default formatting of double precision floating
  ;               point values.   C. Gehman, JPL  September 1998
  ;       Mar 2000, D. Lindler, Modified to use capital E instead of lower case
  ;               e for exponential formats.
  ;       Apr 2000, Make user-supplied format upper-case  W. Landsman
  ;       Oct 2001, Treat COMMENT or blank string like HISTORY keyword W. Landsman
  ;       Jan 2002, Allow BEFORE, AFTER to apply to COMMENT keywords W. Landsman
  ;       June 2003, Added SAVECOMMENT keyword    W. Landsman
  ;       Jan 2004, If END is missing, then add it at the end W. Landsman
  ;
  ;-
  if n_params() lt 3 then begin ; Need at least 3 parameters
    print, 'Syntax - Sxaddpar, Header, Name,  Value, [Comment, Postion'
    print, '                      BEFORE = ,AFTER = , FORMAT =, /SAVECOMMENT]'
    return
  endif

  ; Define a blank line and the END line

  ENDLINE = 'END' + string(replicate(32b, 77)) ; END line
  BLANK = string(replicate(32b, 80)) ; BLANK line
  ;
  ; If Location parameter not defined, set it equal to 'END     '
  ;
  if (n_params() gt 4) then loc = strupcase(Location) else $
    if keyword_set(before) then loc = strupcase(before) else $
      if keyword_set(after) then loc = strupcase(after) else $
        if keyword_set(pdu) then loc = 'BEGIN EX' else $
          loc = 'END'

  while strlen(loc) lt 8 do loc = loc + ' '

  if n_params() lt 4 then Comment = '' ; Is comment field specified?

  n = n_elements(Header) ; # of lines in FITS header
  if (n eq 0) then begin ; header defined?
    Header = strarr(10) ; no, make it.
    Header[0] = ENDLINE
    n = 10
  endif else begin
    s = size(Header) ; check for string type
    if (s[0] ne 1) or (s[2] ne 7) then $
      message, 'FITS Header (first parameter) must be a string array'
  endelse

  ; Make sure Name is 8 characters long

  nn = string(replicate(32b, 8)) ; 8 char name
  strput, nn, strupcase(Name) ; insert name

  ; Extract first 8 characters of each line of header, and locate END line

  keywrd = strmid(Header, 0, 8) ; Header keywords
  iend = where(keywrd eq 'END     ', nfound)
  ;
  ; If no END, then add it.  Either put it after the last non-null string, or
  ; append it to the end.
  ;
  if nfound eq 0 then begin
    ii = where(strtrim(Header) ne '', nfound)
    ii = max(ii) + 1
    if (nfound eq 0) or (ii eq n_elements(Header)) then begin
      Header = [Header, ENDLINE]
      n = n + 1
    endif else Header[ii] = ENDLINE
    keywrd = strmid(Header, 0, 8)
    iend = where(keywrd eq 'END     ', nfound)
  endif
  ;
  iend = iend[0] > 0 ; make scalar

  ; History, comment and "blank" records are treated differently from the
  ; others.  They are simply added to the header array whether there are any
  ; already there or not.

  if (nn eq 'HISTORY ') or (nn eq 'COMMENT ') or $
    (nn eq '        ') then begin ; add history record?
    ;
    ; If the header array needs to grow, then expand it in increments of 5 lines.
    ;

    if iend ge (n - 1) then begin
      Header = [Header, replicate(BLANK, 5)] ; yes, add 5.
      n = n_elements(Header)
    endif

    ; Format the record

    newline = BLANK
    strput, newline, nn + string(Value), 0

    ;
    ; If a history record, then append to the record just before the end.
    ;
    if nn eq 'HISTORY ' then begin
      Header[iend] = newline ; add history rec.
      Header[iend + 1] = ENDLINE
      ;
      ; The comment record is placed immediately after the last previous comment
      ; record, or immediately before the first history record, unless overridden by
      ; either the BEFORE or AFTER keywords.
      ;
    endif else if nn eq 'COMMENT ' then begin
      if loc eq 'END     ' then loc = 'COMMENT '
      iloc = where(keywrd eq loc, nloc)
      if nloc eq 0 then iloc = where(keywrd eq 'HISTORY ', nloc)
      if nloc gt 0 then begin
        i = iloc[nloc - 1]
        if keyword_set(after) or (loc eq 'COMMENT ') then i = i + 1 < iend
        if i gt 0 then Header = [Header[0 : i - 1], newline, Header[i : n - 1]] $
        else Header = [newline, Header[0 : n - 1]]
      endif else begin
        Header[iend] = newline
        Header[iend + 1] = ENDLINE
      endelse

      ;
      ; The "blank" record is placed immediately after the last previous "blank"
      ; record, or immediately before the first comment or history record, unless
      ; overridden by either the BEFORE or AFTER keywords.
      ;
    endif else begin
      if loc eq 'END     ' then loc = '       '
      iloc = where(keywrd[0 : iend] eq loc, nloc)
      if nloc gt 0 then begin
        i = iloc[0]
        if keyword_set(after) and loc ne 'HISTORY ' then i = i + 1 < iend
        if i gt 0 then Header = [Header[0 : i - 1], newline, Header[i : n - 1]] $
        else Header = [newline, Header[0 : n - 1]]
      endif else begin
        iloc = where(keywrd eq 'COMMENT ', nloc)
        if nloc eq 0 then iloc = where(keywrd eq 'HISTORY ', nloc)
        if nloc gt 0 then begin
          i = iloc[0]
          if i gt 0 then Header = [Header[0 : i - 1], newline, Header[i : n - 1]] $
          else Header = [newline, Header[0 : n - 1]]
        endif else begin
          Header[iend] = newline
          Header[iend + 1] = ENDLINE
        endelse
      endelse
    endelse
    RETURN
  endif

  ; Find location to insert keyword.   Save the existing comment if user did
  ; not supply a new one.   Comment starts after column 32 for numeric data,
  ; after the slash (but at least after column 20) for string data.

  ncomment = Comment
  ipos = where(keywrd eq nn, nfound)
  if nfound gt 0 then begin
    i = ipos[nfound - 1]
    if Comment eq '' or keyword_set(savecom) then begin ; save comment?
      if strmid(Header[i], 10, 1) ne '''' then $
        Comment = strmid(Header[i], 32, 48) else begin
        slash = strpos(Header[i], '/', 20)
        if slash ne -1 then $
          ncomment = strmid(Header[i], slash + 1, 80) else $
          ncomment = string(replicate(32b, 80))
      endelse
    endif
    goto, replace
  endif

  if loc ne '' then begin
    iloc = where(keywrd eq loc, nloc)
    if nloc gt 0 then begin
      i = iloc[0]
      if keyword_set(after) and loc ne 'HISTORY ' then i = i + 1 < iend
      if i gt 0 then Header = [Header[0 : i - 1], BLANK, Header[i : n - 1]] $
      else Header = [BLANK, Header[0 : n - 1]]
      goto, replace
    endif
  endif

  ; At this point keyword and location parameters were not found, so a new
  ; line is added at the end of the FITS header

  if iend lt (n - 1) then begin ; Not found, add more?
    Header[iend + 1] = ENDLINE ; no, already long enough.
    i = iend ; position to add.
  endif else begin ; must lengthen.
    Header = [Header, replicate(BLANK, 5)] ; add an element on the end
    Header[n] = ENDLINE ; save "END"
    i = n - 1 ; add to end
  end

  ; Now put value into keyword at line i

  replace:
  h = BLANK ; 80 blanks
  strput, h, nn + '= ' ; insert name and =.
  apost = '''' ; quote a quote
  type = size(Value) ; get type of value parameter
  if type[0] ne 0 then $
    message, 'Keyword Value (third parameter) must be scalar'

  case type[1] of ; which type?

    7: begin
      upval = strupcase(Value) ; force upper case.
      if (upval eq 'T') or (upval eq 'F') then begin
        strput, h, upval, 29 ; insert logical value.
        strput, h, ' /', 30 ; add ' /'
        strput, h, ncomment, 32 ; add comment
      end else begin ; other string?
        strput, h, apost + Value + apost, 10
        strput, h, ' /' + ncomment, 30 > (12 + strlen(Value)) ; insert string val
      endelse
      Header[i] = h ; save line
    endcase

    5: begin
      if (n_elements(format) eq 1) then $ ; use format keyword
        v = string(Value, format = '(' + strupcase(format) + ')') $
      else v = string(Value, format = '(G19.12)')
      s = strlen(v) ; right justify
      strput, h, v, (30 - s) > 10
    end

    else: begin
      if (n_elements(format) eq 1) then $ ; use format keyword
        v = string(Value, format = '(' + strupcase(format) + ')') else $
        v = strtrim(strupcase(Value), 2)
      ; convert to string, default format
      s = strlen(v) ; right justify
      strput, h, v, (30 - s) > 10 ; insert
    end
  endcase

  if type[1] ne 7 then begin
    strput, h, ' /', 30 ; add ' /'
    strput, h, ncomment, 32 ; add comment
    Header[i] = h ; save line
  endif

  return
end