;+
; NAME:
;   CMSET_OP
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;
; PURPOSE:
;   Performs an AND, OR, or XOR operation between two sets
;
; CALLING SEQUENCE:
;   SET      = CMSET_OP(A, OP, B)
;
; DESCRIPTION:
;
;   SET_OP performs three common operations between two sets.  The
;   three supported functions of OP are:
;
;        OP          Meaning
;      'AND' - to find the intersection of A and B;
;      'OR'  - to find the union of A and B;
;      'XOR' - to find the those elements who are members of A or B
;              but not both;
;
;   Sets as defined here are one dimensional arrays composed of
;   numeric or string types.  Comparisons of equality between elements
;   are done using the IDL EQ operator.
;
;   The complements of either set can be taken as well, by using the
;   NOT1 and NOT2 keywords.  For example, it may be desireable to find
;   the elements in A but not B, or B but not A (they are different!).
;   The following IDL expressions achieve each of those effects:
;
;      SET = CMSET_OP(A, 'AND', /NOT2, B)   ; A but not B
;      SET = CMSET_OP(/NOT1, A, 'AND', B)   ; B but not A
;
;   Note the distinction between NOT1 and NOT2.  NOT1 refers to the
;   first set (A) and NOT2 refers to the second (B).  Their ordered
;   placement in the calling sequence is entirely optional, but the
;   above ordering makes the logical meaning explicit.
;
;   NOT1 and NOT2 can only be set for the 'AND' operator, and never
;   simultaneously.  This is because the results of an operation with
;   'OR' or 'XOR' and any combination of NOTs -- or with 'AND' and
;   both NOTs -- formally cannot produce a defined result.
;
;   The implementation depends on the type of operands.  For integer
;   types, a fast technique using HISTOGRAM is used.  However, this
;   algorithm becomes inefficient when the dynamic range in the data
;   is large.  For those cases, and for other data types, a technique
;   based on SORT() is used.  Thus the compute time should scale
;   roughly as (A+B)*ALOG(A+B) or better, rather than (A*B) for the
;   brute force approach.  For large arrays this is a significant
;   benefit.
;
; INPUTS:
;
;   A, B - the two sets to be operated on.  A one dimensional array of
;          either numeric or string type.  A and B must be of the same
;          type.  Empty sets are permitted, and are either represented
;          as an undefined variable, or by setting EMPTY1 or EMPTY2.
;
;   OP - a string, the operation to be performed.  Must be one of
;        'AND', 'OR' or 'XOR' (lower or mixed case is permitted).
;        Other operations will cause an error message to be produced.
;
; KEYWORDS:
;
;   NOT1, NOT2 - if set and OP is 'AND', then the complement of A (for
;                NOT1) or B (for NOT2) will be used in the operation.
;                NOT1 and NOT2 cannot be set simultaneously.
;
;   EMPTY1, EMPTY2 - if set, then A (for EMPTY1) or B (for EMPTY2) are
;                    assumed to be the empty set.  The actual values
;                    passed as A or B are then ignored.
;
;   INDEX - if set, then return a list of indices instead of the array
;           values themselves.  The "slower" set operations are always
;           performed in this case.
;
;           The indices refer to the *combined* array [A,B].  To
;           clarify, in the following call: I = CMSET_OP(..., /INDEX);
;           returned values from 0 to NA-1 refer to A[I], and values
;           from NA to NA+NB-1 refer to B[I-NA].
;
;   COUNT - upon return, the number of elements in the result set.
;           This is only important when the result set is the empty
;           set, in which case COUNT is set to zero.
;
; RETURNS:
;
;   The resulting set as a one-dimensional array.  The set may be
;   represented by either an array of data values (default), or an
;   array of indices (if INDEX is set).  Duplicate elements, if any,
;   are removed, and element order may not be preserved.
;
;   The empty set is represented as a return value of -1L, and COUNT
;   is set to zero.  Note that the only way to recognize the empty set
;   is to examine COUNT.
;
; SEE ALSO:
;
;   SET_UTILS.PRO by RSI
;
; MODIFICATION HISTORY:
;   Written, CM, 23 Feb 2000
;   Added empty set capability, CM, 25 Feb 2000
;   Documentation clarification, CM 02 Mar 2000
;   Incompatible but more consistent reworking of EMPTY keywords, CM,
;     04 Mar 2000
;   Minor documentation clarifications, CM, 26 Mar 2000
;   Corrected bug in empty_arg special case, CM 06 Apr 2000
;   Add INDEX keyword, CM 31 Jul 2000
;   Clarify INDEX keyword documentation, CM 06 Sep 2000
;   Made INDEX keyword always force SLOW_SET_OP, CM 06 Sep 2000
;   Added CMSET_OP_UNIQ, and ability to select FIRST_UNIQUE or
;     LAST_UNIQUE values, CM, 18 Sep 2000
;   Removed FIRST_UNIQUE and LAST_UNIQUE, and streamlined
;     CMSET_OP_UNIQ until problems with SORT can be understood, CM, 20
;     Sep 2000 (thanks to Ben Tupper)
;   Still trying to get documentation of INDEX and NOT right, CM, 28
;     Sep 2000 (no code changes)
;   Correct bug for AND case, when input sets A and B each only have
;     one unique value, and the values are equal.  CM, 04 Mar 2004
;     (thanks to James B. jbattat at cfa dot harvard dot edu)
;   Add support for the cases where the input data types are mixed,
;      but still compatible; also, attempt to return the same data
;      type that was passed in; CM, 05 Feb 2005
;   Fix bug in type checking (thanks to "marit"), CM, 10 Dec 2005
;   Work around a stupidity in the built-in IDL HISTOGRAM routine,
;      which tries to "help" you by restricting the MIN/MAX to the
;      range of the input variable (thanks to Will Maddox), CM, 16 Jan 2006
;
;  $Id: cmset_op.pro,v 1.6 2006/01/16 19:45:22 craigm Exp $
;
;-
; Copyright (C) 2000, 2004, 2005, 2006, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
; -

; ; Utility function, similar to UNIQ, but allowing choice of taking
; ; first or last unique element, or non-unique elements.
; ; Unfortunately this doesn't work because of implementation dependent
; ; versions of the SORT() function.

; function cmset_op_uniq, a, first=first, non=non, count=ct, sort=sortit
; if n_elements(a) LE 1 then return, 0L
; sh = (2L*keyword_set(first)-1L)*(-2L*keyword_set(non)+1)
;
; if keyword_set(sortit) then begin
; ;; Sort it manually
; ii = sort(a) & b = a(ii)
; if keyword_set(non) then wh = where(b EQ shift(b, sh), ct) $
; else                     wh = where(b NE shift(b, sh), ct)
; if ct GT 0 then return, ii(wh)
; endif else begin
; ;; Use the user's values directly
; if keyword_set(non) then wh = where(a EQ shift(a, sh), ct) $
; else                     wh = where(a NE shift(a, sh), ct)
; if ct GT 0 then return, wh
; endelse
;
; if keyword_set(first) then return, 0L else return, n_elements(a)-1
; end

; ; Simplified version of CMSET_OP_UNIQ which sorts, and takes the
; ; "first" value, whatever that may mean.
function cmset_op_uniq, a
  compile_opt idl2
  if n_elements(a) le 1 then return, 0l

  ii = sort(a)
  b = a[ii]
  wh = where(b ne shift(b, + 1l), ct)
  if ct gt 0 then return, ii[wh]

  return, 0l
end

function cmset_op, a, op0, b, not1 = not1, not2 = not2, count = count, $
  empty1 = empty1, empty2 = empty2, maxarray = ma, index = index
  compile_opt idl2

  on_error, 2 ; ; return on error
  count = 0l
  index0 = -1l
  ; ; Histogram technique is used for array sizes < 32,000 elements
  if n_elements(ma) eq 0 then ma = 32l * 1024l

  ; ; Check the number of arguments
  if n_params() lt 3 then begin
    arg_err:
    message, 'USAGE: SET = CMSET_OP(A, OP, B [, COUNT=ct])', /info
    message, '  KEYWORDS: /NOT1, /NOT2, /EMPTY1, /EMPTY2, INDEX', /info
    return, -1l
  endif
  if n_elements(op0) eq 0 then goto, arg_err
  kind = keyword_set(index)
  fst = 1l
  if keyword_set(last) then fst = 0l
  if keyword_set(first) then fst = 1l

  ; ; Check the operation
  sz = size(op0)
  if sz[sz[0] + 1] ne 7 then begin
    op_err:
    message, 'ERROR: OP must be ''AND'', ''OR'' or ''XOR'''
    return, -1l
  endif
  op = strupcase(op0)
  if op ne 'AND' and op ne 'OR' and op ne 'XOR' then goto, op_err

  ; ; Check NOT1 and NOT2
  if keyword_set(not1) and keyword_set(not2) then begin
    message, 'ERROR: NOT1 and NOT2 cannot be set simultaneously'
    return, -1l
  endif
  if (keyword_set(not1) or keyword_set(not2)) and $
    (op eq 'OR' or op eq 'XOR') then begin
    message, 'ERROR: NOT1 and NOT2 cannot be set with ''OR'' or ''XOR'''
    return, -1l
  endif

  ; ; Special cases for empty set
  n1 = n_elements(a)
  n2 = n_elements(b)
  if keyword_set(empty1) then n1 = 0l
  if keyword_set(empty2) then n2 = 0l
  if n1 eq 0 or n2 eq 0 then begin
    ; ; Eliminate duplicates
    if n1 gt 0 then a1 = cmset_op_uniq(a)
    if n2 gt 0 then b1 = cmset_op_uniq(b)
    n1 = n_elements(a1) < n1
    n2 = n_elements(b1) < n2
    case op of
      'OR': if n1 eq 0 then goto, ret_a1 else goto, ret_b1
      'XOR': if n1 eq 0 then goto, ret_b1 else goto, ret_a1
      'AND': begin
        if keyword_set(not1) and n1 eq 0 then goto, ret_b1
        if keyword_set(not2) and n2 eq 0 then goto, ret_a1
        return, -1l
      end
    endcase
    return, -1l
    ret_a1:
    count = n1
    if kind then begin
      if count gt 0 then return, a1 else return, -1l
    endif
    if count gt 0 then return, a[a1] else return, -1l
    ret_b1:
    count = n2
    if kind then begin
      if count gt 0 then return, b1 + n1 else return, -1l
    endif
    if count gt 0 then return, b[b1] else return, -1l
  endif

  ; ; Allow data to have different types, but they must be at least of
  ; ; the same "base" type.  That is, you can't combine a number with a
  ; ; string, etc.
  ; ; basetype 0:undefined 1:real number 6:complex number 7:string
  ; ;     8:structure 10:pointer 11:object

  ; ;          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
  basetype = [0, 1, 1, 1, 1, 1, 6, 7, 8, 6, 10, 11, 1, 1, 1, 1]

  ; ; Check types of operands
  sz1 = size(a)
  tp1 = sz1[sz1[0] + 1]
  sz2 = size(b)
  tp2 = sz2[sz2[0] + 1]
  if tp1 lt 0 or tp1 ge 16 or tp2 lt 0 or tp2 ge 16 then begin
    message, 'ERROR: unrecognized data types for operands'
    return, -1
  endif
  if basetype[tp1] ne basetype[tp2] then begin
    type1_err:
    message, 'ERROR: both A and B must be of the same type'
    return, -1l
  endif
  if tp1 eq 8 or tp1 eq 10 or tp1 eq 11 then begin
    type2_err:
    message, 'ERROR: operands must be a numeric or string type'
    return, -1l
  endif

  ; ; Now use two different kinds of algorithms: a slower but more
  ; ; general algorithm for generic types, and the histogram technique
  ; ; for integer types.  Even for integer types, if there is too much
  ; ; dynamic range, then the slow method is used.

  if tp1 ge 4 and tp1 le 9 then begin
    ; ; String and real types, or large int arrays
    slow_set_op:
    case op of
      'OR': begin
        uu = [a, b] ; ; OR is simple; just take unique values
        index0 = cmset_op_uniq(uu)
        count = n_elements(index0)
        if kind then return, index0
        return, uu[index0]
      end

      'XOR': begin
        ; ; Make ordered list of set union
        ai = cmset_op_uniq(a)
        na = n_elements(ai)
        bi = cmset_op_uniq(b)
        nb = n_elements(bi)
        ui = [ai, bi + n1]
        uu = [a, b]
        uu = uu[ui] ; ; Raw union...
        us = sort(uu)
        uu = uu[us] ; ; ...and sort
        if kind then ui = ui[temporary(us)] else ui = 0

        ; ; Values in one set only will not have duplicates
        wh1 = where(uu ne shift(uu, -1), count1)
        if count1 eq 0 then return, -1l
        wh = where(wh1[1 : *] - wh1 eq 1, count)
        if wh1[0] eq 0 then begin
          if count gt 0 then wh = [-1l, wh] else wh = [-1l]
          count = n_elements(wh)
        endif
        if count eq 0 then return, -1
        if kind then return, ui[wh1[wh + 1]]
        return, uu[wh1[wh + 1]]
      end

      'AND': begin
        ; ; Make ordered list of set union
        ai = cmset_op_uniq(a)
        na = n_elements(ai)
        bi = cmset_op_uniq(b)
        nb = n_elements(bi)
        ui = [ai, bi + n1]
        uu = [a, b]
        uu = uu[ui] ; ; Raw union...
        us = sort(uu)
        uu = uu[us] ; ; ...and sort
        if kind then ui = ui[us] else ui = 0

        if not keyword_set(not1) and not keyword_set(not2) then begin
          ; ; Special case: if there are one in each set, and
          ; ; they are equal, then the SHIFT() technique below
          ; ; fails.  Do this one by hand.
          if na eq 1 and nb eq 1 and uu[0] eq uu[1] then begin
            count = 1l
            if kind then return, 0l
            return, [uu[0]]
          endif

          ; ; If neither "NOT" is set, then find duplicates
          us = 0l ; ; Save memory
          wh = where(uu eq shift(uu, -1l), count) ; ; Find non-unique
          if count eq 0 then return, -1l
          ; ; This should always select the element from A
          ; ; rather than B (the smaller of the two)
          if kind then return, (ui[wh] < ui[wh + 1])
          return, uu[wh]
        endif

        ; ; For "NOT" cases, we need to identify by set
        ii = make_array(na + nb, value = 1b)
        if keyword_set(not1) then ii[0 : na - 1] = 0
        if keyword_set(not2) then ii[na : *] = 0
        ii = ii[temporary(us)]

        ; ; Remove any duplicates
        wh1 = where(uu eq shift(uu, -1l), count1) ; ; Find non-unique
        if count1 gt 0 then ii[[wh1, wh1 + 1]] = 0
        ; ; Remainder is the desired set
        wh = where(ii, count)
        if count eq 0 then return, -1l
        if kind then return, ui[wh]
        return, uu[wh]
      end
    endcase
    return, -1l ; ; DEFAULT CASE
  endif else begin
    ; ; INDEX keyword forces the "slow" operation
    if kind then goto, slow_set_op

    ; ; Integer types - use histogram technique if the data range
    ; ; is small enough, otherwise use the "slow" technique above
    min1 = min(a, max = max1)
    min2 = min(b, max = max2)
    minn = min1 < min2
    maxx = max1 > max2
    nbins = maxx - minn + 1
    if (maxx - minn) gt floor(ma[0]) then goto, slow_set_op

    ; ; Work around a stupidity in the built-in IDL HISTOGRAM routine
    if (tp1 eq 2 or tp2 eq 2) and (minn lt -32768 or maxx gt 32767) then $
      goto, slow_set_op

    ; ; Following operations create a histogram of the integer values.
    ha = histogram(a, min = minn, max = maxx) < 1
    hb = histogram(b, min = minn, max = maxx) < 1

    ; ; Compute NOT cases
    if keyword_set(not1) then ha = 1b - ha
    if keyword_set(not2) then hb = 1b - hb
    case op of
      ; ; Boolean operations
      'AND': mask = temporary(ha) and temporary(hb)
      'OR': mask = temporary(ha) or temporary(hb)
      'XOR': mask = temporary(ha) xor temporary(hb)
    endcase

    wh = where(temporary(mask), count)
    if count eq 0 then return, -1l

    result = temporary(wh + minn)
    if tp1 ne tp2 then return, result
    szr = size(result)
    tpr = szr[szr[0] + 1]

    ; ; Cast to the original type if necessary
    if tpr ne tp1 then begin
      fresult = make_array(n_elements(result), type = tp1)
      fresult[0] = temporary(result)
      result = temporary(fresult)
    endif

    return, result
  endelse

  return, -1l ; ; DEFAULT CASE
end

; Here is how I did the INDEX stuff with fast histogramming.  It
; works, but is complicated, so I forced it to go to SLOW_SET_OP.
; ha = histogram(a, min=minn, max=maxx, reverse=ra) < 1
; rr = ra(0:nbins) & mask = rr NE rr(1:*) & ra = ra(rr)*mask-1L+mask
; hb = histogram(b, min=minn, max=maxx, reverse=rb) < 1
; rr = rb(0:nbins) & mask = rr NE rr(1:*) & rb = rb(rr)*mask-1L+mask
; ...  AND/OR/XOR NOT masking here ...
; ra = ra(wh) & rb = rb(wh)
; return, ra*(ra GE 0) + (rb+n1)*(ra LT 0) ;; is last 'ra' right?