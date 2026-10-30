;+
; NAME:
;   CMPRODUCT
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;
; PURPOSE:
;   CMPRODUCT() is the multiplicative equivalent of TOTAL().
;
; CALLING SEQUENCE:
;   Result = CMPRODUCT(ARRAY)
;
; DESCRIPTION:
;
;   Calculates the product of all the elements of an array.  Vector
;   multiplication in groups of powers of two make this operation
;   faster than a simple FOR loop.  The number of actual
;   multiplications is still N_ELEMENTS(ARRAY).  Double precision
;   should be used for the highest accuracy when multiplying many
;   numbers.
;
; INPUTS:
;
;   ARRAY - Array of elements to multiply together.  For instance,
;           ARRAY could contain the dimensions of another array--then
;           CMPRODUCT(ARRAY) would be the total number of elements of
;           that other array.
;
; RETURNS:
;  The result of the function is the total product of all the elements
;  of ARRAY.
;
; EXAMPLE:
;
; SEE ALSO:
;
;   TOTAL, PRODUCT (from Astronomy User's Library)
;
; MODIFICATION HISTORY:
;   Written, CM, 28 Mar 2000
;     (based on outline of PRODUCT by William Thompson)
;
;  $Id: cmproduct.pro,v 1.2 2001/03/25 18:10:42 craigm Exp $
;
;-
; Copyright (C) 2000, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
; -

function CMPRODUCT, ARRAY
  compile_opt idl2
  on_error, 2
  ;
  ; Check the number of parameters.
  ;
  if n_params() ne 1 then message, 'Syntax:  Result = PRODUCT(ARRAY)'
  ;
  ; Check the type of ARRAY.
  ;
  SZ = size(ARRAY)
  TYPE = SZ[SZ[0] + 1]
  if TYPE eq 0 then message, 'ARRAY not defined'
  if TYPE eq 7 then message, 'Operation illegal with string arrays'
  if TYPE eq 8 then message, 'Operation illegal with structures'
  ;
  ; Calculate the product.
  ;
  X = ARRAY
  N = n_elements(X)
  while N gt 1 do begin
    if (N mod 2) eq 1 then X[0] = X[0] * X[N - 1]
    N2 = floor(N / 2)
    X = X[0 : N2 - 1] * X[N2 : *]
    N = N2
  endwhile
  ;
  RETURN, X[0]
end