;+
; NAME:
;   MPNORMLIM
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE:
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Compute confidence limits for normally distributed variable
;
; MAJOR TOPICS:
;   Curve and Surface Fitting, Statistics
;
; CALLING SEQUENCE:
;   Z = MPNORMLIM(PROB, [/CLEVEL, /SLEVEL ])
;
; DESCRIPTION:
;
;  The function MPNORMLIM() computes confidence limits of a normally
;  distributed variable (with zero mean and unit variance), for a
;  desired probability level.  The returned values, Z, are the
;  limiting values: a the magnitude of a normally distributed value
;  is greater than Z by chance with a probability PROB:
;
;    P_NORM(ABS(X) > Z) = PROB
;
;  In specifying the probability level the user has two choices:
;
;    * give the confidence level (default);
;
;    * give the significance level (i.e., 1 - confidence level) and
;      pass the /SLEVEL keyword; OR
;
;  Note that /SLEVEL and /CLEVEL are mutually exclusive.
;
; INPUTS:
;
;   PROB - scalar or vector number, giving the desired probability
;          level as described above.
;
; RETURNS:
;
;  Returns a scalar or vector of normal confidence limits.
;
; KEYWORD PARAMETERS:
;
;   SLEVEL - if set, then PROB describes the significance level.
;
;   CLEVEL - if set, then PROB describes the confidence level
;            (default).
;
; EXAMPLE:
;
;   print, mpnormlim(0.99d, /clevel)
;
;   Print the 99% confidence limit for a normally distributed
;   variable.  In this case it is about 2.58 sigma.
;
; REFERENCES:
;
;   Algorithms taken from CEPHES special function library, by Stephen
;   Moshier. (http://www.netlib.org/cephes/)
;
; MODIFICATION HISTORY:
;   Completed, 1999, CM
;   Documented, 16 Nov 2001, CM
;   Reduced obtrusiveness of common block and math error handling, 18
;     Nov 2001, CM
;   Convert to IDL 5 array syntax (!), 16 Jul 2006, CM
;   Move STRICTARR compile option inside each function/procedure, 9 Oct 2006
;   Add usage message, 24 Nov 2006, CM
;
;  $Id: mpnormlim.pro,v 1.6 2006/11/25 01:44:13 craigm Exp $
;-
; Copyright (C) 1997-2001, 2006, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
; -

forward_function cephes_polevl, cephes_ndtri, mpnormlim

; ; Set machine constants, once for this session.  Double precision
; ; only.
pro cephes_setmachar
  compile_opt strictarr
  common cephes_machar, cephes_machar_vals
  if n_elements(cephes_machar_vals) gt 0 then return

  if (!version.release) lt 5 then dummy = check_math(1, 1)

  mch = machar(/double)
  machep = mch.eps
  maxnum = mch.xmax
  minnum = mch.xmin
  maxlog = alog(mch.xmax)
  minlog = alog(mch.xmin)
  maxgam = 171.624376956302725d

  cephes_machar_vals = {machep: machep, maxnum: maxnum, minnum: minnum, $
    maxlog: maxlog, minlog: minlog, maxgam: maxgam}

  if (!version.release) lt 5 then dummy = check_math(0, 0)
  return
end

function cephes_polevl, x, coef
  compile_opt strictarr
  ans = coef[0]
  nc = n_elements(coef)
  for i = 1l, nc - 1 do ans = ans * x + coef[i]
  return, ans
end

function cephes_ndtri, y0
  ;
  ; Inverse of Normal distribution function
  ;
  ;
  ;
  ; SYNOPSIS:
  ;
  ; double x, y, ndtri();
  ;
  ; x = ndtri( y );
  ;
  ;
  ;
  ; DESCRIPTION:
  ;
  ; Returns the argument, x, for which the area under the
  ; Gaussian probability density function (integrated from
  ; minus infinity to x) is equal to y.
  ;
  ;
  ; For small arguments 0 < y < exp(-2), the program computes
  ; z = sqrt( -2.0 * log(y) );  then the approximation is
  ; x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).
  ; There are two rational functions P/Q, one for 0 < y < exp(-32)
  ; and the other for y up to exp(-2).  For larger arguments,
  ; w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).
  ;
  ;
  ; ACCURACY:
  ;
  ; Relative error:
  ; arithmetic   domain        # trials      peak         rms
  ; DEC      0.125, 1         5500       9.5e-17     2.1e-17
  ; DEC      6e-39, 0.135     3500       5.7e-17     1.3e-17
  ; IEEE     0.125, 1        20000       7.2e-16     1.3e-16
  ; IEEE     3e-308, 0.135   50000       4.6e-16     9.8e-17
  ;
  ;
  ; ERROR MESSAGES:
  ;
  ; message         condition    value returned
  ; ndtri domain       x <= 0        -MAXNUM
  ; ndtri domain       x >= 1         MAXNUM
  compile_opt strictarr
  common cephes_ndtri_data, s2pi, p0, q0, p1, q1, p2, q2

  if n_elements(s2pi) eq 0 then begin
    s2pi = sqrt(2.d * !dpi)
    p0 = [-5.99633501014107895267d1, 9.80010754185999661536d1, $
      - 5.66762857469070293439d1, 1.39312609387279679503d1, $
      - 1.23916583867381258016d0]
    q0 = [1.d, $
      1.95448858338141759834d0, 4.67627912898881538453d0, $
      8.63602421390890590575d1, -2.25462687854119370527d2, $
      2.00260212380060660359d2, -8.20372256168333339912d1, $
      1.59056225126211695515d1, -1.18331621121330003142d0]
    p1 = [4.05544892305962419923d0, 3.15251094599893866154d1, $
      5.71628192246421288162d1, 4.40805073893200834700d1, $
      1.46849561928858024014d1, 2.18663306850790267539d0, $
      - 1.40256079171354495875d-1, -3.50424626827848203418d-2, $
      - 8.57456785154685413611d-4]
    q1 = [1.d, $
      1.57799883256466749731d1, 4.53907635128879210584d1, $
      4.13172038254672030440d1, 1.50425385692907503408d1, $
      2.50464946208309415979d0, -1.42182922854787788574d-1, $
      - 3.80806407691578277194d-2, -9.33259480895457427372d-4]
    p2 = [3.23774891776946035970d0, 6.91522889068984211695d0, $
      3.93881025292474443415d0, 1.33303460815807542389d0, $
      2.01485389549179081538d-1, 1.23716634817820021358d-2, $
      3.01581553508235416007d-4, 2.65806974686737550832d-6, $
      6.23974539184983293730d-9]
    q2 = [1.d, $
      6.02427039364742014255d0, 3.67983563856160859403d0, $
      1.37702099489081330271d0, 2.16236993594496635890d-1, $
      1.34204006088543189037d-2, 3.28014464682127739104d-4, $
      2.89247864745380683936d-6, 6.79019408009981274425d-9]
  endif

  common cephes_machar, machvals
  MAXNUM = machvals.maxnum

  if y0 le 0 then begin
    message, 'ERROR: domain', /info
    return, -MAXNUM
  endif
  if y0 ge 1 then begin
    message, 'ERROR: domain', /info
    return, MAXNUM
  endif

  code = 1
  y = y0
  exp2 = exp(-2.d)
  if y gt (1.d - exp2) then begin
    y = 1.d - y
    code = 0
  endif
  if y gt exp2 then begin
    y = y - 0.5
    y2 = y * y
    x = y + y * y2 * cephes_polevl(y2, p0) / cephes_polevl(y2, q0)
    x = x * s2pi
    return, x
  endif

  x = sqrt(-2.d * alog(y))
  x0 = x - alog(x) / x
  z = 1.d / x
  if x lt 8. then $
    x1 = z * cephes_polevl(z, p1) / cephes_polevl(z, q1) $
  else $
    x1 = z * cephes_polevl(z, p2) / cephes_polevl(z, q2)

  x = x0 - x1
  if code ne 0 then x = -x
  return, x
end

; MPNORMLIM - given a probability level, return the corresponding
; "sigma" level.
;
; p - Either the significance level (if SLEVEL is set) or the
; confidence level (if CLEVEL is set).  This should be the
; two-tailed level, ie:
;
; * SLEVEL:   p = Prob(|z| > z0)
; * CLEVEL:   p = Prob(|z| < z0)
;

function mpnormlim, p, clevel = clevel, slevel = slevel
  compile_opt strictarr

  if n_params() eq 0 then begin
    message, 'USAGE: Z = MPNORMLIM(PROB, [/CLEVEL, /SLEVEL ])', /info
    return, !values.d_nan
  endif

  cephes_setmachar ; ; Set machine constants

  ; ; Default is to assume the confidence level
  if n_elements(clevel) eq 0 then clevel = 1
  y = 0 * p

  ; ; cephes_ndtri accepts the integrated probability from negative
  ; ; infinity to z, so we have to compute.
  if keyword_set(slevel) then begin
    p1 = 0.5d * p ; ; Take only one of the two tails
    for i = 0l, n_elements(y) - 1 do begin
      y[i] = -cephes_ndtri(p1[i])
    endfor
  endif else if keyword_set(clevel) then begin
    p1 = 0.5d + 0.5d * p ; ; On binary computers this computation is
    ; ; exact (to the machine precision), so don't worry about it.
    ; ; This computation shaves off the top half of the confidence
    ; ; region, and then adds the "negative infinity to zero part.
    for i = 0l, n_elements(y) - 1 do begin
      y[i] = cephes_ndtri(p1[i])
    endfor
  endif else begin
    message, 'ERROR: must specify one of CLEVEL or SLEVEL'
  endelse

  return, y
end