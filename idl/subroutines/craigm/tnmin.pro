;+
; NAME:
;   TNMIN
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE:
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Performs function minimization (Truncated-Newton Method)
;
; MAJOR TOPICS:
;   Optimization and Minimization
;
; CALLING SEQUENCE:
;   parms = TNMIN(MYFUNCT, X, FUNCTARGS=fcnargs, NFEV=nfev,
;                 MAXITER=maxiter, ERRMSG=errmsg, NPRINT=nprint,
;                 QUIET=quiet, XTOL=xtol, STATUS=status,
;                 FGUESS=fguess, PARINFO=parinfo, BESTMIN=bestmin,
;                 ITERPROC=iterproc, ITERARGS=iterargs, niter=niter)
;
; DESCRIPTION:
;
;  TNMIN uses the Truncated-Newton method to minimize an arbitrary IDL
;  function with respect to a given set of free parameters.  The
;  user-supplied function must compute the gradient with respect to
;  each parameter.  TNMIN is based on TN.F (TNBC) by Stephen Nash.
;
;  If you want to solve a least-squares problem, to perform *curve*
;  *fitting*, then you will probably want to use the routines MPFIT,
;  MPFITFUN and MPFITEXPR.  Those routines are specifically optimized
;  for the least-squares problem.  TNMIN is suitable for constrained
;  and unconstrained optimization problems with a medium number of
;  variables.  Function *maximization* can be performed using the
;  MAXIMIZE keyword.
;
;  TNMIN is similar to MPFIT in that it allows parameters to be fixed,
;  simple bounding limits to be placed on parameter values, and
;  parameters to be tied to other parameters.  One major difference
;  between MPFIT and TNMIN is that TNMIN does not compute derivatives
;  automatically by default.  See PARINFO and AUTODERIVATIVE below for
;  more discussion and examples.
;
; USER FUNCTION
;
;  The user must define an IDL function which returns the desired
;  value as a single scalar.  The IDL function must accept a list of
;  numerical parameters, P.  Additionally, keyword parameters may be
;  used to pass more data or information to the user function, via the
;  FUNCTARGS keyword.
;
;  The user function should be declared in the following way:
;
;     FUNCTION MYFUNCT, p, dp [, keywords permitted ]
;       ; Parameter values are passed in "p"
;       f  = ....   ; Compute function value
;       dp = ....   ; Compute partial derivatives (optional)
;       return, f
;     END
;
;  The function *must* accept at least one argument, the parameter
;  list P.  The vector P is implicitly assumed to be a one-dimensional
;  array.  Users may pass additional information to the function by
;  composing and _EXTRA structure and passing it in the FUNCTARGS
;  keyword.
;
;  User functions may also indicate a fatal error condition using the
;  ERROR_CODE common block variable, as described below under the
;  TNMIN_ERROR common block definition (by setting ERROR_CODE to a
;  number between -15 and -1).
;
;  EXPLICIT vs. NUMERICAL DERIVATIVES
;
;  By default, the user must compute the function gradient components
;  explicitly using AUTODERIVATIVE=0.  As explained below, numerical
;  derivatives can also be calculated using AUTODERIVATIVE=1.
;
;  For explicit derivatives, the user should call TNMIN using the
;  default keyword value AUTODERIVATIVE=0.  [ This is different
;  behavior from MPFIT, where AUTODERIVATIVE=1 is the default. ] The
;  IDL user routine should compute the gradient of the function as a
;  one-dimensional array of values, one for each of the parameters.
;  They are passed back to TNMIN via "dp" as shown above.
;
;  The derivatives with respect to fixed parameters are ignored; zero
;  is an appropriate value to insert for those derivatives.  Upon
;  input to the user function, DP is set to a vector with the same
;  length as P, with a value of 1 for a parameter which is free, and a
;  value of zero for a parameter which is fixed (and hence no
;  derivative needs to be calculated).  This input vector may be
;  overwritten as needed.
;
;  For numerical derivatives, a finite differencing approximation is
;  used to estimate the gradient values.  Users can activate this
;  feature by passing the keyword AUTODERIVATIVE=1.  Fine control over
;  this behavior can be achieved using the STEP, RELSTEP and TNSIDE
;  fields of the PARINFO structure.
;
;  When finite differencing is used for computing derivatives (ie,
;  when AUTODERIVATIVE=1), the parameter DP is not passed.  Therefore
;  functions can use N_PARAMS() to indicate whether they must compute
;  the derivatives or not.  However there is no penalty (other than
;  computation time) for computing the gradient values and users may
;  switch between AUTODERIVATIVE=0 or =1 in order to test both
;  scenarios.
;
; CONSTRAINING PARAMETER VALUES WITH THE PARINFO KEYWORD
;
;  The behavior of TNMIN can be modified with respect to each
;  parameter to be optimized.  A parameter value can be fixed; simple
;  boundary constraints can be imposed; limitations on the parameter
;  changes can be imposed; properties of the automatic derivative can
;  be modified; and parameters can be tied to one another.
;
;  These properties are governed by the PARINFO structure, which is
;  passed as a keyword parameter to TNMIN.
;
;  PARINFO should be an array of structures, one for each parameter.
;  Each parameter is associated with one element of the array, in
;  numerical order.  The structure can have the following entries
;  (none are required):
;
;     .VALUE - the starting parameter value (but see the START_PARAMS
;              parameter for more information).
;
;     .FIXED - a boolean value, whether the parameter is to be held
;              fixed or not.  Fixed parameters are not varied by
;              TNMIN, but are passed on to MYFUNCT for evaluation.
;
;     .LIMITED - a two-element boolean array.  If the first/second
;                element is set, then the parameter is bounded on the
;                lower/upper side.  A parameter can be bounded on both
;                sides.  Both LIMITED and LIMITS must be given
;                together.
;
;     .LIMITS - a two-element float or double array.  Gives the
;               parameter limits on the lower and upper sides,
;               respectively.  Zero, one or two of these values can be
;               set, depending on the values of LIMITED.  Both LIMITED
;               and LIMITS must be given together.
;
;     .PARNAME - a string, giving the name of the parameter.  The
;                fitting code of TNMIN does not use this tag in any
;                way.
;
;     .STEP - the step size to be used in calculating the numerical
;             derivatives.  If set to zero, then the step size is
;             computed automatically.  Ignored when AUTODERIVATIVE=0.
;
;     .TNSIDE - the sidedness of the finite difference when computing
;               numerical derivatives.  This field can take four
;               values:
;
;                  0 - one-sided derivative computed automatically
;                  1 - one-sided derivative (f(x+h) - f(x)  )/h
;                 -1 - one-sided derivative (f(x)   - f(x-h))/h
;                  2 - two-sided derivative (f(x+h) - f(x-h))/(2*h)
;
;              Where H is the STEP parameter described above.  The
;              "automatic" one-sided derivative method will chose a
;              direction for the finite difference which does not
;              violate any constraints.  The other methods do not
;              perform this check.  The two-sided method is in
;              principle more precise, but requires twice as many
;              function evaluations.  Default: 0.
;
;     .TIED - a string expression which "ties" the parameter to other
;             free or fixed parameters.  Any expression involving
;             constants and the parameter array P are permitted.
;             Example: if parameter 2 is always to be twice parameter
;             1 then use the following: parinfo(2).tied = '2 * P(1)'.
;             Since they are totally constrained, tied parameters are
;             considered to be fixed; no errors are computed for them.
;             [ NOTE: the PARNAME can't be used in expressions. ]
;
;  Future modifications to the PARINFO structure, if any, will involve
;  adding structure tags beginning with the two letters "MP" or "TN".
;  Therefore programmers are urged to avoid using tags starting with
;  these two combinations of letters; otherwise they are free to
;  include their own fields within the PARINFO structure, and they
;  will be ignored.
;
;  PARINFO Example:
;  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
;                       limits:[0.D,0]}, 5)
;  parinfo(0).fixed = 1
;  parinfo(4).limited(0) = 1
;  parinfo(4).limits(0)  = 50.D
;  parinfo(*).value = [5.7D, 2.2, 500., 1.5, 2000.]
;
;  A total of 5 parameters, with starting values of 5.7,
;  2.2, 500, 1.5, and 2000 are given.  The first parameter
;  is fixed at a value of 5.7, and the last parameter is
;  constrained to be above 50.
;
;
; INPUTS:
;
;   MYFUNCT - a string variable containing the name of the function to
;             be minimized (see USER FUNCTION above).  The IDL routine
;             should return the value of the function and optionally
;             its gradients.
;
;   X - An array of starting values for each of the parameters of the
;       model.
;
;       This parameter is optional if the PARINFO keyword is used.
;       See above.  The PARINFO keyword provides a mechanism to fix or
;       constrain individual parameters.  If both X and PARINFO are
;       passed, then the starting *value* is taken from X, but the
;       *constraints* are taken from PARINFO.
;
;
; RETURNS:
;
;   Returns the array of best-fit parameters.
;
;
; KEYWORD PARAMETERS:
;
;   AUTODERIVATIVE - If this is set, derivatives of the function will
;                    be computed automatically via a finite
;                    differencing procedure.  If not set, then MYFUNCT
;                    must provide the (explicit) derivatives.
;                    Default: 0 (explicit derivatives required)
;
;   BESTMIN - upon return, the best minimum function value that TNMIN
;             could find.
;
;   EPSABS - a nonnegative real variable.  Termination occurs when the
;            absolute error between consecutive iterates is at most
;            EPSABS.  Note that using EPSREL is normally preferable
;            over EPSABS, except in cases where ABS(F) is much larger
;            than 1 at the optimal point.  A value of zero indicates
;            the absolute error test is not applied.  If EPSABS is
;            specified, then both EPSREL and EPSABS tests are applied;
;            if either succeeds then termination occurs.
;            Default: 0 (i.e., only EPSREL is applied).
;
;   EPSREL - a nonnegative input variable. Termination occurs when the
;            relative error between two consecutive iterates is at
;            most EPSREL.  Therefore, EPSREL measures the relative
;            error desired in the function.  An additional, more
;            lenient, stopping condition can be applied by specifying
;            the EPSABS keyword.
;            Default: 100 * Machine Precision
;
;   ERRMSG - a string error or warning message is returned.
;
;   FGUESS - provides the routine a guess to the minimum value.
;            Default: 0
;
;   FUNCTARGS - A structure which contains the parameters to be passed
;               to the user-supplied function specified by MYFUNCT via
;               the _EXTRA mechanism.  This is the way you can pass
;               additional data to your user-supplied function without
;               using common blocks.
;
;               Consider the following example:
;                if FUNCTARGS = { XVAL:[1.D,2,3], YVAL:[1.D,4,9]}
;                then the user supplied function should be declared
;                like this:
;                FUNCTION MYFUNCT, P, XVAL=x, YVAL=y
;
;               By default, no extra parameters are passed to the
;               user-supplied function.
;
;   ITERARGS - The keyword arguments to be passed to ITERPROC via the
;              _EXTRA mechanism.  This should be a structure, and is
;              similar in operation to FUNCTARGS.
;              Default: no arguments are passed.
;
;   ITERDERIV - Intended to print function gradient information.  If
;               set, then the ITERDERIV keyword is set in each call to
;               ITERPROC.  In the default ITERPROC, parameter values
;               and gradient information are both printed when this
;               keyword is set.
;
;   ITERPROC - The name of a procedure to be called upon each NPRINT
;              iteration of the TNMIN routine.  It should be declared
;              in the following way:
;
;              PRO ITERPROC, MYFUNCT, p, iter, fnorm, FUNCTARGS=fcnargs, $
;                PARINFO=parinfo, QUIET=quiet, _EXTRA=extra
;                ; perform custom iteration update
;              END
;
;              ITERPROC must accept the _EXTRA keyword, in case of
;              future changes to the calling procedure.
;
;              MYFUNCT is the user-supplied function to be minimized,
;              P is the current set of model parameters, ITER is the
;              iteration number, and FUNCTARGS are the arguments to be
;              passed to MYFUNCT.  FNORM is should be the function
;              value.  QUIET is set when no textual output should be
;              printed.  See below for documentation of PARINFO.
;
;              In implementation, ITERPROC can perform updates to the
;              terminal or graphical user interface, to provide
;              feedback while the fit proceeds.  If the fit is to be
;              stopped for any reason, then ITERPROC should set the
;              common block variable ERROR_CODE to negative value
;              between -15 and -1 (see TNMIN_ERROR common block
;              below).  In principle, ITERPROC should probably not
;              modify the parameter values, because it may interfere
;              with the algorithm's stability.  In practice it is
;              allowed.
;
;              Default: an internal routine is used to print the
;                       parameter values.
;
;   MAXITER - The maximum number of iterations to perform.  If the
;             number is exceeded, then the STATUS value is set to 5
;             and TNMIN returns.
;             Default: 200 iterations
;
;   MAXIMIZE - If set, the function is maximized instead of
;              minimized.
;
;   MAXNFEV - The maximum number of function evaluations to perform.
;             If the number is exceeded, then the STATUS value is set
;             to -17 and TNMIN returns.  A value of zero indicates no
;             maximum.
;             Default: 0 (no maximum)
;
;   NFEV - upon return, the number of MYFUNCT function evaluations
;          performed.
;
;   NITER - upon return, number of iterations completed.
;
;   NPRINT - The frequency with which ITERPROC is called.  A value of
;            1 indicates that ITERPROC is called with every iteration,
;            while 2 indicates every other iteration, etc.
;            Default value: 1
;
;   PARINFO - Provides a mechanism for more sophisticated constraints
;             to be placed on parameter values.  When PARINFO is not
;             passed, then it is assumed that all parameters are free
;             and unconstrained.  Values in PARINFO are never modified
;             during a call to TNMIN.
;
;             See description above for the structure of PARINFO.
;
;             Default value:  all parameters are free and unconstrained.
;
;   QUIET - set this keyword when no textual output should be printed
;           by TNMIN
;
;   STATUS - an integer status code is returned.  All values greater
;            than zero can represent success (however STATUS EQ 5 may
;            indicate failure to converge).  Gaps in the numbering
;            system below are to maintain compatibility with MPFIT.
;            Upon return, STATUS can have one of the following values:
;
;        -18  a fatal internal error occurred during optimization.
;
;        -17  the maximum number of function evaluations has been
;             reached without convergence.
;
;        -16  a parameter or function value has become infinite or an
;             undefined number.  This is usually a consequence of
;             numerical overflow in the user's function, which must be
;             avoided.
;
;        -15 to -1
;             these are error codes that either MYFUNCT or ITERPROC
;             may return to terminate the fitting process (see
;             description of MPFIT_ERROR common below).  If either
;             MYFUNCT or ITERPROC set ERROR_CODE to a negative number,
;             then that number is returned in STATUS.  Values from -15
;             to -1 are reserved for the user functions and will not
;             clash with MPFIT.
;
;	   0  improper input parameters.
;
;	   1  convergence was reached.
;
;          2-4 (RESERVED)
;
;	   5  the maximum number of iterations has been reached
;
;          6-8 (RESERVED)
;
;
; EXAMPLE:
;
;   FUNCTION F, X, DF, _EXTRA=extra  ;; *** MUST ACCEPT KEYWORDS
;     F  = (X(0)-1)^2 + (X(1)+7)^2
;     DF = [ 2D * (X(0)-1), 2D * (X(1)+7) ] ; Gradient
;     RETURN, F
;   END
;
;   P = TNMIN('F', [0D, 0D], BESTMIN=F0)
;   Minimizes the function F(x0,x1) = (x0-1)^2 + (x1+7)^2.
;
;
; COMMON BLOCKS:
;
;   COMMON TNMIN_ERROR, ERROR_CODE
;
;     User routines may stop the fitting process at any time by
;     setting an error condition.  This condition may be set in either
;     the user's model computation routine (MYFUNCT), or in the
;     iteration procedure (ITERPROC).
;
;     To stop the fitting, the above common block must be declared,
;     and ERROR_CODE must be set to a negative number.  After the user
;     procedure or function returns, TNMIN checks the value of this
;     common block variable and exits immediately if the error
;     condition has been set.  By default the value of ERROR_CODE is
;     zero, indicating a successful function/procedure call.
;
;
; REFERENCES:
;
;   TRUNCATED-NEWTON METHOD, TN.F
;      Stephen G. Nash, Operations Research and Applied Statistics
;      Department
;      http://www.netlib.org/opt/tn
;
;   Nash, S. G. 1984, "Newton-Type Minimization via the Lanczos
;      Method," SIAM J. Numerical Analysis, 21, p. 770-778
;
;
; MODIFICATION HISTORY:
;   Derived from TN.F by Stephen Nash with many changes and additions,
;      to conform to MPFIT paradigm, 19 Jan 1999, CM
;   Changed web address to COW, 25 Sep 1999, CM
;   Alphabetized documented keyword parameters, 02 Oct 1999, CM
;   Changed to ERROR_CODE for error condition, 28 Jan 2000, CM
;   Continued and fairly major improvements (CM, 08 Jan 2001):
;      - calling of user procedure is now concentrated in TNMIN_CALL,
;        and arguments are reduced by storing a large number of them
;        in common blocks;
;      - finite differencing is done within TNMIN_CALL; added
;        AUTODERIVATIVE=1 so that finite differencing can be enabled,
;        both one and two sided;
;      - a new procedure to parse PARINFO fields, borrowed from MPFIT;
;        brought PARINFO keywords up to date with MPFIT;
;      - go through and check for float vs. double discrepancies;
;      - add explicit MAXIMIZE keyword, and support in TNMIN_CALL and
;        TNMIN_DEFITER to print the correct values in that case;
;        TNMIN_DEFITER now prints function gradient values if
;        requested;
;      - convert to common-based system of MPFIT for storing machine
;        constants; revert TNMIN_ENORM to simple sum of squares, at
;        least for now;
;      - remove limit on number of function evaluations, at least for
;        now, and until I can understand what happens when we do
;        numerical derivatives.
;   Further changes: more float vs double; disable TNMINSTEP for now;
;     experimented with convergence test in case of function
;     maximization, 11 Jan 2001, CM
;   TNMINSTEP is parsed but not enabled, 13 Mar 2001
;   Major code cleanups; internal docs; reduced commons, CM, 08 Apr
;     2001
;   Continued code cleanups; documentation; the STATUS keyword
;     actually means something, CM, 10 Apr 2001
;   Added reference to Nash paper, CM, 08 Feb 2002
;   Fixed 16-bit loop indices, D. Schelgel, 22 Apr 2003
;   Changed parens to square brackets because of conflicts with
;     limits function.  K. Tolbert, 23 Feb 2005
;   Some documentation clarifications, CM, 09 Nov 2007
;   Ensure that MY_FUNCT returns a scalar; make it more likely that
;     error messages get back out to the user; fatal CATCH'd error
;     now returns STATUS = -18, CM, 17 Sep 2008
;
; TODO
;  - scale derivatives semi-automatically;
;  - ability to scale and offset parameters;
;
;  $Id: tnmin.pro,v 1.17 2008/10/02 03:39:59 craigm Exp $
;-
; Copyright (C) 1998-2001,2002,2003,2007,2008 Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy and distribute unmodified copies for
; non-commercial purposes, and to modify and use for personal or
; internal use, is granted.  All other rights are reserved.
; -

; %% TRUNCATED-NEWTON METHOD:  SUBROUTINES
; FOR OTHER MACHINES, MODIFY ROUTINE MCHPR1 (MACHINE EPSILON)
; WRITTEN BY:  STEPHEN G. NASH
; OPERATIONS RESEARCH AND APPLIED STATISTICS DEPT.
; GEORGE MASON UNIVERSITY
; FAIRFAX, VA 22030
; ******************************************************************

; ; Routine which declares functions and common blocks
pro tnmin_dummy
  compile_opt idl2
  forward_function tnmin_enorm, tnmin_step1, tnmin
  forward_function tnmin_call, tnmin_autoder
  common tnmin_error, error_code
  common tnmin_machar, tnmin_machar_vals
  common tnmin_config, tnmin_tnconfig
  common tnmin_fcnargs, tnmin_tnfcnargs
  common tnmin_work, lsk, lyk, ldiagb, lsr, lyr
  a = 1
  return
end

; ; Following are machine constants that can be loaded once.  I have
; ; found that bizarre underflow messages can be produced in each call
; ; to MACHAR(), so this structure minimizes the number of calls to
; ; one.
pro tnmin_setmachar, double = isdouble
  compile_opt idl2

  common tnmin_machar, tnmin_machar_vals

  ; ; In earlier versions of IDL, MACHAR itself could produce a load of
  ; ; error messages.  We try to mask some of that out here.
  if (!version.release) lt 5 then dummy = check_math(1, 1)

  mch = 0.
  mch = machar(double = keyword_set(isdouble))
  dmachep = mch.eps
  dmaxnum = mch.xmax
  dminnum = mch.xmin
  dmaxlog = alog(mch.xmax)
  dminlog = alog(mch.xmin)
  if keyword_set(isdouble) then $
    dmaxgam = 171.624376956302725d $
  else $
    dmaxgam = 171.624376956302725
  drdwarf = sqrt(dminnum * 1.5) * 10
  drgiant = sqrt(dmaxnum) * 0.1

  tnmin_machar_vals = {machep: dmachep, maxnum: dmaxnum, minnum: dminnum, $
    maxlog: dmaxlog, minlog: dminlog, maxgam: dmaxgam, $
    rdwarf: drdwarf, rgiant: drgiant}

  if (!version.release) lt 5 then dummy = check_math(0, 0)

  return
end

; ; Procedure to parse the parameter values in PARINFO
pro tnmin_parinfo, parinfo, tnames, tag, values, default = def, status = status, $
  n_param = n
  compile_opt idl2

  status = 0
  if n_elements(n) eq 0 then n = n_elements(parinfo)

  if n eq 0 then begin
    if n_elements(def) eq 0 then return
    values = def
    status = 1
    return
  endif

  if n_elements(parinfo) eq 0 then goto, do_default
  if n_elements(tnames) eq 0 then tnames = tag_names(parinfo)
  wh = where(tnames eq tag, ct)

  if ct eq 0 then begin
    do_default:
    if n_elements(def) eq 0 then return
    values = make_array(n, value = def[0])
    values[0] = def
  endif else begin
    values = parinfo.(wh[0])
  endelse

  status = 1
  return
end

; ; Procedure to tie one parameter to another.
pro tnmin_tie, p, _ptied
  compile_opt idl2
  _wh = where(_ptied ne '', _ct)
  if _ct eq 0 then return
  for _i = 0l, _ct - 1 do begin
    _cmd = 'p(' + strtrim(_wh[_i], 2) + ') = ' + _ptied[_wh[_i]]
    _err = execute(_cmd)
    if _err eq 0 then begin
      message, 'ERROR: Tied expression "' + _cmd + '" failed.'
      return
    endif
  endfor
end

function tnmin_autoder, fcn, x, dx, dside = dside
  compile_opt idl2

  common tnmin_machar, machvals
  common tnmin_config, tnconfig

  MACHEP0 = machvals.machep
  DWARF = machvals.minnum
  if n_elements(dside) ne n_elements(x) then dside = tnconfig.dside

  eps = sqrt(MACHEP0)
  h = eps * (1. + abs(x))

  ; ; if STEP is given, use that
  wh = where(tnconfig.step gt 0, ct)
  if ct gt 0 then h[wh] = tnconfig.step(wh)

  ; ; if relative step is given, use that
  wh = where(tnconfig.dstep gt 0, ct)
  if ct gt 0 then h[wh] = abs(tnconfig.dstep(wh) * x[wh])

  ; ; In case any of the step values are zero
  wh = where(h eq 0, ct)
  if ct gt 0 then h[wh] = eps

  ; ; Reverse the sign of the step if we are up against the parameter
  ; ; limit, or if the user requested it.
  mask = (dside eq -1 or (tnconfig.ulimited and (x gt tnconfig.ulimit - h)))
  wh = where(mask, ct)
  if ct gt 0 then h[wh] = -h[wh]

  dx = x * 0.
  f = tnmin_call(fcn, x)
  for i = 0l, n_elements(x) - 1 do begin
    if tnconfig.pfixed(i) eq 1 then goto, next_par
    hh = h[i]

    restart_par:
    xp = x
    xp[i] = xp[i] + hh

    fp = tnmin_call(fcn, xp)

    if abs(dside[i]) le 1 then begin
      ; ; COMPUTE THE ONE-SIDED DERIVATIVE
      dx[i] = (fp - f) / hh
    endif else begin
      ; ; COMPUTE THE TWO-SIDED DERIVATIVE
      xp[i] = x[i] - hh

      fm = tnmin_call(fcn, xp)
      dx[i] = (fp - fm) / (2 * hh)
    endelse
    next_par:
  endfor

  return, f
end

; ; Call user function or procedure, with _EXTRA or not, with
; ; derivatives or not.
function tnmin_call, fcn, x1, dx, fullparam_ = xall
  compile_opt idl2

  ; on_error, 2
  common tnmin_config, tnconfig
  common tnmin_fcnargs, fcnargs

  if keyword_set(tnconfig.qanytied) then tnmin_tie, x, tnconfig.ptied
  ifree = tnconfig.ifree
  ; ; Following promotes the byte array to a floating point array so
  ; ; that users who simply re-fill the array aren't surprised when
  ; ; their gradient comes out as bytes. :-)
  dx = tnconfig.pfixed + x1[0] * 0.

  if n_elements(xall) gt 0 then begin
    x = xall
    x[ifree] = x1
  endif else begin
    x = x1
  endelse

  ; ; Decide whether we are calling a procedure or function
  if tnconfig.proc then proc = 1 else proc = 0
  tnconfig.nfev = tnconfig.nfev + 1

  if n_params() eq 3 then begin
    if tnconfig.autoderiv then $
      f = tnmin_autoder(fcn, x, dx) $
    else if n_elements(fcnargs) gt 0 then $
      f = call_function(fcn, x, dx, _extra = fcnargs) $
    else $
      f = call_function(fcn, x, dx)

    dx = dx[ifree]
    if tnconfig.max then begin
      dx = -dx
      f = -f
    endif
  endif else begin
    if n_elements(fcnargs) gt 0 then $
      f = call_function(fcn, x, _extra = fcnargs) $
    else $
      f = call_function(fcn, x)

    if n_elements(f) ne 1 then begin
      message, 'ERROR: function "' + fcn + '" returned a vector when ' + $
        'a scalar was expected.'
    endif
  endelse

  if n_elements(f) gt 1 then return, f $
  else return, f[0]
end

function tnmin_enorm, vec
  compile_opt idl2

  common tnmin_config, tnconfig
  ; ; Very simple-minded sum-of-squares
  if n_elements(tnconfig) gt 0 then if tnconfig.fastnorm then begin
    ans = sqrt(total(vec ^ 2, 1))
    goto, terminate
  endif

  common tnmin_machar, machvals

  agiant = machvals.rgiant / n_elements(vec)
  adwarf = machvals.rdwarf * n_elements(vec)

  ; ; This is hopefully a compromise between speed and robustness.
  ; ; Need to do this because of the possibility of over- or underflow.
  mx = max(vec, min = mn)
  mx = max(abs([mx, mn]))
  if mx eq 0 then return, vec[0] * 0.

  if mx gt agiant or mx lt adwarf then ans = mx * sqrt(total((vec / mx) ^ 2)) $
  else ans = sqrt(total(vec ^ 2))

  terminate:
  return, ans
end

;
; ROUTINES TO INITIALIZE PRECONDITIONER
;
pro tnmin_initpc, diagb, emat, n, upd1, yksk, gsk, yrsr, lreset
  compile_opt idl2
  ; ; Rename common variables as they appear in INITP3.  Those
  ; ; indicated in all caps are not used or renamed here.
  ; common tnmin_work, lsk, lyk, ldiagb, lsr, lyr
  common tnmin_work, sk, yk, LDIAGB, sr, yr
  ; I    I            I    I

  ; ; From INITP3
  if keyword_set(upd1) then begin
    emat = diagb
  endif else if keyword_set(lreset) then begin
    BSK = diagb * sk
    SDS = total(sk * BSK)
    emat = diagb - diagb * diagb * sk * sk / SDS + yk * yk / yksk
  endif else begin
    BSK = diagb * sr
    SDS = total(sr * BSK)
    SRDS = total(sk * BSK)
    YRSK = total(yr * sk)
    BSK = diagb * sk - BSK * SRDS / SDS + yr * YRSK / yrsr
    emat = diagb - diagb * diagb * sr * sr / SDS + yr * yr / yrsr
    SDS = total(sk * BSK)
    emat = emat - BSK * BSK / SDS + yk * yk / yksk
  endelse

  return
end

pro tnmin_ssbfgs, n, gamma, sj, yj, hjv, hjyj, yjsj, yjhyj, $
  vsj, vhyj, hjp1v
  compile_opt idl2
  ;
  ; SELF-SCALED BFGS
  ;
  DELTA = (1. + gamma * yjhyj / yjsj) * vsj / yjsj - gamma * vhyj / yjsj
  BETA = -gamma * vsj / yjsj
  hjp1v = gamma * hjv + DELTA * sj + BETA * hjyj
  RETURN
end

;
; THIS ROUTINE ACTS AS A PRECONDITIONING STEP FOR THE
; LINEAR CONJUGATE-GRADIENT ROUTINE.  IT IS ALSO THE
; METHOD OF COMPUTING THE SEARCH DIRECTION FROM THE
; GRADIENT FOR THE NON-LINEAR CONJUGATE-GRADIENT CODE.
; IT REPRESENTS A TWO-STEP SELF-SCALED BFGS FORMULA.
;
pro tnmin_msolve, g, y, n, upd1, yksk, gsk, yrsr, lreset, first, $
  hyr, hyk, ykhyk, yrhyr
  compile_opt idl2
  ; ; Rename common variables as they appear in MSLV
  ; common tnmin_work, lsk, lyk, ldiagb, lsr, lyr
  common tnmin_work, sk, yk, diagb, sr, yr
  ; I    I      I     I    I

  ; ; From MSLV
  if keyword_set(upd1) then begin
    y = g / diagb
    RETURN
  endif

  ONE = g[0] * 0 + 1.
  gsk = total(g * sk)

  if keyword_set(lreset) then begin
    ;
    ; COMPUTE GH AND HY WHERE H IS THE INVERSE OF THE DIAGONALS
    ;
    HG = g / diagb
    if keyword_set(first) then begin
      hyk = yk / diagb
      ykhyk = total(yk * hyk)
    endif
    GHYK = total(g * hyk)
    tnmin_ssbfgs, n, ONE, sk, yk, HG, hyk, yksk, ykhyk, gsk, GHYK, y
    RETURN
  endif

  ;
  ; COMPUTE HG AND HY WHERE H IS THE INVERSE OF THE DIAGONALS
  ;
  HG = g / diagb
  if keyword_set(first) then begin
    hyk = yk / diagb
    hyr = yr / diagb
    YKSR = total(yk * sr)
    YKHYR = total(yk * hyr)
  endif
  GSR = total(g * sr)
  GHYR = total(g * hyr)
  if keyword_set(first) then begin
    yrhyr = total(yr * hyr)
  endif

  tnmin_ssbfgs, n, ONE, sr, yr, HG, hyr, yrsr, yrhyr, GSR, GHYR, HG
  if keyword_set(first) then begin
    tnmin_ssbfgs, n, ONE, sr, yr, hyk, hyr, yrsr, yrhyr, YKSR, YKHYR, hyk
  endif
  ykhyk = total(hyk * yk)
  GHYK = total(hyk * g)
  tnmin_ssbfgs, n, ONE, sk, yk, HG, hyk, yksk, ykhyk, gsk, GHYK, y

  RETURN
end

;
; THIS ROUTINE COMPUTES THE PRODUCT OF THE MATRIX G TIMES THE VECTOR
; V AND STORES THE RESULT IN THE VECTOR GV (FINITE-DIFFERENCE VERSION)
;
pro tnmin_gtims, v, gv, n, x, g, fcn, first, delta, accrcy, xnorm, $
  xnew
  compile_opt idl2

  if keyword_set(first) then begin
    ; ; Extra factor of ten is to avoid clashing with the finite
    ; ; difference scheme which computes the derivatives
    delta = sqrt(100 * accrcy) * (1. + xnorm) ; ; XXX diff than TN.F
    ; DELTA = SQRT(ACCRCY)*(1.+XNORM)
    first = 0
  endif
  DINV = 1. / delta

  F = tnmin_call(fcn, x + delta * v, gv, fullparam_ = xnew)
  gv = (gv - g) * DINV
  return
end

;
; UPDATE THE PRECONDITIOING MATRIX BASED ON A DIAGONAL VERSION
; OF THE BFGS QUASI-NEWTON UPDATE.
;
pro tnmin_ndia3, n, e, v, gv, r, vgv
  compile_opt idl2
  VR = total(v * r)
  e = e - r * r / VR + gv * gv / vgv
  wh = where(e le 1d-6, ct)
  if ct gt 0 then e[wh] = 1
  return
end

pro tnmin_fix, whlpeg, whupeg, z
  compile_opt idl2
  if whlpeg[0] ne -1 then z[whlpeg] = 0
  if whupeg[0] ne -1 then z[whupeg] = 0
end

;
; THIS ROUTINE PERFORMS A PRECONDITIONED CONJUGATE-GRADIENT
; ITERATION IN ORDER TO SOLVE THE NEWTON EQUATIONS FOR A SEARCH
; DIRECTION FOR A TRUNCATED-NEWTON ALGORITHM.  WHEN THE VALUE OF THE
; QUADRATIC MODEL IS SUFFICIENTLY REDUCED,
; THE ITERATION IS TERMINATED.
;
; PARAMETERS
;
; ZSOL        - COMPUTED SEARCH DIRECTION
; G           - CURRENT GRADIENT
; GV,GZ1,V    - SCRATCH VECTORS
; R           - RESIDUAL
; DIAGB,EMAT  - DIAGONAL PRECONDITONING MATRIX
; NITER       - NONLINEAR ITERATION #
; FEVAL       - VALUE OF QUADRATIC FUNCTION
pro tnmin_modlnp, zsol, gv, r, v, diagb, emat, $
  x, g, zk, n, niter, maxit, nmodif, nlincg, $
  upd1, yksk, gsk, yrsr, lreset, fcn, whlpeg, whupeg, $
  accrcy, gtp, gnorm, xnorm, xnew
  compile_opt idl2

  ;
  ; GENERAL INITIALIZATION
  ;
  zero = x[0] * 0.
  one = zero + 1
  if (maxit eq 0) then RETURN
  FIRST = 1
  RHSNRM = gnorm
  TOL = zero + 1.e-12
  QOLD = zero

  ;
  ; INITIALIZATION FOR PRECONDITIONED CONJUGATE-GRADIENT ALGORITHM
  ;
  tnmin_initpc, diagb, emat, n, upd1, yksk, gsk, yrsr, lreset

  r = -g
  v = g * 0.
  zsol = v

  ;
  ; ************************************************************
  ; MAIN ITERATION
  ; ************************************************************
  ;
  for K = 1l, maxit do begin
    nlincg = nlincg + 1

    ;
    ; CG ITERATION TO SOLVE SYSTEM OF EQUATIONS
    ;
    tnmin_fix, whlpeg, whupeg, r
    tnmin_msolve, r, zk, n, upd1, yksk, gsk, yrsr, lreset, FIRST, $
      HYR, HYK, YKHYK, YRHYR
    tnmin_fix, whlpeg, whupeg, zk
    RZ = total(r * zk)
    if (RZ / RHSNRM lt TOL) then goto, modlnp_80
    if (K eq 1) then BETA = zero $
    else BETA = RZ / RZOLD
    v = zk + BETA * v
    tnmin_fix, whlpeg, whupeg, v
    tnmin_gtims, v, gv, n, x, g, fcn, FIRST, DELTA, accrcy, xnorm, xnew
    tnmin_fix, whlpeg, whupeg, gv
    VGV = total(v * gv)
    if (VGV / RHSNRM lt TOL) then goto, modlnp_50
    tnmin_ndia3, n, emat, v, gv, r, VGV
    ;
    ; COMPUTE LINEAR STEP LENGTH
    ;
    ALPHA = RZ / VGV
    ;
    ; COMPUTE CURRENT SOLUTION AND RELATED VECTORS
    ;
    zsol = zsol + ALPHA * v
    r = r - ALPHA * gv
    ;
    ; TEST FOR CONVERGENCE
    ;
    gtp = total(zsol * g)
    PR = total(r * zsol)
    QNEW = 5.e-1 * (gtp + PR)
    QTEST = K * (1.e0 - QOLD / QNEW)
    if (QTEST lt 0.d0) then goto, modlnp_70
    QOLD = QNEW
    if (QTEST le 5.d-1) then goto, modlnp_70
    ;
    ; PERFORM CAUTIONARY TEST
    ;
    if (gtp gt 0) then goto, modlnp_40
    RZOLD = RZ
  endfor

  ;
  ; TERMINATE ALGORITHM
  ;
  K = K - 1
  goto, modlnp_70

  modlnp_40:
  zsol = zsol - ALPHA * v
  gtp = total(zsol * g)
  goto, modlnp_90

  modlnp_50:
  ; ; printed output
  modlnp_60:
  if (K gt 1) then goto, modlnp_70
  tnmin_msolve, g, zsol, n, upd1, yksk, gsk, yrsr, lreset, FIRST, $
    HYR, HYK, YKHYK, YRHYR
  zsol = -zsol
  tnmin_fix, whlpeg, whupeg, zsol
  gtp = total(zsol * g)
  modlnp_70:
  goto, modlnp_90
  modlnp_80:
  if (K gt 1) then goto, modlnp_70
  zsol = -g
  tnmin_fix, whlpeg, whupeg, zsol
  gtp = total(zsol * g)
  goto, modlnp_70

  ;
  ; STORE (OR RESTORE) DIAGONAL PRECONDITIONING
  ;
  modlnp_90:
  diagb = emat
  return
end

function tnmin_step1, fnew, fm, gtp, smax, epsmch
  compile_opt idl2

  ; ********************************************************
  ; STEP1 RETURNS THE LENGTH OF THE INITIAL STEP TO BE TAKEN ALONG THE
  ; VECTOR P IN THE NEXT LINEAR SEARCH.
  ; ********************************************************

  D = abs(fnew - fm)
  ALPHA = fnew[0] * 0 + 1.
  if (2.d0 * D le (-gtp) and D ge epsmch) then $
    ALPHA = -2. * D / gtp
  if (ALPHA ge smax) then ALPHA = smax
  return, ALPHA
end

;
; ************************************************************
; GETPTC, AN ALGORITHM FOR FINDING A STEPLENGTH, CALLED REPEATEDLY BY
; ROUTINES WHICH REQUIRE A STEP LENGTH TO BE COMPUTED USING CUBIC
; INTERPOLATION. THE PARAMETERS CONTAIN INFORMATION ABOUT THE INTERVAL
; IN WHICH A LOWER POINT IS TO BE FOUND AND FROM THIS GETPTC COMPUTES A
; POINT AT WHICH THE FUNCTION CAN BE EVALUATED BY THE CALLING PROGRAM.
; THE VALUE OF THE INTEGER PARAMETERS IENTRY DETERMINES THE PATH TAKEN
; THROUGH THE CODE.
; ************************************************************
pro tnmin_getptc, big, small, rtsmll, reltol, abstol, tnytol, $
  fpresn, eta, rmu, xbnd, u, fu, gu, xmin, fmin, gmin, $
  xw, fw, gw, a, b, oldf, b1, scxbnd, e, step, factor, $
  braktd, gtest1, gtest2, tol, ientry, itest
  compile_opt idl2

  ; ; This chicanery is so that we get the data types right
  ZERO = fu[0] * 0.
  ; a1 = zero & scale = zero & chordm = zero
  ; chordu = zero & d1 = zero & d2 = zero
  ; denom = zero
  POINT1 = ZERO + 0.1
  HALF = ZERO + 0.5
  ONE = ZERO + 1
  THREE = ZERO + 3
  FIVE = ZERO + 5
  ELEVEN = ZERO + 11

  if ientry eq 1 then begin ; ; else clause = 20 (OK)
    ;
    ; IENTRY=1
    ; CHECK INPUT PARAMETERS
    ;
    ; ; GETPTC_10:
    itest = 2
    if (u le ZERO or xbnd le tnytol or gu gt ZERO) then RETURN
    itest = 1
    if (xbnd lt abstol) then abstol = xbnd
    tol = abstol
    TWOTOL = tol + tol
    ;
    ; A AND B DEFINE THE INTERVAL OF UNCERTAINTY, X AND XW ARE POINTS
    ; WITH LOWEST AND SECOND LOWEST FUNCTION VALUES SO FAR OBTAINED.
    ; INITIALIZE A,SMIN,XW AT ORIGIN AND CORRESPONDING VALUES OF
    ; FUNCTION AND PROJECTION OF THE GRADIENT ALONG DIRECTION OF SEARCH
    ; AT VALUES FOR LATEST ESTIMATE AT MINIMUM.
    ;
    a = ZERO
    xw = ZERO
    xmin = ZERO
    oldf = fu
    fmin = fu
    fw = fu
    gw = gu
    gmin = gu
    step = u
    factor = FIVE
    ;
    ; THE MINIMUM HAS NOT YET BEEN BRACKETED.
    ;
    braktd = 0
    ;
    ; SET UP XBND AS A BOUND ON THE STEP TO BE TAKEN. (XBND IS NOT COMPUTED
    ; EXPLICITLY BUT SCXBND IS ITS SCALED VALUE.)  SET THE UPPER BOUND
    ; ON THE INTERVAL OF UNCERTAINTY INITIALLY TO XBND + TOL(XBND).
    ;
    scxbnd = xbnd
    b = scxbnd + reltol * abs(scxbnd) + abstol
    e = b + b
    b1 = b
    ;
    ; COMPUTE THE CONSTANTS REQUIRED FOR THE TWO CONVERGENCE CRITERIA.
    ;
    gtest1 = -rmu * gu
    gtest2 = -eta * gu
    ;
    ; SET IENTRY TO INDICATE THAT THIS IS THE FIRST ITERATION
    ;
    ientry = 2
    goto, getptc_210
  endif

  ;
  ; IENTRY = 2
  ;
  ; UPDATE A,B,XW, AND XMIN
  ;
  ; ; GETPTC_20:
  if (fu gt fmin) then goto, getptc_60
  ;
  ; IF FUNCTION VALUE NOT INCREASED, NEW POINT BECOMES NEXT
  ; ORIGIN AND OTHER POINTS ARE SCALED ACCORDINGLY.
  ;
  CHORDU = oldf - (xmin + u) * gtest1
  if not (fu le CHORDU) then begin
    ;
    ; THE NEW FUNCTION VALUE DOES NOT SATISFY THE SUFFICIENT DECREASE
    ; CRITERION. PREPARE TO MOVE THE UPPER BOUND TO THIS POINT AND
    ; FORCE THE INTERPOLATION SCHEME TO EITHER BISECT THE INTERVAL OF
    ; UNCERTAINTY OR TAKE THE LINEAR INTERPOLATION STEP WHICH ESTIMATES
    ; THE ROOT OF F(ALPHA)=CHORD(ALPHA).
    ;
    CHORDM = oldf - xmin * gtest1
    gu = -gmin
    DENOM = CHORDM - fmin
    if (abs(DENOM) lt 1.d-15) then begin
      DENOM = ZERO + 1.e-15
      if (CHORDM - fmin lt 0.d0) then DENOM = -DENOM
    endif
    if (xmin ne ZERO) then gu = gmin * (CHORDU - fu) / DENOM
    fu = (HALF * u * (gmin + gu) + fmin) > fmin
    ;
    ; IF FUNCTION VALUE INCREASED, ORIGIN REMAINS UNCHANGED
    ; BUT NEW POINT MAY NOW QUALIFY AS W.
    ;
    getptc_60:
    if (u ge ZERO) then begin
      b = u
      braktd = 1
    endif else begin
      a = u
    endelse
    xw = u
    fw = fu
    gw = gu
  endif else begin
    ; ; GETPTC_30:
    fw = fmin
    fmin = fu
    gw = gmin
    gmin = gu
    xmin = xmin + u
    a = a - u
    b = b - u
    xw = -u
    scxbnd = scxbnd - u
    if (gu gt ZERO) then begin
      b = ZERO
      braktd = 1
    endif else begin
      a = ZERO
    endelse
    tol = abs(xmin) * reltol + abstol
  endelse

  TWOTOL = tol + tol
  XMIDPT = HALF * (a + b)

  ;
  ; CHECK TERMINATION CRITERIA
  ;
  CONVRG = abs(XMIDPT) le TWOTOL - HALF * (b - a) or $
    abs(gmin) le gtest2 and fmin lt oldf and $
    (abs(xmin - xbnd) gt tol or not braktd)
  if CONVRG then begin
    itest = 0
    if (xmin ne ZERO) then RETURN
    ;
    ; IF THE FUNCTION HAS NOT BEEN REDUCED, CHECK TO SEE THAT THE RELATIVE
    ; CHANGE IN F(X) IS CONSISTENT WITH THE ESTIMATE OF THE DELTA-
    ; UNIMODALITY CONSTANT, TOL.  IF THE CHANGE IN F(X) IS LARGER THAN
    ; EXPECTED, REDUCE THE VALUE OF TOL.
    ;
    itest = 3
    if (abs(oldf - fw) le fpresn * (ONE + abs(oldf))) then RETURN
    tol = POINT1 * tol
    if (tol lt tnytol) then RETURN
    reltol = POINT1 * reltol
    abstol = POINT1 * abstol
    TWOTOL = POINT1 * TWOTOL
  endif

  ;
  ; CONTINUE WITH THE COMPUTATION OF A TRIAL STEP LENGTH
  ;
  ; ; GETPTC_100:
  R = ZERO
  Q = ZERO
  S = ZERO
  if (abs(e) gt tol) then begin
    ;
    ; FIT CUBIC THROUGH XMIN AND XW
    ;
    R = THREE * (fmin - fw) / xw + gmin + gw
    ABSR = abs(R)
    Q = ABSR
    if (gw eq ZERO or gmin eq ZERO) eq 0 then begin ; ; else clause = 140 (OK)
      ;
      ; COMPUTE THE SQUARE ROOT OF (R*R - GMIN*GW) IN A WAY
      ; WHICH AVOIDS UNDERFLOW AND OVERFLOW.
      ;
      ABGW = abs(gw)
      ABGMIN = abs(gmin)
      S = sqrt(ABGMIN) * sqrt(ABGW)
      if ((gw / ABGW) * gmin le ZERO) then begin
        ;
        ; COMPUTE THE SQUARE ROOT OF R*R + S*S.
        ;
        SUMSQ = ONE
        P = ZERO
        if (ABSR lt S) then begin ; ; else clause = 110 (OK)
          ;
          ; THERE IS A POSSIBILITY OF OVERFLOW.
          ;
          if (S gt rtsmll) then P = S * rtsmll
          if (ABSR ge P) then SUMSQ = ONE + (ABSR / S) ^ 2
          SCALE = S
        endif else begin ; ; flow to 120 (OK)
          ;
          ; THERE IS A POSSIBILITY OF UNDERFLOW.
          ;
          ; ; GETPTC_110:
          if (ABSR gt rtsmll) then P = ABSR * rtsmll
          if (S ge P) then SUMSQ = ONE + (S / ABSR) ^ 2
          SCALE = ABSR
        endelse ; ; flow to 120 (OK)
        ; ; GETPTC_120:
        SUMSQ = sqrt(SUMSQ)
        Q = big
        if (SCALE lt big / SUMSQ) then Q = SCALE * SUMSQ
      endif else begin ; ; flow to 140
        ;
        ; COMPUTE THE SQUARE ROOT OF R*R - S*S
        ;
        ; ; GETPTC_130:
        Q = sqrt(abs(R + S)) * sqrt(abs(R - S))
        if (R ge S or R le (-S)) eq 0 then begin
          R = ZERO
          Q = ZERO
          goto, getptc_150
        endif
      endelse
    endif
    ;
    ; COMPUTE THE MINIMUM OF FITTED CUBIC
    ;
    ; ; GETPTC_140:
    if (xw lt ZERO) then Q = -Q
    S = xw * (gmin - R - Q)
    Q = gw - gmin + Q + Q
    if (Q gt ZERO) then S = -S
    if (Q le ZERO) then Q = -Q
    R = e
    if (b1 ne step or braktd) then e = step
  endif

  ;
  ; CONSTRUCT AN ARTIFICIAL BOUND ON THE ESTIMATED STEPLENGTH
  ;
  getptc_150:
  A1 = a
  b1 = b
  step = XMIDPT
  if (braktd) eq 0 then begin ; ; else flow to 160 (OK)
    step = -factor * xw
    if (step gt scxbnd) then step = scxbnd
    if (step ne scxbnd) then factor = FIVE * factor
    ; ; flow to 170 (OK)
  endif else begin
    ;
    ; IF THE MINIMUM IS BRACKETED BY 0 AND XW THE STEP MUST LIE
    ; WITHIN (A,B).
    ;
    ; ; GETPTC_160:
    if (a ne ZERO or xw ge ZERO) and (b ne ZERO or xw le ZERO) then $
      goto, getptc_180
    ;
    ; IF THE MINIMUM IS NOT BRACKETED BY 0 AND XW THE STEP MUST LIE
    ; WITHIN (A1,B1).
    ;
    D1 = xw
    D2 = a
    if (a eq ZERO) then D2 = b
    ; THIS LINE MIGHT BE
    ; IF (A EQ ZERO) THEN D2 = E
    u = -D1 / D2
    step = FIVE * D2 * (POINT1 + ONE / u) / ELEVEN
    if (u lt ONE) then step = HALF * D2 * sqrt(u)
  endelse
  ; ; GETPTC_170:
  if (step le ZERO) then A1 = step
  if (step gt ZERO) then b1 = step
  ;
  ; REJECT THE STEP OBTAINED BY INTERPOLATION IF IT LIES OUTSIDE THE
  ; REQUIRED INTERVAL OR IT IS GREATER THAN HALF THE STEP OBTAINED
  ; DURING THE LAST-BUT-ONE ITERATION.
  ;
  getptc_180:
  if not (abs(S) le abs(HALF * Q * R) or S le Q * A1 or S ge Q * b1) then begin
    ; ; else clause = 200 (OK)
    ;
    ; A CUBIC INTERPOLATION STEP
    ;
    step = S / Q
    ;
    ; THE FUNCTION MUST NOT BE EVALUTATED TOO CLOSE TO A OR B.
    ;
    if not (step - a ge TWOTOL and b - step ge TWOTOL) then begin
      ; ; else clause = 210 (OK)
      if (XMIDPT le ZERO) then step = -tol else step = tol
    endif ; ; flow to 210 (OK)
  endif else begin
    ; ; GETPTC_200:
    e = b - a
  endelse

  ;
  ; IF THE STEP IS TOO LARGE, REPLACE BY THE SCALED BOUND (SO AS TO
  ; COMPUTE THE NEW POINT ON THE BOUNDARY).
  ;
  getptc_210:
  if (step ge scxbnd) then begin ; ; else clause = 220 (OK)
    step = scxbnd
    ;
    ; MOVE SXBD TO THE LEFT SO THAT SBND + TOL(XBND) = XBND.
    ;
    scxbnd = scxbnd - (reltol * abs(xbnd) + abstol) / (ONE + reltol)
  endif
  ; ; GETPTC_220:
  u = step
  if (abs(step) lt tol and step lt ZERO) then u = -tol
  if (abs(step) lt tol and step ge ZERO) then u = tol
  itest = 1
  RETURN
end

;
; LINE SEARCH ALGORITHMS OF GILL AND MURRAY
;
pro tnmin_linder, n, fcn, small, epsmch, reltol, abstol, $
  tnytol, eta, sftbnd, xbnd, p, gtp, x, f, alpha, g, $
  iflag, xnew
  compile_opt idl2

  zero = f[0] * 0.
  one = zero + 1.

  LSPRNT = 0l
  NPRNT = 10000l
  RTSMLL = sqrt(small)
  BIG = 1. / small
  ITCNT = 0l

  ;
  ; SET THE ESTIMATED RELATIVE PRECISION IN F(X).
  ;
  FPRESN = 10. * epsmch
  U = alpha
  FU = f
  FMIN = f
  GU = gtp
  RMU = zero + 1e-4

  ;
  ; FIRST ENTRY SETS UP THE INITIAL INTERVAL OF UNCERTAINTY.
  ;
  IENTRY = 1l

  linder_10:
  ;
  ; TEST FOR TOO MANY ITERATIONS
  ;
  ITCNT = ITCNT + 1
  if (ITCNT gt 30) then begin
    ; ; deviation from Nash: allow optimization to continue in outer
    ; ; loop even if we fail to converge, if IFLAG EQ 0.  A value of
    ; ; 1 indicates failure.  I believe that I tried IFLAG=0 once and
    ; ; there was some problem, but I forget what it was.
    iflag = 1
    f = FMIN
    alpha = XMIN
    x = x + alpha * p
    RETURN
  endif

  iflag = 0
  tnmin_getptc, BIG, small, RTSMLL, reltol, abstol, tnytol, $
    FPRESN, eta, RMU, xbnd, U, FU, GU, XMIN, FMIN, GMIN, $
    XW, FW, GW, A, B, OLDF, B1, SCXBND, E, STEP, FACTOR, $
    BRAKTD, GTEST1, GTEST2, TOL, IENTRY, ITEST

  ;
  ; IF ITEST=1, THE ALGORITHM REQUIRES THE FUNCTION VALUE TO BE
  ; CALCULATED.
  ;
  if (ITEST eq 1) then begin
    UALPHA = XMIN + U
    FU = tnmin_call(fcn, x + UALPHA * p, LG, fullparam_ = xnew)
    GU = total(LG * p)
    ;
    ; THE GRADIENT VECTOR CORRESPONDING TO THE BEST POINT IS
    ; OVERWRITTEN IF FU IS LESS THAN FMIN AND FU IS SUFFICIENTLY
    ; LOWER THAN F AT THE ORIGIN.
    ;
    if (FU le FMIN and FU le OLDF - UALPHA * GTEST1) then $
      g = LG
    ; print, 'fu = ', fu
    goto, linder_10
  endif
  ;
  ; IF ITEST=2 OR 3 A LOWER POINT COULD NOT BE FOUND
  ;
  iflag = 1
  if (ITEST ne 0) then RETURN

  ;
  ; IF ITEST=0 A SUCCESSFUL SEARCH HAS BEEN MADE
  ;
  ; print, 'itcnt = ', itcnt
  iflag = 0
  f = FMIN
  alpha = XMIN
  x = x + alpha * p
  RETURN
end

pro tnmin_defiter, fcn, x, iter, fnorm, fmt = fmt, functargs = fcnargs, $
  quiet = quiet, deriv = df, dprint = dprint, pfixed = pfixed, $
  maximize = maximize, _extra = iterargs
  compile_opt idl2

  if keyword_set(quiet) then return
  if n_params() eq 3 then begin
    fnorm = tnmin_call(fcn, x, df)
  endif

  if keyword_set(maximize) then f = -fnorm else f = fnorm
  print, iter, f, format = '("Iter ",I6,"   FUNCTION = ",G20.8)'
  if n_elements(fmt) gt 0 then begin
    print, x, format = fmt
  endif else begin
    n = n_elements(x)
    ii = lindgen(n)
    p = '     P(' + strtrim(ii, 2) + ') = ' + string(x, format = '(G)')
    if keyword_set(dprint) then begin
      p1 = strarr(n)
      wh = where(pfixed eq 0, ct)
      if ct gt 0 and n_elements(df) ge ct then begin
        if keyword_set(maximize) then df1 = -df else df1 = df
        p1[wh] = string(df1, format = '(G)')
      endif
      wh = where(pfixed eq 1, ct)
      if ct gt 0 then $
        p1[wh] = '          (FIXED)'
      p = p + '  :  DF/DP(' + strtrim(ii, 2) + ') = ' + p1
    endif
    print, p, format = '(A)'
  endelse

  return
end

; SUBROUTINE TNBC (IERROR, N, X, F, G, W, LW, SFUN, LOW, UP, IPIVOT)
; IMPLICIT          DOUBLE PRECISION (A-H,O-Z)
; INTEGER           IERROR, N, LW, IPIVOT(N)
; DOUBLE PRECISION  X(N), G(N), F, W(LW), LOW(N), UP(N)
;
; THIS ROUTINE SOLVES THE OPTIMIZATION PROBLEM
;
; MINIMIZE     F(X)
; X
; SUBJECT TO   LOW <= X <= UP
;
; WHERE X IS A VECTOR OF N REAL VARIABLES.  THE METHOD USED IS
; A TRUNCATED-NEWTON ALGORITHM (SEE "NEWTON-TYPE MINIMIZATION VIA
; THE LANCZOS ALGORITHM" BY S.G. NASH (TECHNICAL REPORT 378, MATH.
; THE LANCZOS METHOD" BY S.G. NASH (SIAM J. NUMER. ANAL. 21 (1984),
; PP. 770-778).  THIS ALGORITHM FINDS A LOCAL MINIMUM OF F(X).  IT DOES
; NOT ASSUME THAT THE FUNCTION F IS CONVEX (AND SO CANNOT GUARANTEE A
; GLOBAL SOLUTION), BUT DOES ASSUME THAT THE FUNCTION IS BOUNDED BELOW.
; IT CAN SOLVE PROBLEMS HAVING ANY NUMBER OF VARIABLES, BUT IT IS
; ESPECIALLY USEFUL WHEN THE NUMBER OF VARIABLES (N) IS LARGE.
;
; SUBROUTINE PARAMETERS:
;
; IERROR  - (INTEGER) ERROR CODE
; ( 0 => NORMAL RETURN
; ( 2 => MORE THAN MAXFUN EVALUATIONS
; ( 3 => LINE SEARCH FAILED TO FIND LOWER
; (          POINT (MAY NOT BE SERIOUS)
; (-1 => ERROR IN INPUT PARAMETERS
; N       - (INTEGER) NUMBER OF VARIABLES
; X       - (REAL*8) VECTOR OF LENGTH AT LEAST N; ON INPUT, AN INITIAL
; ESTIMATE OF THE SOLUTION; ON OUTPUT, THE COMPUTED SOLUTION.
; G       - (REAL*8) VECTOR OF LENGTH AT LEAST N; ON OUTPUT, THE FINAL
; VALUE OF THE GRADIENT
; F       - (REAL*8) ON INPUT, A ROUGH ESTIMATE OF THE VALUE OF THE
; OBJECTIVE FUNCTION AT THE SOLUTION; ON OUTPUT, THE VALUE
; OF THE OBJECTIVE FUNCTION AT THE SOLUTION
; W       - (REAL*8) WORK VECTOR OF LENGTH AT LEAST 14*N
; LW      - (INTEGER) THE DECLARED DIMENSION OF W
; SFUN    - A USER-SPECIFIED SUBROUTINE THAT COMPUTES THE FUNCTION
; AND GRADIENT OF THE OBJECTIVE FUNCTION.  IT MUST HAVE
; THE CALLING SEQUENCE
; SUBROUTINE SFUN (N, X, F, G)
; INTEGER           N
; DOUBLE PRECISION  X(N), G(N), F
; LOW, UP - (REAL*8) VECTORS OF LENGTH AT LEAST N CONTAINING
; THE LOWER AND UPPER BOUNDS ON THE VARIABLES.  IF
; THERE ARE NO BOUNDS ON A PARTICULAR VARIABLE, SET
; THE BOUNDS TO -1.D38 AND 1.D38, RESPECTIVELY.
; IPIVOT  - (INTEGER) WORK VECTOR OF LENGTH AT LEAST N, USED
; TO RECORD WHICH VARIABLES ARE AT THEIR BOUNDS.
;
; THIS IS AN EASY-TO-USE DRIVER FOR THE MAIN OPTIMIZATION ROUTINE
; LMQNBC.  MORE EXPERIENCED USERS WHO WISH TO CUSTOMIZE PERFORMANCE
; OF THIS ALGORITHM SHOULD CALL LMQBC DIRECTLY.
;
; ----------------------------------------------------------------------
; THIS ROUTINE SETS UP ALL THE PARAMETERS FOR THE TRUNCATED-NEWTON
; ALGORITHM.  THE PARAMETERS ARE:
;
; ETA    - SEVERITY OF THE LINESEARCH
; MAXFUN - MAXIMUM ALLOWABLE NUMBER OF FUNCTION EVALUATIONS
; XTOL   - DESIRED ACCURACY FOR THE SOLUTION X*
; STEPMX - MAXIMUM ALLOWABLE STEP IN THE LINESEARCH
; ACCRCY - ACCURACY OF COMPUTED FUNCTION VALUES
; MSGLVL - CONTROLS QUANTITY OF PRINTED OUTPUT
; 0 = NONE, 1 = ONE LINE PER MAJOR ITERATION.
; MAXIT  - MAXIMUM NUMBER OF INNER ITERATIONS PER STEP
;
function tnmin, fcn, xall, fguess = fguess, functargs = fcnargs, parinfo = parinfo, $
  epsrel = epsrel0, epsabs = epsabs0, fastnorm = fastnorm, $
  nfev = nfev, maxiter = maxiter0, maxnfev = maxfun0, maximize = fmax, $
  errmsg = errmsg, nprint = nprint, status = status, nocatch = nocatch, $
  iterproc = iterproc, iterargs = iterargs, niter = niter, quiet = quiet, $
  autoderivative = autoderiv, iterderiv = iterderiv, bestmin = f
  compile_opt idl2

  if n_elements(nprint) eq 0 then nprint = 1
  if n_elements(iterproc) eq 0 then iterproc = 'TNMIN_DEFITER'
  if n_elements(autoderiv) eq 0 then autoderiv = 0
  if n_elements(msglvl) eq 0 then msglvl = 0
  if n_params() eq 0 then begin
    message, 'USAGE: PARMS = TNMIN(''MYFUNCT'', START_PARAMS, ... )', /info
    return, !values.d_nan
  endif
  iterd = keyword_set(iterderiv)
  maximize = keyword_set(fmax)
  status = 0l
  nfev = 0l
  errmsg = ''
  catch_msg = 'in TNMIN'

  common tnmin_config, tnconfig
  tnconfig = {fastnorm: keyword_set(fastnorm), proc: 0, nfev: 0l, $
    autoderiv: keyword_set(autoderiv), max: maximize}

  ; ; Handle error conditions gracefully
  if not keyword_set(nocatch) then begin
    catch, catcherror
    if catcherror ne 0 then begin
      catch, /cancel
      err_string = '' + !error_state.msg
      message, /cont, 'Error detected while ' + catch_msg + ':'
      message, /cont, err_string
      message, /cont, 'Error condition detected. Returning to MAIN level.'
      if errmsg eq '' then $
        errmsg = 'Error detected while ' + catch_msg + ': ' + err_string
      if status eq 0 then status = -18
      return, !values.d_nan
    endif
  endif

  ; ; Parinfo:
  ; ; --------------- STARTING/CONFIG INFO (passed in to routine, not changed)
  ; ; .value   - starting value for parameter
  ; ; .fixed   - parameter is fixed
  ; ; .limited - a two-element array, if parameter is bounded on
  ; ;            lower/upper side
  ; ; .limits - a two-element array, lower/upper parameter bounds, if
  ; ;           limited vale is set
  ; ; .step   - step size in Jacobian calc, if greater than zero

  catch_msg = 'parsing input parameters'
  ; ; Parameters can either be stored in parinfo, or x.  Parinfo takes
  ; ; precedence if it exists.
  if n_elements(xall) eq 0 and n_elements(parinfo) eq 0 then begin
    errmsg = 'ERROR: must pass parameters in X or PARINFO'
    goto, terminate
  endif

  ; ; Be sure that PARINFO is of the right type
  if n_elements(parinfo) gt 0 then begin
    parinfo_size = size(parinfo)
    if parinfo_size[parinfo_size[0] + 1] ne 8 then begin
      errmsg = 'ERROR: PARINFO must be a structure.'
      goto, terminate
    endif
    if n_elements(xall) gt 0 and n_elements(xall) ne n_elements(parinfo) $
    then begin
      errmsg = 'ERROR: number of elements in PARINFO and X must agree'
      goto, terminate
    endif
  endif

  ; ; If the parameters were not specified at the command line, then
  ; ; extract them from PARINFO
  if n_elements(xall) eq 0 then begin
    tnmin_parinfo, parinfo, tagnames, 'VALUE', xall, status = stx
    if stx eq 0 then begin
      errmsg = 'ERROR: either X or PARINFO(*).VALUE must be supplied.'
      goto, terminate
    endif

    sz = size(xall)
    ; ; Convert to double if parameters are not float or double
    if sz[sz[0] + 1] ne 4 and sz[sz[0] + 1] ne 5 then $
      xall = double(xall)
  endif
  npar = n_elements(xall)
  zero = xall[0] * 0.
  one = zero + 1
  ten = zero + 10
  twothird = (zero + 2) / (zero + 3)
  quarter = zero + 0.25
  half = zero + 0.5

  ; ; Extract machine parameters
  sz = size(xall)
  tp = sz[sz[0] + 1]
  if tp ne 4 and tp ne 5 then begin
    if not keyword_set(quiet) then begin
      message, 'WARNING: input parameters must be at least FLOAT', /info
      message, '         (converting parameters to FLOAT)', /info
    endif
    xall = float(xall)
    sz = size(xall)
  endif
  isdouble = (sz[sz[0] + 1] eq 5)

  common tnmin_machar, machvals
  tnmin_setmachar, double = isdouble
  MCHPR1 = machvals.machep

  ; ; TIED parameters?
  tnmin_parinfo, parinfo, tagnames, 'TIED', ptied, default = '', n = npar
  ptied = strtrim(ptied, 2)
  wh = where(ptied ne '', qanytied)
  qanytied = qanytied gt 0
  tnconfig = create_struct(tnconfig, 'QANYTIED', qanytied, 'PTIED', ptied)

  ; ; FIXED parameters ?
  tnmin_parinfo, parinfo, tagnames, 'FIXED', pfixed, default = 0, n = npar
  pfixed = pfixed eq 1
  pfixed = pfixed or (ptied ne '') ; ; Tied parameters are also effectively fixed

  ; ; Finite differencing step, absolute and relative, and sidedness of derivative
  tnmin_parinfo, parinfo, tagnames, 'STEP', step, default = zero, n = npar
  tnmin_parinfo, parinfo, tagnames, 'RELSTEP', dstep, default = zero, n = npar
  tnmin_parinfo, parinfo, tagnames, 'TNSIDE', dside, default = 2, n = npar

  ; ; Maximum and minimum steps allowed to be taken in one iteration
  tnmin_parinfo, parinfo, tagnames, 'TNMAXSTEP', maxstep, default = zero, n = npar
  tnmin_parinfo, parinfo, tagnames, 'TNMINSTEP', minstep, default = zero, n = npar
  qmin = minstep * 0 ; ; Disable minstep for now
  qmax = maxstep ne 0
  wh = where(qmin and qmax and maxstep lt minstep, ct)
  if ct gt 0 then begin
    errmsg = 'ERROR: TNMINSTEP is greater than TNMAXSTEP'
    goto, terminate
  endif
  wh = where(qmin and qmax, ct)
  qminmax = ct gt 0

  ; ; Finish up the free parameters
  ifree = where(pfixed ne 1, ct)
  if ct eq 0 then begin
    errmsg = 'ERROR: no free parameters'
    goto, terminate
  endif

  ; ; Compose only VARYING parameters
  xnew = xall ; ; xnew is the set of parameters to be returned
  x = xnew[ifree] ; ; x is the set of free parameters

  ; ; LIMITED parameters ?
  tnmin_parinfo, parinfo, tagnames, 'LIMITED', limited, status = st1
  tnmin_parinfo, parinfo, tagnames, 'LIMITS', limits, status = st2
  if st1 eq 1 and st2 eq 1 then begin
    ; ; Error checking on limits in parinfo
    wh = where((limited[0, *] and xall lt limits[0, *]) or $
      (limited[1, *] and xall gt limits[1, *]), ct)
    if ct gt 0 then begin
      errmsg = 'ERROR: parameters are not within PARINFO limits'
      goto, terminate
    endif
    wh = where(limited[0, *] and limited[1, *] and $
      limits[0, *] ge limits[1, *] and pfixed eq 0, ct)
    if ct gt 0 then begin
      errmsg = 'ERROR: PARINFO parameter limits are not consistent'
      goto, terminate
    endif

    ; ; Transfer structure values to local variables
    qulim = limited[1, ifree]
    ulim = limits[1, ifree]
    qllim = limited[0, ifree]
    llim = limits[0, ifree]

    wh = where(qulim or qllim, ct)
    if ct gt 0 then qanylim = 1 else qanylim = 0
  endif else begin
    ; ; Fill in local variables with dummy values
    qulim = lonarr(n_elements(ifree))
    ulim = x * 0.
    qllim = qulim
    llim = x * 0.
    qanylim = 0
  endelse

  tnconfig = create_struct(tnconfig, $
    'PFIXED', pfixed, 'IFREE', ifree, $
    'STEP', step, 'DSTEP', dstep, 'DSIDE', dside, $
    'ULIMITED', qulim, 'ULIMIT', ulim)

  common tnmin_fcnargs, tnfcnargs
  tnfcnargs = 0
  dummy = temporary(tnfcnargs)
  if n_elements(fcnargs) gt 0 then tnfcnargs = fcnargs

  ; ; SET UP CUSTOMIZING PARAMETERS
  ; ; ETA    - SEVERITY OF THE LINESEARCH
  ; ; MAXFUN - MAXIMUM ALLOWABLE NUMBER OF FUNCTION EVALUATIONS
  ; ; XTOL   - DESIRED ACCURACY FOR THE SOLUTION X*
  ; ; STEPMX - MAXIMUM ALLOWABLE STEP IN THE LINESEARCH
  ; ; ACCRCY - ACCURACY OF COMPUTED FUNCTION VALUES
  ; ; MSGLVL - DETERMINES QUANTITY OF PRINTED OUTPUT
  ; ;          0 = NONE, 1 = ONE LINE PER MAJOR ITERATION.
  ; ; MAXIT  - MAXIMUM NUMBER OF INNER ITERATIONS PER STEP

  n = n_elements(x)
  if n_elements(maxit) eq 0 then begin
    maxit = (n / 2) < 50 > 2 ; ; XXX diff than TN.F
  endif

  if n_elements(maxfun0) eq 0 then $
    maxfun = 0l $
  else $
    maxfun = floor(maxfun0[0]) > 1
  ; maxfun = 150L*n
  ; if keyword_set(autoderiv) then $
  ; maxfun = maxfun*(1L + round(total(abs(dside)> 1 < 2)))
  eta = zero + 0.25
  stepmx = zero + 10

  if n_elements(maxiter0) eq 0 then $
    maxiter = 200l $
  else $
    maxiter = floor(maxiter0[0]) > 1

  g = replicate(x[0] * 0., n)

  ; ; call minimizer
  ;
  ; THIS ROUTINE IS A BOUNDS-CONSTRAINED TRUNCATED-NEWTON METHOD.
  ; THE TRUNCATED-NEWTON METHOD IS PRECONDITIONED BY A LIMITED-MEMORY
  ; QUASI-NEWTON METHOD (THIS PRECONDITIONING STRATEGY IS DEVELOPED
  ; IN THIS ROUTINE) WITH A FURTHER DIAGONAL SCALING (SEE ROUTINE NDIA3).
  ; FOR FURTHER DETAILS ON THE PARAMETERS, SEE ROUTINE TNBC.
  ;

  ;
  ; initialize variables
  ;
  common tnmin_work, lsk, lyk, ldiagb, lsr, lyr
  ; I/O  I/O     I/O  I/O  I/O
  lsk = 0
  lyk = 0
  ldiagb = 0
  lsr = 0
  lyr = 0

  zero = x[0] * 0.
  one = zero + 1

  if n_elements(fguess) eq 0 then fguess = one
  if maximize then f = -fguess else f = fguess
  conv = 0
  lreset = 0
  upd1 = 0
  newcon = 0
  gsk = zero
  yksk = zero
  gtp = zero
  gtpnew = zero
  yrsr = zero

  upd1 = 1
  ireset = 0l
  nmodif = 0l
  nlincg = 0l
  fstop = f
  conv = 0
  nm1 = n - 1

  ; ; From CHKUCP
  ;
  ; CHECKS PARAMETERS AND SETS CONSTANTS WHICH ARE COMMON TO BOTH
  ; DERIVATIVE AND NON-DERIVATIVE ALGORITHMS
  ;
  EPSMCH = MCHPR1
  SMALL = EPSMCH * EPSMCH
  TINY = SMALL
  NWHY = -1l
  ;
  ; SET CONSTANTS FOR LATER
  ;
  ; ; Some of these constants have been moved around for clarity (!)
  if n_elements(epsrel0) eq 0 then epsrel = 100 * MCHPR1 $
  else epsrel = epsrel0[0] + 0.
  if n_elements(epsabs0) eq 0 then epsabs = zero $
  else epsabs = abs(epsabs0[0]) + 0.

  ACCRCY = epsrel
  XTOL = sqrt(ACCRCY)

  RTEPS = sqrt(EPSMCH)
  RTOL = XTOL
  if (abs(RTOL) lt ACCRCY) then RTOL = 10. * RTEPS

  FTOL2 = 0
  FTOL1 = RTOL ^ 2 + EPSMCH ; ; For func chg convergence test (U1a)
  if epsabs ne 0 then $
    FTOL2 = epsabs + EPSMCH ; ; For absolute func convergence test (U1b)
  PTOL = RTOL + RTEPS ; ; For parm chg convergence test (U2)
  GTOL1 = ACCRCY ^ twothird ; ; For gradient convergence test (U3, squared)
  GTOL2 = (1d-2 * XTOL) ^ 2 ; ; For gradient convergence test (U4, squared)

  ;
  ; CHECK FOR ERRORS IN THE INPUT PARAMETERS
  ;
  if (eta lt 0.d0 or stepmx lt RTOL) then begin
    errmsg = 'ERROR: input keywords are inconsistent'
    goto, terminate
  endif
  ; ; Check input parameters for errors
  if (n le 0) or (XTOL le 0) or (maxit le 0) then begin
    errmsg = 'ERROR: input keywords are inconsistent'
    goto, terminate
  endif
  NWHY = 0l

  XNORM = tnmin_enorm(x)
  ALPHA = zero
  TEST = zero

  ; From SETUCR
  ;
  ; CHECK INPUT PARAMETERS, COMPUTE THE INITIAL FUNCTION VALUE, SET
  ; CONSTANTS FOR THE SUBSEQUENT MINIMIZATION
  ;
  fm = f
  ;
  ; COMPUTE THE INITIAL FUNCTION VALUE
  ;
  catch_msg = 'calling TNMIN_CALL'
  fnew = tnmin_call(fcn, x, g, fullparam_ = xnew)

  ;
  ; SET CONSTANTS FOR LATER
  ;
  niter = 0l
  OLDF = fnew
  GTG = total(g * g)

  common tnmin_error, tnerr

  if nprint gt 0 and iterproc ne '' then begin
    iflag = 0l
    if (niter - 1) mod nprint eq 0 then begin
      catch_msg = 'calling ' + iterproc
      tnerr = 0
      call_procedure, iterproc, fcn, xnew, niter, fnew, $
        functargs = fcnargs, parinfo = parinfo, quiet = quiet, $
        dprint = iterd, deriv = g, pfixed = pfixed, maximize = maximize, $
        _extra = iterargs
      iflag = tnerr
      if iflag lt 0 then begin
        errmsg = 'WARNING: premature termination by "' + iterproc + '"'
        NWHY = 4l
        goto, cleanup
      endif
    endif
  endif

  fold = fnew
  flast = fnew

  ; From LMQNBC
  ;
  ; TEST THE LAGRANGE MULTIPLIERS TO SEE IF THEY ARE NON-NEGATIVE.
  ; BECAUSE THE CONSTRAINTS ARE ONLY LOWER BOUNDS, THE COMPONENTS
  ; OF THE GRADIENT CORRESPONDING TO THE ACTIVE CONSTRAINTS ARE THE
  ; LAGRANGE MULTIPLIERS.  AFTERWORDS, THE PROJECTED GRADIENT IS FORMED.
  ;
  catch_msg = 'zeroing derivatives of pegged parameters'
  lmask = qllim and (x eq llim) and (g ge 0)
  umask = qulim and (x eq ulim) and (g le 0)
  whlpeg = where(lmask, nlpeg)
  whupeg = where(umask, nupeg)
  tnmin_fix, whlpeg, whupeg, g
  GTG = total(g * g)

  ;
  ; CHECK IF THE INITIAL POINT IS A LOCAL MINIMUM.
  ;
  FTEST = one + abs(fnew)
  if (GTG lt GTOL2 * FTEST * FTEST) then goto, cleanup

  ;
  ; SET INITIAL VALUES TO OTHER PARAMETERS
  ;
  ICYCLE = nm1
  GNORM = sqrt(GTG)
  DIFNEW = zero
  EPSRED = half / ten
  FKEEP = fnew

  ;
  ; SET THE DIAGONAL OF THE APPROXIMATE HESSIAN TO UNITY.
  ;
  ldiagb = replicate(one, n)

  ;
  ; ..................START OF MAIN ITERATIVE LOOP..........
  ;
  ; COMPUTE THE NEW SEARCH DIRECTION
  ;
  catch_msg = 'calling TNMIN_MODLNP'
  tnmin_modlnp, lpk, lgv, lz1, lv, ldiagb, lemat, $
    x, g, lzk, n, niter, maxit, nmodif, nlincg, upd1, yksk, $
    gsk, yrsr, lreset, fcn, whlpeg, whupeg, ACCRCY, gtpnew, GNORM, XNORM, $
    xnew

  iter_loop:
  catch_msg = 'computing step length'
  LOLDG = g
  PNORM = tnmin_enorm(lpk)
  OLDF = fnew
  OLDGTP = gtpnew

  ;
  ; PREPARE TO COMPUTE THE STEP LENGTH
  ;
  PE = PNORM + EPSMCH

  ;
  ; COMPUTE THE ABSOLUTE AND RELATIVE TOLERANCES FOR THE LINEAR SEARCH
  ;
  RELTOL = RTEPS * (XNORM + one) / PE
  ABSTOL = -EPSMCH * FTEST / (OLDGTP - EPSMCH)

  ;
  ; COMPUTE THE SMALLEST ALLOWABLE SPACING BETWEEN POINTS IN
  ; THE LINEAR SEARCH
  ;
  TNYTOL = EPSMCH * (XNORM + one) / PE

  ; ; From STPMAX
  SPE = stepmx / PE
  mmask = (not lmask and not umask)
  wh = where(mmask and (lpk gt 0) and qulim and (ulim - x lt SPE * lpk), ct)
  if ct gt 0 then begin
    SPE = min((ulim[wh] - x[wh]) / lpk[wh])
  endif
  wh = where(mmask and (lpk lt 0) and qllim and (llim - x gt SPE * lpk), ct)
  if ct gt 0 then begin
    SPE = min((llim[wh] - x[wh]) / lpk[wh])
  endif

  ; ; From LMQNBC
  ;
  ; SET THE INITIAL STEP LENGTH.
  ;
  ALPHA = tnmin_step1(fnew, fm, OLDGTP, SPE, EPSMCH)

  ;
  ; PERFORM THE LINEAR SEARCH
  ;
  catch_msg = 'performing linear search'
  tnmin_linder, n, fcn, SMALL, EPSMCH, RELTOL, ABSTOL, TNYTOL, $
    eta, zero, SPE, lpk, OLDGTP, x, fnew, ALPHA, g, NWHY, xnew

  newcon = 0
  if (abs(ALPHA - SPE) gt 1.d1 * EPSMCH) eq 0 then begin
    newcon = 1
    NWHY = 0l

    ; ; From MODZ
    mmask = (not lmask and not umask)
    wh = where(mmask and (lpk lt 0) and qllim $
      and (x - llim le 10 * EPSMCH * (abs(llim) + one)), ct)
    if ct gt 0 then begin
      flast = fnew
      lmask[wh] = 1
      x[wh] = llim[wh]
      whlpeg = where(lmask, nlpeg)
    endif
    wh = where(mmask and (lpk gt 0) and qulim $
      and (ulim - x le 10 * EPSMCH * (abs(ulim) + one)), ct)
    if ct gt 0 then begin
      flast = fnew
      umask[wh] = 1
      x[wh] = ulim[wh]
      whupeg = where(umask, nupeg)
    endif
    xnew[ifree] = x

    ; ; From LMQNBC
    flast = fnew
  endif

  fold = fnew
  niter = niter + 1

  ;
  ; IF REQUIRED, PRINT THE DETAILS OF THIS ITERATION
  ;
  if nprint gt 0 and iterproc ne '' then begin
    iflag = 0l
    xx = xnew
    xx[ifree] = x
    if (niter - 1) mod nprint eq 0 then begin
      catch_msg = 'calling ' + iterproc
      tnerr = 0
      call_procedure, iterproc, fcn, xx, niter, fnew, $
        functargs = fcnargs, parinfo = parinfo, quiet = quiet, $
        dprint = iterd, deriv = g, pfixed = pfixed, maximize = maximize, $
        _extra = iterargs
      iflag = tnerr
      if iflag lt 0 then begin
        errmsg = 'WARNING: premature termination by "' + iterproc + '"'
        NWHY = 4l
        goto, cleanup
      endif
    endif
  endif

  catch_msg = 'testing for convergence'
  if (NWHY lt 0) then begin
    NWHY = -2l
    goto, cleanup
  endif
  if (NWHY ne 0 and NWHY ne 2) then begin
    ; ; THE LINEAR SEARCH HAS FAILED TO FIND A LOWER POINT
    NWHY = 3l
    goto, cleanup
  endif
  if NWHY gt 1 then begin
    fnew = tnmin_call(fcn, x, g, fullparam_ = xnew)
  endif
  wh = where(finite(x) eq 0, ct)
  if ct gt 0 or finite(fnew) eq 0 then begin
    NWHY = -3l
    goto, cleanup
  endif

  ;
  ; TERMINATE IF MORE THAN MAXFUN EVALUATIONS HAVE BEEN MADE
  ;
  NWHY = 2l
  if maxfun gt 0 and tnconfig.nfev gt maxfun then goto, cleanup
  if niter gt maxiter then goto, cleanup
  NWHY = 0l

  ;
  ; SET UP PARAMETERS USED IN CONVERGENCE AND RESETTING TESTS
  ;
  DIFOLD = DIFNEW
  DIFNEW = OLDF - fnew

  ;
  ; IF THIS IS THE FIRST ITERATION OF A NEW CYCLE, COMPUTE THE
  ; PERCENTAGE REDUCTION FACTOR FOR THE RESETTING TEST.
  ;
  if (ICYCLE eq 1) then begin
    if (DIFNEW gt 2.d0 * DIFOLD) then EPSRED = EPSRED + EPSRED
    if (DIFNEW lt 5.0d-1 * DIFOLD) then EPSRED = half * EPSRED
  endif
  lgv = g
  tnmin_fix, whlpeg, whupeg, lgv
  GTG = total(lgv * lgv)
  GNORM = sqrt(GTG)
  FTEST = one + abs(fnew)
  XNORM = tnmin_enorm(x)

  ; ; From CNVTST
  LTEST = (flast - fnew) le (-5.d-1 * gtpnew)
  wh = where((lmask and g lt 0) or (umask and g gt 0), ct)
  if ct gt 0 then begin
    conv = 0
    if not LTEST then begin
      mx = max(abs(g[wh]), wh2)
      lmask[wh[wh2]] = 0
      umask[wh[wh2]] = 0
      flast = fnew
      goto, cnvtst_done
    endif
  endif
  ; ; Gill Murray and Wright tests are listed to the right.
  ; ; Modifications due to absolute function value test are done here.

  fconv = abs(DIFNEW) lt FTOL1 * FTEST ; ; U1a
  if FTOL2 eq 0 then begin
    pconv = ALPHA * PNORM lt PTOL * (one + XNORM) ; ; U2
    gconv = GTG lt GTOL1 * FTEST * FTEST ; ; U3
  endif else begin
    ; ; Absolute tolerance implies a loser constraint on parameters
    fconv = fconv or (abs(DIFNEW) lt FTOL2) ; ; U1b
    acc2 = (FTOL2 / FTEST + EPSMCH)
    pconv = ALPHA * PNORM lt sqrt(acc2) * (one + XNORM) ; ; U2
    gconv = GTG lt (acc2 ^ twothird) * FTEST * FTEST ; ; U3
  endelse
  if ((pconv and fconv and gconv) $ ; ; U1 + U2 + U3
    or (GTG lt GTOL2 * FTEST * FTEST)) then begin ; ; U4
    conv = 1
  endif else begin
    ; ; Convergence failed
    conv = 0
  endelse

  ;
  ; FOR DETAILS, SEE GILL, MURRAY, AND WRIGHT (1981, P. 308) AND
  ; FLETCHER (1981, P. 116).  THE MULTIPLIER TESTS (HERE, TESTING
  ; THE SIGN OF THE COMPONENTS OF THE GRADIENT) MAY STILL NEED TO
  ; MODIFIED TO INCORPORATE TOLERANCES FOR ZERO.
  ;

  cnvtst_done:
  ; ; From LMQNBC

  if (conv) then goto, cleanup
  tnmin_fix, whlpeg, whupeg, g

  ;
  ; COMPUTE THE CHANGE IN THE ITERATES AND THE CORRESPONDING CHANGE
  ; IN THE GRADIENTS
  ;
  if newcon eq 0 then begin
    lyk = g - LOLDG
    lsk = ALPHA * lpk
    ;
    ; SET UP PARAMETERS USED IN UPDATING THE PRECONDITIONING STRATEGY.
    ;
    yksk = total(lyk * lsk)
    lreset = 0

    if (ICYCLE eq nm1 or DIFNEW lt EPSRED * (FKEEP - fnew)) then lreset = 1
    if (lreset eq 0) then begin
      yrsr = total(lyr * lsr)
      if (yrsr le zero) then lreset = 1
    endif
    upd1 = 0
  endif

  ;
  ; COMPUTE THE NEW SEARCH DIRECTION
  ;
  ; ; TNMIN_90:
  catch_msg = 'calling TNMIN_MODLNP'

  tnmin_modlnp, lpk, lgv, lz1, lv, ldiagb, lemat, $
    x, g, lzk, n, niter, maxit, nmodif, nlincg, upd1, yksk, $
    gsk, yrsr, lreset, fcn, whlpeg, whupeg, ACCRCY, gtpnew, GNORM, XNORM, $
    xnew

  if (newcon) then goto, iter_loop
  ; IF (NOT LRESET) OR ICYCLE EQ 1 AND n_elements(LSR) GT 0 THEN BEGIN   ;; For testing
  if (lreset eq 0) then begin
    ;
    ; COMPUTE THE ACCUMULATED STEP AND ITS CORRESPONDING
    ; GRADIENT DIFFERENCE.
    ;
    lsr = lsr + lsk
    lyr = lyr + lyk
    ICYCLE = ICYCLE + 1
    goto, iter_loop
  endif

  ;
  ; RESET
  ;
  ; ; TNMIN_110:
  ireset = ireset + 1
  ;
  ; INITIALIZE THE SUM OF ALL THE CHANGES IN X.
  ;
  lsr = lsk
  lyr = lyk
  FKEEP = fnew
  ICYCLE = 1l
  goto, iter_loop

  ;
  ; ...............END OF MAIN ITERATION.......................
  ;
  cleanup:

  nfev = tnconfig.nfev
  tnfcnargs = 0
  catch, /cancel
  case NWHY of
    -3: begin
      ; ; INDEFINITE VALUE
      status = -16l
      if errmsg eq '' then $
        errmsg = ('ERROR: parameter or function value(s) have become ' + $
          'infinite; check model function for over- ' + $
          'and underflow')
      return, !values.d_nan
    end
    - 2: begin
      ; ; INTERNAL ERROR IN LINE SEARCH
      status = -18l
      if errmsg eq '' then $
        errmsg = 'ERROR: Fatal error during line search'
      return, !values.d_nan
    end
    - 1: begin
      terminate:
      ; ; FATAL TERMINATION
      status = 0l
      if errmsg eq '' then errmsg = 'ERROR: Invalid inputs'
      return, !values.d_nan
    end
    0: begin
      converged:
      status = 1l
    end
    2: begin
      ; ; MAXIMUM NUMBER of FUNC EVALS or ITERATIONS REACHED
      if maxfun gt 0 and nfev gt maxfun then begin
        status = -17l
        if errmsg eq '' then $
          errmsg = ('WARNING: no convergence within maximum ' + $
            'number of function calls')
      endif else begin
        status = 5l
        if errmsg eq '' then $
          errmsg = ('WARNING: no convergence within maximum ' + $
            'number of iterations')
      endelse
      fnew = OLDF
    end
    3: begin
      status = -18l
      if errmsg eq '' then errmsg = 'ERROR: Line search failed to converge'
    end
    4: begin
      ; ; ABNORMAL TERMINATION BY USER ROUTINE
      status = iflag
    end
  endcase

  ; ; Successful return
  f = fnew
  xnew[ifree] = x
  return, xnew
end