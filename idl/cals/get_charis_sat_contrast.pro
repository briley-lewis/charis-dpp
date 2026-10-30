function get_charis_sat_contrast, modulation_size, wavelengths, mjd, errors = errors, manual = manual0, help = help
  compile_opt idl2

  if (n_params() eq 0 or keyword_set(help)) then begin
    print, 'function get_charis_sat_contrast,modulation_size,wavelengths,mjd,errors=errors,manual=manual0'
    print, ''
    print, '***Keywords***'
    print, ''
    print, '*modulation - the modulation amplitude used'
    print, '*wavelengths - wavelength array'
    print, '*mjd - modified Julien day (pulled from fits headers)'
    print, '*errors - return the errors?'
    print, '*manual - manually set the contrast at a central wavelength of 1550 nm'

    goto, skiptotheend
  endif

  ; Given the modulation amplitude, returns the attenuation of the satellite spots/wavelength.  Formally only directly tested for broadband data but probably fine for highres.
  ; ****Based on internal source testing 07/20/2017 and earlier in summer 2017; note: values from Lozi+18 tests give a 25nm contrast ~8% higher.

  ; Decision tree.
  ; 1. You read the mod size from the fits headers
  ; --If the modulation size is 50nm, then you know it precisely
  ; --If the modulation size is 25nm, then you know it precisely
  ; --If the modulation size is something else, then extrapolate, knowing that the precision is good to ~a few percent or so.

  atten_fac = fltarr(n_elements(wavelengths))
  central_wvlh = 1.550
  conv_fac = (central_wvlh / wavelengths) ^ 2.
  if ~keyword_set(manual0) then begin
    ; [BUGS.md #68, #122] the contrast table lives in charis_get_constant (epochs
    ; keyed on MJD); the caller's modulation_size is no longer modified
    modsize = long(round(modulation_size[0]))
    mag25 = charis_get_constant(name = 'satcon25', mjd = mjd)
    mag50 = charis_get_constant(name = 'satcon50', mjd = mjd)
    magerr = charis_get_constant(name = 'satconerr', mjd = mjd)
    case modsize of
      50: atten_fac0 = 10 ^ (-0.4 * mag50)
      25: atten_fac0 = 10 ^ (-0.4 * mag25)
      else: atten_fac0 = (modsize / 25.) ^ 2. * 10 ^ (-0.4 * mag25)
    endcase
    atten_fac = atten_fac0 * conv_fac
    errors = magerr / 1.0857 * atten_fac
  endif else begin
    ; manual attenuation: a scalar (scaled with wavelength) or one value per channel
    if n_elements(manual0) eq n_elements(wavelengths) then atten_fac = manual0 else atten_fac = manual0[0] * conv_fac
    errors = 0.05 / 1.0857 * atten_fac
  endelse
  return, atten_fac

  skiptotheend:
end