pro get_charis_wvlh, header, wavelengths, wavecalfile = wavecalfile, manual = manual, filtname = filtname, help = help
  compile_opt idl2
  ; ,filtname=filtname
  ; v2.

  if (n_params() eq 0 or keyword_set(help)) then begin
    print, 'get_charis_wvlh,header,wavelengths,wavecalfile=wavecalfile,manual=manual'
    print, ''
    print, 'simple program to get the wavelength array for CHARIS data'
    print, ''
    print, '***Keywords***'
    print, ''
    print, '*manual - [values: ''lowres'', ''J'',''H'',''K''] - set this to return the wavelength array for filter you manually set'
    print, '*header - the 0th extension of the CHARIS fits header: determine filter from here'
    print, '*wavecalfile - override the choice of file read in [rarely used]'
    goto, skiptotheend
  endif

  ; will need to change for each user!!!
  wavecalpath = charis_path(pathname = 'wavecalpath')

  if ~keyword_set(manual) then begin
    ; Takes the fits header from Charis and returns the wavelengths
    ; Note: must change wavecal path in line 9

    filtname = strtrim(sxpar(header, 'CAL_BAND', count = filtcount)) ; for right now, using the calibration band to define the mode since FILTNAME gets unpopulated.

    if filtcount eq 0 then read, 'Select the filter to use (lowres,J,H,K)', filtname
  endif else begin
    filtname = manual
  endelse

  if ~keyword_set(wavecalfile) then begin
    ; broadband
    case filtname of
      'J': wavecalfile = 'highres_jwvlh.txt'
      'H': wavecalfile = 'highres_hwvlh.txt'
      'K': wavecalfile = 'highres_kwvlh.txt'
      'Broadband': wavecalfile = 'lowres_wvlh.txt'
      'lowres': wavecalfile = 'lowres_wvlh.txt'
      else: wavecalfile = 'lowres_wvlh.txt'
    endcase
  endif
  ; override
  ; wavelengths=readfits(wavecalpath+wavecalfile,h1)
  readcol, wavecalpath + wavecalfile, wavelengths, /silent

  skiptotheend:
end