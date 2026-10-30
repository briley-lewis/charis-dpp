function wollaston_efficiency, lambda
    compile_opt idl2
    ; linearly interpolates fitted wollaston efficiency from McIntosh+ 2026 


    ; Read the CSV file with wavelength-fit wollaston eff
  csv_path = file_dirname(routine_filepath()) + '/mcintosh_2026_spie_wavelength_fits.csv'
  data = read_csv(csv_path, header = headers)
  headers = strlowcase(headers)

  wavelength_idx = where(headers eq 'wavelength_bin', n_wavelength_idx)
  wollaston_epsilon_idx = where(headers eq 'wollaston_eta', n_wollaston_epsilon_idx)

  wavelengths = data.(wavelength_idx[0])  
  eps = data.(wollaston_epsilon_idx[0])

  ; Interpolate via wavelength
  epsilon_bs = interpol(eps, wavelengths, lambda)

  return, epsilon_bs
end